"""Run-scoped output layout for ref-build.

Layout (locked in (internal issue review):

    <output_root>/
    └── <sample_id>/
        └── <sample_id>_<run_id>/
            ├── spatial_adata/                        (xenium-preprocess)
            ├── rctd/
            │   ├── <sample_id>_test_object.rds       (xenium-preprocess)
            │   ├── <sample_id>_reference_post_rules.h5ad   (THIS step)
            │   ├── <sample_id>_reference.rds               (THIS step)
            │   └── <sample_id>_rctd_results.rds      (rctd-split)
            ├── config.yaml                  (merged across steps)
            └── logs/{xenium-preprocess,ref-build,rctd-split}.log

Intermediate outputs (`loaded/`, `census/`, `mtx_bundle/`) also land
under `<run_dir>/` while the pipeline runs, and are DROPPED by
`pipeline._drop_intermediate_outputs` after `rctd_reference_build`
succeeds.

`run_id` precedence (bound in cli._resolve_run_id): `--run-id` >
`$SLURM_JOB_ID` > `YYYYMMDD_HHMMSS` timestamp fallback. Every terminal
path this pipeline writes is derived from `run_dir()`; nothing else
knows the folder shape.
"""
from __future__ import annotations

import os
from pathlib import Path


_RCTD_BASENAMES = {
    "test_object": "{sample_id}_test_object.rds",
    "reference": "{sample_id}_reference.rds",
    "reference_post_rules": "{sample_id}_reference_post_rules.h5ad",
    "rctd_results": "{sample_id}_rctd_results.rds",
}


def run_dir(output_root: Path, sample_id: str, run_id: str) -> Path:
    """`<output_root>/<sample_id>/<sample_id>_<run_id>/` — the anchor
    every other path is derived from."""
    return Path(output_root) / sample_id / f"{sample_id}_{run_id}"


def rctd_dir(output_root: Path, sample_id: str, run_id: str) -> Path:
    return run_dir(output_root, sample_id, run_id) / "rctd"


def logs_dir(output_root: Path, sample_id: str, run_id: str) -> Path:
    return run_dir(output_root, sample_id, run_id) / "logs"


def resolved_config_path(output_root: Path, sample_id: str, run_id: str) -> Path:
    return run_dir(output_root, sample_id, run_id) / "config.yaml"


def rctd_path(
    output_root: Path, sample_id: str, run_id: str, key: str
) -> Path:
    """Return the path to a file under `rctd/`. `key` is one of
    `_RCTD_BASENAMES`."""
    basename = _RCTD_BASENAMES[key].format(sample_id=sample_id)
    return rctd_dir(output_root, sample_id, run_id) / basename


def _sanitize_frame_for_h5ad(df):
    """Convert pandas nullable-string extension columns/index to numpy
    object so anndata's h5ad `write_elem` registry can serialize them.

    Pandas 3.x ships `future.infer_string=True` as the default, so
    `pd.read_csv` / `ad.read_h5ad` yield columns and indices backed by
    `pd.arrays.ArrowStringArray` (dtype `pd.StringDtype`). anndata
    0.11's writer registry has no `ArrowStringArray` handler on either
    the index code path or the `write_categorical` code path (which
    serializes a Categorical's `.categories._values` separately).
    Failures observed:
    - `/obs` index (internal issue review) — string index of `obs` from `ad.read_h5ad(...)` or
      after `ad.concat(...)`.
    - `/obs/<col>` categorical categories
      (internal issue review) — obs column is
      Categorical (either explicit `.astype('category')` or anndata's
      auto-categorization of string obs columns) and its `.categories`
      Index is backed by ArrowStringArray.

    This is a verbatim replica of `xenium_preprocess._internal.layout.
    _sanitize_frame_for_h5ad` (internal issue review). Kept in
    sync manually — if either implementation changes, update both. See
    the report on ref-build's arrow-string fix for the rationale on
    replication vs. cross-package import.
    """
    import numpy as np
    import pandas as pd

    def _is_string_ext(dtype):
        return isinstance(dtype, pd.StringDtype)

    def _cat_needs_categories_fix(dtype):
        return (
            isinstance(dtype, pd.CategoricalDtype)
            and _is_string_ext(dtype.categories.dtype)
        )

    def _needs(dtype):
        return _is_string_ext(dtype) or _cat_needs_categories_fix(dtype)

    def _fix_series(s):
        if _is_string_ext(s.dtype):
            return s.astype(object)
        # Categorical with arrow-string categories: rebuild with
        # object-dtype categories, preserving codes and ordered flag so
        # the write path stays on `write_categorical` (which downstream
        # readers rely on for `.cat` semantics + storage efficiency).
        new_cats = pd.Index(
            np.asarray(s.cat.categories, dtype=object),
            dtype=object,
        )
        new_cat = pd.Categorical.from_codes(
            s.cat.codes.to_numpy(),
            categories=new_cats,
            ordered=s.cat.ordered,
        )
        return pd.Series(new_cat, index=s.index, name=s.name)

    cols_to_fix = [c for c in df.columns if _needs(df[c].dtype)]
    index_needs_fix = _is_string_ext(df.index.dtype)
    if not cols_to_fix and not index_needs_fix:
        return df

    df = df.copy()
    for c in cols_to_fix:
        df[c] = _fix_series(df[c])
    if index_needs_fix:
        # `dtype=object` is load-bearing: with `future.infer_string=True`
        # (pandas 3.x default), a bare `pd.Index([...])` re-infers strings
        # to Arrow, defeating the sanitize.
        import pandas as pd  # noqa: F811 — keep name local
        df.index = pd.Index(
            np.asarray(df.index, dtype=object),
            dtype=object,
            name=df.index.name,
        )
    return df


def atomic_write_h5ad(adata, out_path: Path, **write_kwargs) -> None:
    """Write via `<out_path>.tmp` + `os.replace` so a mid-write crash
    leaves the pre-existing file intact.

    Under pandas 3.x (`future.infer_string=True`) `anndata.write_h5ad`
    runs `adata.strings_to_categoricals()` on entry as its documented
    default, converting string obs/var columns into Categoricals whose
    `.categories` are `ArrowStringArray`-backed — the exact shape
    anndata 0.11's writer registry can't serialize
    (internal issue review)
    (internal issue review)). Do the categorization ourselves BEFORE sanitize (on
    frame copies so the caller's DataFrames stay intact);
    `_sanitize_frame_for_h5ad` then rebuilds each affected Categorical
    with object-backed categories, and anndata's re-invocation inside
    `write_h5ad` is a no-op — the sanitize survives to the writer.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")

    old_obs, old_var = adata.obs, adata.var
    # `strings_to_categoricals` mutates the frame in place; work on
    # copies so the caller's obs/var references stay untouched.
    adata.obs = old_obs.copy()
    adata.var = old_var.copy()
    try:
        adata.strings_to_categoricals()
        fixed_obs = _sanitize_frame_for_h5ad(adata.obs)
        fixed_var = _sanitize_frame_for_h5ad(adata.var)
        if fixed_obs is not adata.obs:
            adata.obs = fixed_obs
        if fixed_var is not adata.var:
            adata.var = fixed_var
        adata.write_h5ad(tmp, **write_kwargs)
    finally:
        adata.obs = old_obs
        adata.var = old_var

    os.replace(tmp, out_path)
