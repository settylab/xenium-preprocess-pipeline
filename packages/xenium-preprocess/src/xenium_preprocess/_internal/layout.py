"""Run-scoped output layout for xenium-preprocess.

Layout (locked in (internal issue review):

    <output_root>/
    └── <sample_id>/
        └── <sample_id>_<run_id>/
            ├── spatial_adata/
            │   ├── <sample_id>_proseg_raw.h5ad
            │   ├── <sample_id>_xenium_ranger.h5ad
            │   └── <sample_id>_proseg_purified.h5ad     (written by step 4)
            ├── rctd/
            │   ├── <sample_id>_test_object.rds          (this pipeline)
            │   ├── <sample_id>_reference.rds            (step 3)
            │   ├── <sample_id>_reference_post_rules.h5ad (step 3)
            │   └── <sample_id>_rctd_results.rds         (step 4)
            ├── config.yaml                     (merged across steps)
            └── logs/{xenium-preprocess,ref-build,rctd-split}.log

`run_id` precedence (bound in cli._resolve_run_id): `--run-id` >
`$SLURM_JOB_ID` > `YYYYMMDD_HHMMSS` timestamp fallback. Every path this
package writes is derived from `run_dir()`; nothing else knows the
folder shape.
"""
from __future__ import annotations

import os
from pathlib import Path


# Keys index which h5ad this pipeline knows how to place. Each maps to
# the basename inside `spatial_adata/`. Step 3 / step 4 outputs live
# under `rctd/` and are placed by their own pipelines — the constants
# for those are captured in RCTD_BASENAMES below purely for downstream
# stages inside THIS pipeline that need to READ them.
_SPATIAL_ADATA_BASENAMES = {
    "proseg_raw": "{sample_id}_proseg_raw.h5ad",
    "xenium_ranger": "{sample_id}_xenium_ranger.h5ad",
    "proseg_purified": "{sample_id}_proseg_purified.h5ad",
}
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


def spatial_adata_dir(output_root: Path, sample_id: str, run_id: str) -> Path:
    return run_dir(output_root, sample_id, run_id) / "spatial_adata"


def rctd_dir(output_root: Path, sample_id: str, run_id: str) -> Path:
    return run_dir(output_root, sample_id, run_id) / "rctd"


def logs_dir(output_root: Path, sample_id: str, run_id: str) -> Path:
    return run_dir(output_root, sample_id, run_id) / "logs"


def resolved_config_path(output_root: Path, sample_id: str, run_id: str) -> Path:
    return run_dir(output_root, sample_id, run_id) / "config.yaml"


def spatial_adata_path(
    output_root: Path, sample_id: str, run_id: str, key: str
) -> Path:
    """Return the path to a h5ad under `spatial_adata/`. `key` is one
    of `_SPATIAL_ADATA_BASENAMES`."""
    basename = _SPATIAL_ADATA_BASENAMES[key].format(sample_id=sample_id)
    return spatial_adata_dir(output_root, sample_id, run_id) / basename


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
    `pd.read_csv` yields columns backed by `pd.arrays.ArrowStringArray`
    (dtype `pd.StringDtype`). anndata's writer registry has no
    `ArrowStringArray` handler on either the index code path or the
    `write_categorical` code path (which serializes a Categorical's
    `.categories._values` separately). Failures observed:
    - `/obs` index (internal issue review) —
      string index of `obs` from `pd.Index(series, str, True, ...)`.
    - `/obs/fov` categorical categories
      (internal issue review) — `obs['fov']` is
      Categorical (either explicit `.astype('category')` or anndata's
      auto-categorization of string obs columns) and its `.categories`
      Index is backed by ArrowStringArray.
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
    default (`convert_strings_to_categoricals=True`), which converts
    string obs/var columns into Categoricals whose `.categories` are
    `ArrowStringArray`-backed — the exact shape anndata's writer
    registry then can't serialize (internal issue review). To keep our sanitizer effective,
    do the categorization ourselves BEFORE sanitize (on frame copies
    so the caller's DataFrames stay intact); `_sanitize_frame_for_h5ad`
    then rebuilds each affected Categorical with object-backed
    categories, and anndata's re-invocation inside `write_h5ad` is a
    no-op (columns are already Categorical) — the sanitize survives
    to the writer.
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
