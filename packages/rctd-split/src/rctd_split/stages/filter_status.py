"""Stage 5: per-cell filter-status columns for the intermediate unpurified h5ad.

Adds per-cell provenance telling you WHICH filter (RCTD reject vs. SPLIT
purification drop) dropped each cell.

Two outputs from one computation:

  1. Sidecar CSV at
       <run_dir>/intermediate/adata/<sample_id>_filter_status.csv
     with columns: cell_id, passed_rctd, passed_split_purify,
     filtered_by_purification.
  2. In-place augmentation of the intermediate ``<sample_id>_unpurified.h5ad``
     — adds the same three columns to .obs. Only runs when ``write_inplace=True``.

Column semantics (indexed by unpurified.obs_names):

    passed_rctd              = ~unpurified.obs[first_type_col].isna()
    passed_split_purify      = unpurified.obs_names.isin(purified.obs_names)
    filtered_by_purification = passed_rctd & ~passed_split_purify

Idempotent: recomputes and overwrites every run when the sentinel is
absent (or ``force_rerun=True``). No "already up to date" detection inside
the stage; the sidecar CSV is the pipeline-level sentinel.
"""
from __future__ import annotations

import os
from pathlib import Path

from rctd_split._internal.compat import sentinel_exists
from rctd_split._internal.layout import (
    atomic_write_h5ad,
    intermediate_path,
    spatial_adata_path,
)
from rctd_split._internal.logging import log

DEFAULT_INPLACE_COLUMNS = (
    "passed_rctd",
    "passed_split_purify",
    "filtered_by_purification",
)


def _compute_columns(
    unpurified_obs_names,
    unpurified_first_type,
    purified_obs_names,
):
    """Return the three columns as a dict of ``pandas.Series`` indexed
    by ``unpurified_obs_names``."""
    import pandas as pd

    idx = pd.Index(unpurified_obs_names, name="cell_id")
    passed_rctd = pd.Series(~pd.isna(unpurified_first_type), index=idx,
                            name="passed_rctd").astype(bool)
    passed_split_purify = pd.Series(
        idx.isin(pd.Index(purified_obs_names)),
        index=idx, name="passed_split_purify",
    ).astype(bool)
    filtered_by_purification = pd.Series(
        passed_rctd.to_numpy() & ~passed_split_purify.to_numpy(),
        index=idx, name="filtered_by_purification",
    ).astype(bool)
    return {
        "passed_rctd": passed_rctd,
        "passed_split_purify": passed_split_purify,
        "filtered_by_purification": filtered_by_purification,
    }


def _write_sidecar_csv(cols: dict, out_csv: Path) -> None:
    import pandas as pd

    df = pd.DataFrame(
        {
            "passed_rctd": cols["passed_rctd"].to_numpy(),
            "passed_split_purify": cols["passed_split_purify"].to_numpy(),
            "filtered_by_purification": cols["filtered_by_purification"].to_numpy(),
        },
        index=cols["passed_rctd"].index,
    )
    df.index.name = "cell_id"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_csv.with_suffix(out_csv.suffix + ".tmp")
    df.to_csv(tmp, index=True)
    os.replace(tmp, out_csv)


def _augment_h5ad_inplace(
    unpurified_h5ad: Path,
    cols: dict,
    inplace_columns,
    h5ad_compression: str | None,
) -> None:
    """Add the requested filter-status columns to ``.obs`` and atomically
    rewrite. Overwrites the columns if they already exist (idempotent
    re-run)."""
    import anndata as ad

    log(f"[filter_status]   reading {unpurified_h5ad} for in-place augmentation")
    adata = ad.read_h5ad(unpurified_h5ad)

    for name in inplace_columns:
        if name not in cols:
            raise SystemExit(
                f"[filter_status] inplace_columns entry {name!r} is not one of "
                f"the computed columns {list(cols)}"
            )
        series = cols[name].reindex(adata.obs_names)
        if series.isna().any():
            missing = int(series.isna().sum())
            raise SystemExit(
                f"[filter_status] {missing} unpurified cells missing from the "
                f"computed {name!r} series — obs_names mismatch bug"
            )
        adata.obs[name] = series.astype(bool).to_numpy()

    write_kwargs = {"compression": h5ad_compression} if h5ad_compression else {}
    atomic_write_h5ad(adata, unpurified_h5ad, **write_kwargs)
    log(f"[filter_status]   wrote augmented {unpurified_h5ad} "
        f"(+{len(inplace_columns)} .obs columns)")


def run_filter_status(
    sample_id: str,
    run_id: str,
    output_root: Path,
    write_inplace: bool,
    inplace_columns,
    first_type_column: str,
    h5ad_compression: str | None,
    force_rerun: bool,
) -> Path:
    """Compute the three filter-status columns and write the sidecar CSV
    (always) + optionally augment the intermediate unpurified h5ad in-place.

    Returns the sidecar CSV path.
    """
    import anndata as ad

    unpurified_h5ad = intermediate_path(
        output_root, sample_id, run_id, "unpurified_h5ad",
    )
    purified_h5ad = spatial_adata_path(
        output_root, sample_id, run_id, "proseg_purified",
    )
    sidecar_csv = intermediate_path(
        output_root, sample_id, run_id, "filter_status_csv",
    )
    sidecar_csv.parent.mkdir(parents=True, exist_ok=True)

    if sentinel_exists(sidecar_csv, force_rerun):
        if not write_inplace:
            log(f"[filter_status] sentinel exists: {sidecar_csv} — skipping "
                "(pass --force-rerun to re-run).")
            return sidecar_csv
        try:
            adata = ad.read_h5ad(unpurified_h5ad, backed="r")
            present = set(adata.obs.columns)
            adata.file.close()
        except Exception:
            present = set()
        if all(c in present for c in inplace_columns):
            log(f"[filter_status] sentinels exist: {sidecar_csv} + inplace "
                f"columns on {unpurified_h5ad} — skipping (pass --force-rerun "
                "to re-run).")
            return sidecar_csv

    if not unpurified_h5ad.exists():
        raise SystemExit(
            f"[filter_status] unpurified h5ad not found: {unpurified_h5ad}. "
            "Run the mtx_to_h5ad stage first."
        )
    if not purified_h5ad.exists():
        raise SystemExit(
            f"[filter_status] purified h5ad not found: {purified_h5ad}. "
            "Run the mtx_to_h5ad stage first."
        )

    log(f"[filter_status] reading {unpurified_h5ad}")
    unp = ad.read_h5ad(unpurified_h5ad)
    log(f"[filter_status] reading {purified_h5ad}")
    pur = ad.read_h5ad(purified_h5ad)

    if first_type_column not in unp.obs.columns:
        raise SystemExit(
            f"[filter_status] unpurified.obs is missing the RCTD first-type "
            f"column {first_type_column!r}. Present columns: "
            f"{list(unp.obs.columns)}"
        )

    cols = _compute_columns(
        unpurified_obs_names=list(unp.obs_names),
        unpurified_first_type=unp.obs[first_type_column].to_numpy(),
        purified_obs_names=list(pur.obs_names),
    )

    total = len(cols["passed_rctd"])
    n_rctd = int(cols["passed_rctd"].sum())
    n_split = int(cols["passed_split_purify"].sum())
    n_filt = int(cols["filtered_by_purification"].sum())

    def _pct(n: int) -> str:
        return f"{100.0 * n / total:.1f}" if total else "0.0"

    log(
        f"[filter_status] total={total} "
        f"passed_rctd={n_rctd} ({_pct(n_rctd)}%) "
        f"passed_split_purify={n_split} ({_pct(n_split)}%) "
        f"filtered_by_purification={n_filt} ({_pct(n_filt)}%)"
    )

    _write_sidecar_csv(cols, sidecar_csv)
    log(f"[filter_status] wrote sidecar {sidecar_csv}")

    if write_inplace:
        _augment_h5ad_inplace(
            unpurified_h5ad=unpurified_h5ad,
            cols=cols,
            inplace_columns=tuple(inplace_columns),
            h5ad_compression=h5ad_compression,
        )
    else:
        log("[filter_status] write_inplace=false — unpurified.h5ad NOT modified")

    return sidecar_csv
