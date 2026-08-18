"""Stage: xenium-ranger bundle → `<S>_xenium_ranger.h5ad`.

Adapted from an internal recipe (reference fixture; `gex_only=True`
gives 5001 features, `False` gives 10029 including controls).

Inputs (under `<xenium_ranger_dir>`):
    - `cell_feature_matrix.h5`    — 10x-format h5 read by `sc.read_10x_h5`.
    - `cells.csv.gz`              — Xenium per-cell centroids + QC.

Outputs:
    - `<spatial_adata>/<sample_id>_xenium_ranger.h5ad`  (via atomic_write_h5ad)

Content:
    - `.X`        — raw counts from `cell_feature_matrix.h5` (gex_only=True by default).
    - `.obs_names` — Xenium cell UUIDs from the matrix.
    - `.obsm['spatial']` — (n_cells, 2) float64 x_centroid / y_centroid from `cells.csv.gz`.
    - `.obs`      — QC-relevant columns from `cells.csv.gz`
                   (`transcript_counts`, `control_probe_counts`, `total_counts`,
                    `cell_area`, `nucleus_area`, `nucleus_count`, `segmentation_method`).

The stage is minimal-at-build (no cell-type, no NN mapping). The
post-rctd-split `celltype_writeback` stage augments this file in place with
`.obs['celltype']` derived from rctd-split outputs — see
(internal issue review) §2.4. That writeback lands in a
follow-up worker.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from xenium_preprocess._internal.compat import sentinel_exists
from xenium_preprocess._internal.layout import atomic_write_h5ad
from xenium_preprocess._internal.logging import log


# Columns pulled from cells.csv.gz onto `.obs`. Missing columns are
# skipped with a WARN log rather than a hard fail — Xenium's cells.csv
# schema has drifted across ranger versions and we want the stage to
# survive a missing optional column.
_OBS_COLS = (
    "transcript_counts",
    "control_probe_counts",
    "total_counts",
    "cell_area",
    "nucleus_area",
    "nucleus_count",
    "segmentation_method",
)


def _read_cells_csv(path: Path) -> pd.DataFrame:
    """Xenium's `cells.csv.gz` — `cell_id` may be the DataFrame index
    (older exports) or a regular column (newer exports). Normalise to
    `cell_id`-indexed.
    """
    df = pd.read_csv(path)
    if "cell_id" in df.columns:
        df = df.set_index("cell_id")
    else:
        # Fall back: assume the first column is the id (older csv.gz).
        df2 = pd.read_csv(path, index_col=0)
        df2.index.name = "cell_id"
        df = df2
    return df


def run_xenium_ranger_to_anndata(
    sample_id: str,
    run_id: str,
    xenium_ranger_dir: Path,
    out_h5ad: Path,
    gex_only: bool,
    force_rerun: bool,
    cells_csv_name: str = "cells.csv.gz",
    cell_feature_matrix_name: str = "cell_feature_matrix.h5",
) -> Path:
    """Read the xenium-ranger bundle, write `<sample_id>_xenium_ranger.h5ad`.

    Returns the output h5ad path.
    """
    import anndata
    import scanpy as sc

    if xenium_ranger_dir is None:
        raise SystemExit(
            "[xenium_ranger_to_anndata] xenium_ranger_dir is null — pass "
            "--xenium-ranger-dir <path/to/xenium_ranger/<sample>/> or set "
            "config.xenium_ranger_to_anndata.xenium_ranger_dir. To skip this "
            "stage entirely, drop it from --stages."
        )
    xenium_ranger_dir = Path(xenium_ranger_dir)
    if not xenium_ranger_dir.exists():
        raise SystemExit(
            f"[xenium_ranger_to_anndata] xenium_ranger_dir does not exist: "
            f"{xenium_ranger_dir}"
        )

    h5_path = xenium_ranger_dir / cell_feature_matrix_name
    cells_csv = xenium_ranger_dir / cells_csv_name
    if not h5_path.exists():
        raise SystemExit(
            f"[xenium_ranger_to_anndata] cell_feature_matrix.h5 not found: "
            f"{h5_path}"
        )
    if not cells_csv.exists():
        raise SystemExit(
            f"[xenium_ranger_to_anndata] cells.csv.gz not found: {cells_csv}"
        )

    if sentinel_exists(out_h5ad, force_rerun):
        log(f"[xenium_ranger_to_anndata] sentinel exists: {out_h5ad} — "
            f"skipping (pass --force-rerun to re-run).")
        return out_h5ad

    log(f"[xenium_ranger_to_anndata] cell_feature_matrix: {h5_path} "
        f"(gex_only={gex_only})")
    adata = sc.read_10x_h5(h5_path, gex_only=gex_only)
    # anndata unique-var warning on non-unique gene symbols is loud but
    # non-fatal; scanpy already logs it. We accept whatever names come
    # out of read_10x_h5 as-is.
    adata.var_names_make_unique()
    log(f"[xenium_ranger_to_anndata] adata: {adata.n_obs} cells × "
        f"{adata.n_vars} genes")

    log(f"[xenium_ranger_to_anndata] cells csv: {cells_csv}")
    cells = _read_cells_csv(cells_csv)
    log(f"[xenium_ranger_to_anndata] cells.csv: {len(cells)} rows, "
        f"columns={list(cells.columns)[:12]}"
        f"{'...' if len(cells.columns) > 12 else ''}")

    # Sanity: every adata.obs_name should exist in cells.csv.
    missing = adata.obs_names.difference(cells.index)
    if len(missing) > 0:
        raise SystemExit(
            f"[xenium_ranger_to_anndata] {len(missing)} adata cell_ids "
            f"missing from cells.csv (first 5: {list(missing[:5])}). "
            f"Refusing to build the h5ad — the two bundle files are "
            f"supposed to describe the same cells."
        )

    # Align cells.csv rows to adata order.
    cells_aligned = cells.reindex(adata.obs_names)

    # Spatial obsm: (n_cells, 2) float64 x_centroid / y_centroid.
    for col in ("x_centroid", "y_centroid"):
        if col not in cells_aligned.columns:
            raise SystemExit(
                f"[xenium_ranger_to_anndata] cells.csv missing required "
                f"column {col!r}. Present: {list(cells_aligned.columns)}"
            )
    xy = np.stack(
        [
            cells_aligned["x_centroid"].astype("float64").to_numpy(),
            cells_aligned["y_centroid"].astype("float64").to_numpy(),
        ],
        axis=1,
    )
    n_nan = int(np.isnan(xy).any(axis=1).sum())
    if n_nan > 0:
        raise SystemExit(
            f"[xenium_ranger_to_anndata] .obsm['spatial'] has {n_nan} rows "
            f"with NaN centroid — refusing to write. Inspect cells.csv."
        )
    adata.obsm["spatial"] = xy

    # QC obs columns.
    for col in _OBS_COLS:
        if col in cells_aligned.columns:
            adata.obs[col] = cells_aligned[col].to_numpy()
        else:
            log(f"[xenium_ranger_to_anndata] WARN: cells.csv missing "
                f"optional obs col {col!r}; skipping.")

    log(f"[xenium_ranger_to_anndata] .obsm['spatial']: shape={adata.obsm['spatial'].shape} "
        f"x=[{xy[:,0].min():.1f},{xy[:,0].max():.1f}] "
        f"y=[{xy[:,1].min():.1f},{xy[:,1].max():.1f}]")

    # Identity in .uns — load-bearing for H&E-integration downstream:
    # reports/he-reg-xenium-integration_2026-08-11_021430_*.md locks
    # `.uns['sample_id']` / `.uns['run_id']` as the discovery mechanism
    # so registration doesn't have to parse the enclosing folder name.
    adata.uns["sample_id"] = sample_id
    adata.uns["run_id"] = run_id

    atomic_write_h5ad(adata, out_h5ad)
    log(f"[xenium_ranger_to_anndata] wrote {out_h5ad}")
    return out_h5ad
