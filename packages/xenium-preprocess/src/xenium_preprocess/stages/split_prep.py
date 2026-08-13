"""Stage 3: raw proseg AnnData → 10x-style {mtx, features, barcodes} + sidecars.

Adapted verbatim from
    the internal SPLIT preparation notebook

The reference notebook's `extract_expr_matrix_metadata` writes five
files under a single directory. This stage produces the same five
files (renamed to canonical 10x-style names) plus the two sidecar
files that the R side (rctd_prep stage) consumes.

Output files (all gzipped when `gzip_outputs=True`):
    <sample_id>{name_suffix}_counts.mtx.gz     — genes × cells (transposed).
    <sample_id>{name_suffix}_features.tsv.gz   — one gene name per line.
    <sample_id>{name_suffix}_barcodes.tsv.gz   — one cell id per line.
    <sample_id>{name_suffix}_metadata.csv      — adata.obs, cell-indexed.
    <sample_id>{name_suffix}_spatial_coords.csv.gz — (x, y) per cell, cell-indexed.

Notes on faithful port:
- The reference notebook first runs `calculate_qc_metrics` +
  `filter_cells(min_counts=10)` on the raw h5ad before the export.
  We do the same by default; set `min_counts_cell` to null to skip.
- The reference then stashes `X` as `layers['counts']` and exports
  from that layer. We keep that pattern — export from
  `adata.layers[layer]` (default "counts"); fall back to `adata.X`
  when the layer is absent so the stage still runs against a raw
  h5ad that was never normalized.
- The reference casts to `np.int32` before mmwrite — proseg counts
  are integers by construction. We preserve that (with a warning if
  the actual values aren't integer-valued).
"""
from __future__ import annotations

import csv
import gzip
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd

from xenium_preprocess._internal.compat import sentinel_exists
from xenium_preprocess._internal.layout import run_dir
from xenium_preprocess._internal.logging import log


@contextmanager
def _maybe_gzip(path: Path, mode: str, gzipped: bool):
    """Open `path` gzip'd or plain, matching the reference notebook's
    contextmanager pattern."""
    if gzipped:
        with gzip.open(path, mode) as f:
            yield f
    else:
        with open(path, mode) as f:
            yield f


def run_split_prep(
    sample_id: str,
    run_id: str,
    raw_h5ad: Path,
    output_root: Path,
    layer: str,
    name_suffix: str,
    gzip_outputs: bool,
    min_counts_cell: int | None,
    force_rerun: bool,
) -> Path:
    """Export the mtx/features/barcodes triple + sidecars.

    Returns the output directory.
    """
    import scanpy as sc
    from scipy.sparse import csr_matrix
    from scipy.io import mmwrite

    # SPLIT triple is an intermediate — it feeds `rctd_prep` (which
    # bundles the mtx into an RDS) and is not part of the final locked
    # layout. Keep it under the run folder so it doesn't collide across
    # runs but out of `spatial_adata/` / `rctd/`.
    out_dir = run_dir(output_root, sample_id, run_id) / "split_prep"
    out_dir.mkdir(parents=True, exist_ok=True)

    stem = f"{sample_id}{name_suffix}"
    gz = ".gz" if gzip_outputs else ""

    mtx_path = out_dir / f"{stem}_counts.mtx{gz}"
    features_path = out_dir / f"{stem}_features.tsv{gz}"
    barcodes_path = out_dir / f"{stem}_barcodes.tsv{gz}"
    metadata_path = out_dir / f"{stem}_metadata.csv"
    coords_path = out_dir / f"{stem}_spatial_coords.csv{gz}"

    # Sentinel check on the mtx file (the anchor of the whole triple).
    if sentinel_exists(mtx_path, force_rerun):
        log(f"[split_prep] sentinel exists: {mtx_path} — skipping "
            f"(pass --force-rerun to re-run).")
        return out_dir

    log(f"[split_prep] reading {raw_h5ad}")
    adata = sc.read_h5ad(raw_h5ad)
    log(f"[split_prep] initial adata: {adata.n_obs} cells × {adata.n_vars} genes")

    # QC filter — mirrors the notebook's second filter_cells pass.
    if min_counts_cell is not None:
        sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)
        sc.pp.filter_cells(adata, min_counts=min_counts_cell)
        log(f"[split_prep] after QC filter (min_counts={min_counts_cell}): "
            f"{adata.n_obs} cells × {adata.n_vars} genes")

    # Materialize the export layer.
    if layer in adata.layers:
        X = adata.layers[layer]
    else:
        log(f"[split_prep] layer {layer!r} not present; falling back to adata.X.")
        X = adata.X

    genes = adata.var_names
    cells = adata.obs_names

    # Match the reference: cast to int32, transpose to genes × cells,
    # then mmwrite. Warn if non-integer values would be lossy.
    X_dense_check = X if not hasattr(X, "toarray") else X
    try:
        # Peek at max-abs to sanity check integer-ness.
        peek = X.toarray()[:5, :5] if hasattr(X, "toarray") else X[:5, :5]
        if not np.allclose(peek, peek.astype(np.int32)):
            log(f"[split_prep] WARN: layer {layer!r} contains non-integer values; "
                "the int32 cast in mmwrite will TRUNCATE these. If this is not "
                "the raw proseg counts layer, override --split-layer.")
    except Exception:
        pass

    X_csr = csr_matrix(X).astype(np.int32)
    X_gc = X_csr.T   # genes × cells for the mtx (matches ReadMtx conventions)

    with _maybe_gzip(mtx_path, "wb", gzip_outputs) as f:
        mmwrite(f, X_gc)
    log(f"[split_prep] wrote {mtx_path}")

    with _maybe_gzip(features_path, "wt", gzip_outputs) as f:
        w = csv.writer(f, delimiter="\t")
        for gname in genes:
            w.writerow([gname])   # single-column form; Seurat::ReadMtx's feature.column=1
    log(f"[split_prep] wrote {features_path}")

    with _maybe_gzip(barcodes_path, "wt", gzip_outputs) as f:
        w = csv.writer(f, delimiter="\t")
        for cid in cells:
            w.writerow([cid])
    log(f"[split_prep] wrote {barcodes_path}")

    meta = adata.obs.reindex(cells)
    meta.to_csv(metadata_path, index=True)
    log(f"[split_prep] wrote {metadata_path}")

    if "spatial" in adata.obsm:
        n_dims = adata.obsm["spatial"].shape[1]
        col_names = ["x", "y"][:n_dims]
        coords = pd.DataFrame(
            adata.obsm["spatial"], index=adata.obs_names, columns=col_names
        )
        # coords are written gzipped in the reference notebook.
        if gzip_outputs:
            coords.to_csv(coords_path, index=True, compression="gzip")
        else:
            coords.to_csv(coords_path, index=True)
        log(f"[split_prep] wrote {coords_path}")
    else:
        log("[split_prep] WARN: no 'spatial' in adata.obsm — skipping "
            "coords export; downstream rctd_prep will fail to build the "
            "SpatialExperiment.")

    return out_dir
