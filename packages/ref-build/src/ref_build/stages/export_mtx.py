"""Stage 4: assembled reference AnnData → 10X-style bundle.

Writes the four files that both the Stage C R script and any downstream
RCTD workflow can read:

    <sample>_counts.mtx.gz     — MatrixMarket, gzipped, int32,
                                 transposed to genes × cells.
    <sample>_features.tsv.gz   — one gene name per line
                                 (single-column form; matches
                                  Seurat::ReadMtx(feature.column=1)).
    <sample>_barcodes.tsv.gz   — one cell id per line.
    <sample>_metadata.csv      — adata.obs, cell-indexed, plain CSV.

File naming convention mirrors xenium-preprocess's `split_prep` output shape so
downstream R code sees a consistent layout across pipelines.
"""
from __future__ import annotations

import csv
import gzip
from contextlib import contextmanager
from pathlib import Path

import numpy as np

from ref_build._internal.compat import sentinel_exists
from ref_build._internal.logging import log


@contextmanager
def _maybe_gzip(path: Path, mode: str, gzipped: bool):
    """Open `path` gzip'd or plain."""
    if gzipped:
        with gzip.open(path, mode) as f:
            yield f
    else:
        with open(path, mode) as f:
            yield f


def _pick_export_matrix(adata, layer_hint: str):
    """Return (matrix, used_layer_name) for the export.

    Tries `layer_hint` (default `counts`), then `raw_count`, then `.X`.
    Same shape as `assemble._pick_counts`.
    """
    for name in (layer_hint, "raw_count"):
        if name and name in adata.layers:
            return adata.layers[name], name
    return adata.X, "X"


def run_export_mtx(
    sample_id: str,
    reference_h5ad: Path,
    output_root: Path,
    layer: str,
    gzip_outputs: bool,
    force_rerun: bool,
    run_id: str,
) -> Path:
    """Write the 10X-style bundle for `sample_id` at
    `<run_dir>/mtx_bundle/` (intermediate — dropped after
    `rctd_reference_build` succeeds).

    Returns the output directory.
    """
    import anndata as ad
    from scipy.sparse import csr_matrix
    from scipy.io import mmwrite

    from ref_build._internal.layout import run_dir as _run_dir

    out_dir = _run_dir(output_root, sample_id, run_id) / "mtx_bundle"
    out_dir.mkdir(parents=True, exist_ok=True)

    gz = ".gz" if gzip_outputs else ""
    mtx_path = out_dir / f"{sample_id}_counts.mtx{gz}"
    features_path = out_dir / f"{sample_id}_features.tsv{gz}"
    barcodes_path = out_dir / f"{sample_id}_barcodes.tsv{gz}"
    metadata_path = out_dir / f"{sample_id}_metadata.csv"

    # Sentinel check on the mtx (anchor of the bundle).
    if sentinel_exists(mtx_path, force_rerun):
        log(f"[export_mtx] sentinel exists: {mtx_path} — skipping "
            f"(pass --force-rerun to re-run).")
        return out_dir

    log(f"[export_mtx] reading {reference_h5ad}")
    adata = ad.read_h5ad(reference_h5ad)
    log(f"[export_mtx]   reference shape: {adata.n_obs} cells × {adata.n_vars} genes")

    X, used_layer = _pick_export_matrix(adata, layer)
    log(f"[export_mtx]   counts source layer: {used_layer!r}")

    genes = adata.var_names
    cells = adata.obs_names

    # Cast to int32; transpose to genes × cells (matches ReadMtx conventions).
    X_csr = csr_matrix(X).astype(np.int32)
    X_gc = X_csr.T

    with _maybe_gzip(mtx_path, "wb", gzip_outputs) as f:
        mmwrite(f, X_gc)
    log(f"[export_mtx] wrote {mtx_path}")

    # Single-column features.tsv — matches Seurat::ReadMtx(feature.column=1).
    with _maybe_gzip(features_path, "wt", gzip_outputs) as f:
        w = csv.writer(f, delimiter="\t")
        for gname in genes:
            w.writerow([gname])
    log(f"[export_mtx] wrote {features_path}")

    with _maybe_gzip(barcodes_path, "wt", gzip_outputs) as f:
        w = csv.writer(f, delimiter="\t")
        for cid in cells:
            w.writerow([cid])
    log(f"[export_mtx] wrote {barcodes_path}")

    meta = adata.obs.reindex(cells)
    meta.to_csv(metadata_path, index=True)
    log(f"[export_mtx] wrote {metadata_path}")

    return out_dir
