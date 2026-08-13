"""Stage: annotate raw proseg h5ad with `.obs['qc_filtered']`.

Adds a boolean per-cell column to the raw h5ad in place so downstream
consumers (`split_prep`, `rctd_prep`, and Tracy's external analysis
scripts) can filter without re-running QC. Strictly additive: no cells
are removed at this stage — only a marker column is written.

Default rule (matches the pre-refactor `preprocess` / `split_prep`
`sc.pp.filter_cells(min_counts=10)` gate): a cell is `qc_filtered=True`
when `n_counts_raw >= min_counts_cell`. `n_counts_raw` is derived
column-summing `.X` in this stage (NOT read from
`.obs['n_counts']`, which is written by scanpy's `calculate_qc_metrics`
but not guaranteed to be present on a fresh raw h5ad).

Sentinel: a sibling `.qc_filter_done.sentinel` (leading dot so it's
hidden from `ls`) next to the h5ad so a re-run of stage 1 alone
doesn't force this stage to re-run.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from xenium_preprocess._internal.compat import sentinel_exists
from xenium_preprocess._internal.layout import atomic_write_h5ad
from xenium_preprocess._internal.logging import log


def run_qc_filter(
    sample_id: str,
    raw_h5ad: Path,
    min_counts_cell: int,
    force_rerun: bool,
    qc_filtered_col: str = "qc_filtered",
) -> Path:
    """Compute `.obs[qc_filtered_col]` on the raw h5ad in place.

    Returns the h5ad path.
    """
    import anndata
    from scipy.sparse import issparse

    if not raw_h5ad.exists():
        raise SystemExit(
            f"[qc_filter] raw h5ad not found: {raw_h5ad}. "
            f"Run the proseg_to_anndata stage first."
        )

    sentinel = raw_h5ad.parent / ".qc_filter_done.sentinel"
    if sentinel_exists(sentinel, force_rerun):
        log(f"[qc_filter] sentinel exists: {sentinel} — skipping "
            f"(pass --force-rerun to re-run).")
        return raw_h5ad

    log(f"[qc_filter] reading {raw_h5ad}")
    adata = anndata.read_h5ad(raw_h5ad)
    log(f"[qc_filter] adata: {adata.n_obs} cells × {adata.n_vars} genes")

    if qc_filtered_col in adata.obs.columns and not force_rerun:
        log(f"[qc_filter] WARN: adata.obs already carries {qc_filtered_col!r}; "
            f"recomputing (sentinel was missing).")

    X = adata.X
    if issparse(X):
        # np.asarray(sparse.sum(axis=1)) returns (n, 1); squeeze to (n,).
        n_counts_raw = np.asarray(X.sum(axis=1)).squeeze(axis=1)
    else:
        n_counts_raw = np.asarray(X).sum(axis=1)
    # Cast for safety — `n_counts_raw` may be float (proseg expected-counts
    # matrix) or int (maxpost); the >= comparison works either way.
    qc_mask = n_counts_raw >= float(min_counts_cell)
    adata.obs[qc_filtered_col] = qc_mask.astype(bool)

    n_kept = int(qc_mask.sum())
    n_total = int(qc_mask.size)
    pct = 100.0 * n_kept / n_total if n_total else 0.0
    log(f"[qc_filter] min_counts_cell={min_counts_cell}: "
        f"kept={n_kept}/{n_total} ({pct:.1f}%) cells (qc_filtered=True)")

    atomic_write_h5ad(adata, raw_h5ad)
    log(f"[qc_filter] wrote {raw_h5ad}")
    sentinel.write_text("ok\n")
    log(f"[qc_filter] sentinel -> {sentinel.name}")
    return raw_h5ad
