"""`qc_filter` stage: `.obs['qc_filtered']` matches
`sc.pp.filter_cells(min_counts=…)`'s barcode set."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest


def _write_synth_h5ad(tmp_path: Path) -> Path:
    """Build a small AnnData with known row-sums, write to h5ad."""
    try:
        import anndata
        from scipy.sparse import csr_matrix
    except ImportError as e:
        pytest.skip(f"missing dep: {e.name}")

    # 5 cells × 4 genes. Row sums: 5, 12, 3, 20, 10.
    X = np.array([
        [1, 2, 1, 1],
        [3, 3, 3, 3],
        [1, 1, 1, 0],
        [5, 5, 5, 5],
        [2, 3, 2, 3],
    ], dtype=np.int32)
    adata = anndata.AnnData(X=csr_matrix(X))
    adata.obs_names = [f"cell_{i}" for i in range(5)]
    adata.var_names = [f"g{i}" for i in range(4)]
    path = tmp_path / "raw.h5ad"
    adata.write_h5ad(path)
    return path


def test_qc_filter_writes_boolean_column_matching_min_counts_gate(tmp_path: Path):
    from xenium_preprocess.stages.qc_filter import run_qc_filter
    import anndata

    raw = _write_synth_h5ad(tmp_path)
    run_qc_filter(sample_id="SYNTH", raw_h5ad=raw, min_counts_cell=10,
                  force_rerun=False)
    a = anndata.read_h5ad(raw)
    assert "qc_filtered" in a.obs.columns
    # Row sums were [5, 12, 3, 20, 10]; gate is >=10 → [F, T, F, T, T].
    assert a.obs["qc_filtered"].tolist() == [False, True, False, True, True]


def test_qc_filter_matches_scanpy_filter_cells(tmp_path: Path):
    """The barcodes that pass the sc.pp.filter_cells gate must be
    exactly the ones marked qc_filtered=True."""
    try:
        import scanpy as sc
    except ImportError as e:
        pytest.skip(f"missing dep: {e.name}")
    import anndata
    from xenium_preprocess.stages.qc_filter import run_qc_filter

    raw = _write_synth_h5ad(tmp_path)

    # Reference: sc.pp.filter_cells against a fresh read.
    a_ref = anndata.read_h5ad(raw)
    sc.pp.filter_cells(a_ref, min_counts=10)
    ref_barcodes = set(a_ref.obs_names.tolist())

    # Our qc_filter run.
    run_qc_filter(sample_id="SYNTH", raw_h5ad=raw, min_counts_cell=10,
                  force_rerun=False)
    a = anndata.read_h5ad(raw)
    ours_barcodes = set(a.obs_names[a.obs["qc_filtered"].to_numpy()])
    assert ours_barcodes == ref_barcodes


def test_qc_filter_sentinel_skips_second_run(tmp_path: Path):
    from xenium_preprocess.stages.qc_filter import run_qc_filter
    import anndata

    raw = _write_synth_h5ad(tmp_path)
    run_qc_filter(sample_id="SYNTH", raw_h5ad=raw, min_counts_cell=10,
                  force_rerun=False)
    assert (raw.parent / ".qc_filter_done.sentinel").exists()

    # Overwrite qc_filtered externally; re-run without force should skip.
    a = anndata.read_h5ad(raw)
    a.obs["qc_filtered"] = False  # sabotage
    a.write_h5ad(raw)

    run_qc_filter(sample_id="SYNTH", raw_h5ad=raw, min_counts_cell=10,
                  force_rerun=False)
    a2 = anndata.read_h5ad(raw)
    # Because sentinel exists → stage skipped → the sabotaged False stays.
    assert (a2.obs["qc_filtered"] == False).all()  # noqa: E712
