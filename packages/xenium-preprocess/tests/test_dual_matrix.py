"""Dual-matrix support tests (Tracy request 2026-07-10, TracyY123-nexus#14).

Two tests:

1. `proseg_to_anndata` loads BOTH matrices into layers.
2. `run_preprocess(dual_matrix_mode=True)` populates both suffixed
   PCA/UMAP/Leiden outputs.

Both tests skip cleanly on any missing science-stack dependency
(`anndata`, `scipy`, `scanpy`, `umap-learn`, `leidenalg`, `igraph`) so
the smoke suite still runs on a stock Python-plus-pytest env.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest


# ---------------------------------------------------------------------
# proseg_to_anndata dual-matrix round-trip
# ---------------------------------------------------------------------

def _synth_proseg_dir(tmp_path: Path, n_cells: int = 12, n_genes: int = 8) -> Path:
    """Write a fake proseg output tree.

    - `expected-counts.csv` — cells × genes, float (uniform random 0..3).
    - `maxpost_counts.csv` — cells × genes, int (Poisson-style small ints),
      SAME gene column names + order as expected-counts.
    - `cell-metadata.csv` — carries `cell`, `centroid_x`, `centroid_y`.

    Matches the shape proseg emits (rows = cells; no explicit cell-id
    column in the count matrix; cell ordering matches cell-metadata).
    """
    rng = np.random.default_rng(0)
    proseg_dir = tmp_path / "proseg"
    proseg_dir.mkdir()

    gene_names = [f"GENE{i:03d}" for i in range(n_genes)]

    expected = rng.uniform(0.0, 3.0, size=(n_cells, n_genes)).astype(np.float32)
    maxpost = rng.integers(0, 5, size=(n_cells, n_genes)).astype(np.int32)

    pd.DataFrame(expected, columns=gene_names).to_csv(
        proseg_dir / "expected-counts.csv", index=False
    )
    pd.DataFrame(maxpost, columns=gene_names).to_csv(
        proseg_dir / "maxpost_counts.csv", index=False
    )
    pd.DataFrame({
        "cell": [f"cell_{i:04d}" for i in range(n_cells)],
        "centroid_x": rng.uniform(0.0, 1000.0, size=n_cells),
        "centroid_y": rng.uniform(0.0, 1000.0, size=n_cells),
    }).to_csv(proseg_dir / "cell-metadata.csv", index=False)

    return proseg_dir


def test_proseg_to_anndata_loads_both_matrices(tmp_path: Path):
    """Both `expected_counts` and `maxpost_counts` land as layers with
    matching shape, and `.X` mirrors the DEFAULT source (maxpost_counts,
    per the 2026-07-23 default flip — Tracy request TracyY123-nexus#14
    comment 5064106232)."""
    try:
        import anndata
        import scipy  # noqa: F401
    except ImportError as e:
        pytest.skip(f"missing dep: {e.name}")

    from xenium_preprocess.stages.proseg_to_anndata import run_proseg_to_anndata

    n_cells, n_genes = 12, 8
    proseg_dir = _synth_proseg_dir(tmp_path, n_cells=n_cells, n_genes=n_genes)
    out_root = tmp_path / "out"

    out_h5ad = run_proseg_to_anndata(
        sample_id="SYNTH",
        run_id="test_run",
        proseg_dir=proseg_dir,
        output_root=out_root,
        count_matrix_path=None,
        cell_metadata_path=None,
        maxpost_matrix_path=None,
        count_matrix_glob="expected-counts*",
        cell_metadata_glob="cell-metadata*",
        maxpost_matrix_glob="maxpost_counts*",
        centroid_x_col="centroid_x",
        centroid_y_col="centroid_y",
        cell_id_col="cell",
        force_rerun=False,
    )
    assert out_h5ad.exists()

    adata = anndata.read_h5ad(out_h5ad)
    assert adata.n_obs == n_cells and adata.n_vars == n_genes
    # Both layers must be present.
    assert "expected_counts" in adata.layers
    assert "maxpost_counts" in adata.layers
    # Shape parity across `.X` + both layers.
    assert adata.X.shape == (n_cells, n_genes)
    assert adata.layers["expected_counts"].shape == (n_cells, n_genes)
    assert adata.layers["maxpost_counts"].shape == (n_cells, n_genes)
    # `.X` mirrors maxpost_counts by default (2026-07-23 flip).
    Xe = adata.layers["expected_counts"]
    Xm = adata.layers["maxpost_counts"]
    Xe_arr = Xe.toarray() if hasattr(Xe, "toarray") else np.asarray(Xe)
    Xm_arr = Xm.toarray() if hasattr(Xm, "toarray") else np.asarray(Xm)
    X_arr = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    assert np.array_equal(Xm_arr, X_arr), (
        "`.X` should mirror layers['maxpost_counts'] under the 2026-07-23 default"
    )
    # The two layers describe DIFFERENT numbers (not a copy-paste bug).
    assert not np.array_equal(Xe_arr, Xm_arr)


def test_proseg_to_anndata_x_source_expected_flips_back(tmp_path: Path):
    """`x_source='expected_counts'` restores the pre-2026-07-23 layout
    (`.X` = expected). Reverse-knob back-compat lever."""
    try:
        import anndata
        import scipy  # noqa: F401
    except ImportError as e:
        pytest.skip(f"missing dep: {e.name}")

    from xenium_preprocess.stages.proseg_to_anndata import run_proseg_to_anndata

    n_cells, n_genes = 12, 8
    proseg_dir = _synth_proseg_dir(tmp_path, n_cells=n_cells, n_genes=n_genes)
    out_root = tmp_path / "out"

    out_h5ad = run_proseg_to_anndata(
        sample_id="SYNTH",
        run_id="test_run",
        proseg_dir=proseg_dir,
        output_root=out_root,
        count_matrix_path=None,
        cell_metadata_path=None,
        maxpost_matrix_path=None,
        count_matrix_glob="expected-counts*",
        cell_metadata_glob="cell-metadata*",
        maxpost_matrix_glob="maxpost_counts*",
        centroid_x_col="centroid_x",
        centroid_y_col="centroid_y",
        cell_id_col="cell",
        force_rerun=False,
        x_source="expected_counts",
    )
    adata = anndata.read_h5ad(out_h5ad)
    # Both layers still present — only `.X` mirror changes.
    assert "expected_counts" in adata.layers
    assert "maxpost_counts" in adata.layers
    Xe = adata.layers["expected_counts"]
    Xm = adata.layers["maxpost_counts"]
    Xe_arr = Xe.toarray() if hasattr(Xe, "toarray") else np.asarray(Xe)
    Xm_arr = Xm.toarray() if hasattr(Xm, "toarray") else np.asarray(Xm)
    X_arr = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    assert np.array_equal(Xe_arr, X_arr), (
        "with x_source='expected_counts', `.X` should mirror expected_counts"
    )
    # Sanity: `.X` should NOT be maxpost under the reverse knob.
    assert not np.array_equal(Xm_arr, X_arr)


def test_proseg_to_anndata_x_source_invalid_fails_loud(tmp_path: Path):
    """`x_source='bogus'` → SystemExit — no silent surprise layout."""
    try:
        import anndata  # noqa: F401
        import scipy  # noqa: F401
    except ImportError as e:
        pytest.skip(f"missing dep: {e.name}")

    from xenium_preprocess.stages.proseg_to_anndata import run_proseg_to_anndata

    proseg_dir = _synth_proseg_dir(tmp_path, n_cells=5, n_genes=4)
    out_root = tmp_path / "out"

    with pytest.raises(SystemExit):
        run_proseg_to_anndata(
            sample_id="SYNTH",
            run_id="test_run",
            proseg_dir=proseg_dir,
            output_root=out_root,
            count_matrix_path=None,
            cell_metadata_path=None,
            maxpost_matrix_path=None,
            count_matrix_glob="expected-counts*",
            cell_metadata_glob="cell-metadata*",
            maxpost_matrix_glob="maxpost_counts*",
            centroid_x_col="centroid_x",
            centroid_y_col="centroid_y",
            cell_id_col="cell",
            force_rerun=False,
            x_source="bogus",
        )


def test_proseg_to_anndata_maxpost_default_null_glob_falls_back_to_expected(
    tmp_path: Path,
):
    """Default `x_source='maxpost_counts'` + `maxpost_matrix_glob=None` →
    `.X` silently falls back to expected_counts with a WARN log.
    Preserves back-compat for pre-dual-mode user YAMLs."""
    try:
        import anndata
        import scipy  # noqa: F401
    except ImportError as e:
        pytest.skip(f"missing dep: {e.name}")

    from xenium_preprocess.stages.proseg_to_anndata import run_proseg_to_anndata

    proseg_dir = _synth_proseg_dir(tmp_path, n_cells=5, n_genes=4)
    out_root = tmp_path / "out"

    out_h5ad = run_proseg_to_anndata(
        sample_id="SYNTH",
        run_id="test_run",
        proseg_dir=proseg_dir,
        output_root=out_root,
        count_matrix_path=None,
        cell_metadata_path=None,
        maxpost_matrix_path=None,
        count_matrix_glob="expected-counts*",
        cell_metadata_glob="cell-metadata*",
        maxpost_matrix_glob=None,       # <-- disable maxpost loading
        centroid_x_col="centroid_x",
        centroid_y_col="centroid_y",
        cell_id_col="cell",
        force_rerun=False,
        # x_source omitted → defaults to "maxpost_counts"
    )
    adata = anndata.read_h5ad(out_h5ad)
    assert "expected_counts" in adata.layers
    assert "maxpost_counts" not in adata.layers
    # `.X` falls back to expected.
    Xe = adata.layers["expected_counts"]
    Xe_arr = Xe.toarray() if hasattr(Xe, "toarray") else np.asarray(Xe)
    X_arr = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    assert np.array_equal(Xe_arr, X_arr)


def test_proseg_to_anndata_null_maxpost_reverts_to_single(tmp_path: Path):
    """`maxpost_matrix_glob=None` skips the maxpost load — back-compat."""
    try:
        import anndata
        import scipy  # noqa: F401
    except ImportError as e:
        pytest.skip(f"missing dep: {e.name}")

    from xenium_preprocess.stages.proseg_to_anndata import run_proseg_to_anndata

    proseg_dir = _synth_proseg_dir(tmp_path, n_cells=5, n_genes=4)
    out_root = tmp_path / "out"

    out_h5ad = run_proseg_to_anndata(
        sample_id="SYNTH",
        run_id="test_run",
        proseg_dir=proseg_dir,
        output_root=out_root,
        count_matrix_path=None,
        cell_metadata_path=None,
        maxpost_matrix_path=None,
        count_matrix_glob="expected-counts*",
        cell_metadata_glob="cell-metadata*",
        maxpost_matrix_glob=None,   # <-- disable
        centroid_x_col="centroid_x",
        centroid_y_col="centroid_y",
        cell_id_col="cell",
        force_rerun=False,
    )
    adata = anndata.read_h5ad(out_h5ad)
    assert "maxpost_counts" not in adata.layers
    assert "expected_counts" in adata.layers


def test_proseg_to_anndata_shape_mismatch_fails_loud(tmp_path: Path):
    """Bogus maxpost matrix with wrong shape → SystemExit, not a
    silent, mis-aligned h5ad."""
    try:
        import anndata  # noqa: F401
        import scipy  # noqa: F401
    except ImportError as e:
        pytest.skip(f"missing dep: {e.name}")

    from xenium_preprocess.stages.proseg_to_anndata import run_proseg_to_anndata

    proseg_dir = _synth_proseg_dir(tmp_path, n_cells=10, n_genes=6)
    # Overwrite maxpost with a shape-mismatched CSV.
    pd.DataFrame(
        np.zeros((10, 5), dtype=int),
        columns=[f"GENE{i:03d}" for i in range(5)],
    ).to_csv(proseg_dir / "maxpost_counts.csv", index=False)

    out_root = tmp_path / "out"

    with pytest.raises(SystemExit):
        run_proseg_to_anndata(
            sample_id="SYNTH",
            run_id="test_run",
            proseg_dir=proseg_dir,
            output_root=out_root,
            count_matrix_path=None,
            cell_metadata_path=None,
            maxpost_matrix_path=None,
            count_matrix_glob="expected-counts*",
            cell_metadata_glob="cell-metadata*",
            maxpost_matrix_glob="maxpost_counts*",
            centroid_x_col="centroid_x",
            centroid_y_col="centroid_y",
            cell_id_col="cell",
            force_rerun=False,
        )


# ---------------------------------------------------------------------
# preprocess dual-mode: both suffixed output keys land.
# ---------------------------------------------------------------------

def _dual_layer_adata(n_cells: int = 30, n_genes: int = 20):
    """Build a synthetic AnnData with two count layers — the shape
    `run_preprocess` expects after stage 1 with dual-matrix support.
    """
    import anndata
    from scipy.sparse import csr_matrix

    rng = np.random.default_rng(1)
    # Draw two DIFFERENT count matrices so leiden partitions plausibly
    # differ. Keep the totals well above the min_counts_cell=10 filter.
    # Use non-Neg / non-Unassigned gene names so remove_control_probes
    # doesn't drop them.
    genes = [f"GN{i:03d}" for i in range(n_genes)]
    expected = rng.poisson(2.0, size=(n_cells, n_genes)).astype(np.float32) + 0.1
    maxpost = rng.poisson(2.0, size=(n_cells, n_genes)).astype(np.float32)

    var = pd.DataFrame(index=genes)
    obs = pd.DataFrame({"dummy": range(n_cells)},
                       index=[f"c_{i:04d}" for i in range(n_cells)])

    adata = anndata.AnnData(
        X=csr_matrix(expected), var=var, obs=obs,
        layers={
            "expected_counts": csr_matrix(expected),
            "maxpost_counts": csr_matrix(maxpost),
        },
        obsm={"spatial": rng.uniform(0, 100, size=(n_cells, 2))},
    )
    return adata


def test_run_preprocess_dual_mode_populates_suffixed_keys(tmp_path: Path):
    """dual_matrix_mode=True → both `_expected` and `_maxpost` suffixed
    obsm/obs keys land."""
    try:
        import anndata  # noqa: F401
        import scanpy  # noqa: F401
        import umap  # noqa: F401
        import leidenalg  # noqa: F401
        import igraph  # noqa: F401
    except ImportError as e:
        pytest.skip(f"missing dep: {e.name}")

    from xenium_preprocess.stages.preprocess import run_preprocess

    adata = _dual_layer_adata(n_cells=30, n_genes=20)
    # Mirror the pipeline's on-disk shape: raw h5ad inside a `spatial_adata/`
    # subdir so `raw_h5ad.parent.parent` — where the deprecated preprocess
    # stage puts `legacy_preprocess/` — resolves to the run-dir level.
    run_dir = tmp_path / "SYNTH" / "SYNTH_test_run"
    (run_dir / "spatial_adata").mkdir(parents=True)
    raw_h5ad = run_dir / "spatial_adata" / "raw.h5ad"
    adata.write_h5ad(raw_h5ad)

    out_root = tmp_path / "out"
    resolution = 0.4

    out_h5ad = run_preprocess(
        sample_id="SYNTH",
        run_id="test_run",
        raw_h5ad=raw_h5ad,
        output_root=out_root,
        qc_percent_top=[10, 20],
        min_counts_cell=1,   # keep every synthetic cell
        min_prop=1e-3,
        positive_x=True,
        n_neighbors=5,
        random_state=0,
        leiden_resolution=resolution,
        use_leiden_weights=False,
        global_non_tumor_json=None,
        global_tumor_json=None,
        tumor_type=None,
        threshold_global=-0.005,
        low_count_label="low count",
        global_level1_celltype_col="global_level1_celltype",
        umap_dpi=100,
        umap_figsize=(4, 3),
        force_rerun=False,
        dual_matrix_mode=True,
    )
    assert out_h5ad.exists()

    import anndata
    ad = anndata.read_h5ad(out_h5ad)

    # PCA + UMAP + KNN keys for both passes.
    for suffix in ("_expected", "_maxpost"):
        assert f"X_pca{suffix}" in ad.obsm, (
            f"missing obsm key X_pca{suffix} — dual pass didn't complete")
        assert f"X_umap{suffix}" in ad.obsm
        assert f"knn_indices{suffix}" in ad.obsm
        assert f"X_clipped{suffix}" in ad.layers
        leiden_key = f"leiden_{resolution}{suffix}"
        assert leiden_key in ad.obs, (
            f"missing obs key {leiden_key} — Leiden didn't run for this pass")

    # UMAP plots for both passes should exist. Preprocess is DEPRECATED
    # and now writes into a `<run_dir>/legacy_preprocess/plots/` folder —
    # under the run folder but out of the locked `spatial_adata/` shape.
    plots_dir = run_dir / "legacy_preprocess" / "plots"
    assert (plots_dir / "umap_leiden_expected.png").exists()
    assert (plots_dir / "umap_leiden_maxpost.png").exists()


def test_run_preprocess_single_mode_falls_back(tmp_path: Path):
    """dual_matrix_mode=False → unsuffixed keys populate; no `_expected`
    / `_maxpost` collision. Back-compat with pre-2026-07-10."""
    try:
        import anndata  # noqa: F401
        import scanpy  # noqa: F401
        import umap  # noqa: F401
        import leidenalg  # noqa: F401
        import igraph  # noqa: F401
    except ImportError as e:
        pytest.skip(f"missing dep: {e.name}")

    from xenium_preprocess.stages.preprocess import run_preprocess

    adata = _dual_layer_adata(n_cells=25, n_genes=15)
    run_dir = tmp_path / "SYNTH" / "SYNTH_test_run"
    (run_dir / "spatial_adata").mkdir(parents=True)
    raw_h5ad = run_dir / "spatial_adata" / "raw.h5ad"
    adata.write_h5ad(raw_h5ad)

    out_root = tmp_path / "out"

    out_h5ad = run_preprocess(
        sample_id="SYNTH",
        run_id="test_run",
        raw_h5ad=raw_h5ad,
        output_root=out_root,
        qc_percent_top=[5, 10],
        min_counts_cell=1,
        min_prop=1e-3,
        positive_x=True,
        n_neighbors=5,
        random_state=0,
        leiden_resolution=0.4,
        use_leiden_weights=False,
        global_non_tumor_json=None,
        global_tumor_json=None,
        tumor_type=None,
        threshold_global=-0.005,
        low_count_label="low count",
        global_level1_celltype_col="global_level1_celltype",
        umap_dpi=100,
        umap_figsize=(4, 3),
        force_rerun=False,
        dual_matrix_mode=False,
    )

    import anndata
    ad = anndata.read_h5ad(out_h5ad)

    # Unsuffixed keys populated.
    assert "X_pca" in ad.obsm
    assert "X_umap" in ad.obsm
    assert "leiden_0.4" in ad.obs
    # No dual-mode keys.
    assert "X_pca_expected" not in ad.obsm
    assert "X_pca_maxpost" not in ad.obsm
