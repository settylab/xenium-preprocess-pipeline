"""Unit tests for the postprocess stage (Stage 6).

Skipped if scanpy/anndata/scipy/umap-learn/leidenalg/igraph aren't
importable.
"""
from __future__ import annotations

from pathlib import Path

import pytest


pytest.importorskip("anndata")
pytest.importorskip("scipy")
pytest.importorskip("pandas")


RUN_ID = "42"


def _build_purified_h5ad(tmp_dir: Path, sample_id: str, seed: int = 0):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    from rctd_split._internal.layout import spatial_adata_path

    n, m = 40, 30
    rng = np.random.default_rng(seed)
    X = csr_matrix(rng.poisson(3.0, size=(n, m)).astype("float32"))
    unp = csr_matrix(rng.poisson(4.0, size=(n, m)).astype("float32"))

    cell_ids = [f"C{i:03d}" for i in range(n)]
    obs = pd.DataFrame(
        {
            "first_type": ["tumor"] * (n // 2) + ["stroma"] * (n - n // 2),
            "nCount_Proseg": np.asarray(X.sum(axis=1)).ravel().astype(int),
            "purification_status": ["purified"] * n,
            "centroid_x": rng.uniform(0, 1000, size=n),
            "centroid_y": rng.uniform(0, 1000, size=n),
        },
        index=pd.Index(cell_ids, name="cell_id"),
    )
    var = pd.DataFrame(index=[f"gene{i:02d}" for i in range(m)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["unpurified_counts"] = unp
    adata.obsm["spatial"] = obs[["centroid_x", "centroid_y"]].to_numpy(dtype=float)

    output_root = tmp_dir / "runs"
    out = spatial_adata_path(output_root, sample_id, RUN_ID, "proseg_purified")
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out, compression="gzip")
    return output_root, out


def _run_stage(output_root, sample_id, **overrides):
    pytest.importorskip("scanpy")
    pytest.importorskip("umap")
    pytest.importorskip("leidenalg")
    pytest.importorskip("igraph")
    from rctd_split.stages.postprocess import run_postprocess

    kwargs = dict(
        sample_id=sample_id,
        run_id=RUN_ID,
        output_root=output_root,
        min_counts=3,
        minprop=1.0e-3,
        positive_x=True,
        n_neighbors=5,
        random_state=0,
        leiden_resolutions=[0.3],
        process_unpurified_layer=True,
        h5ad_compression="gzip",
        force_rerun=False,
    )
    kwargs.update(overrides)
    return run_postprocess(**kwargs)


def test_postprocess_augments_purified_h5ad(tmp_path):
    import anndata as ad
    from rctd_split._internal.layout import intermediate_path

    sample_id = "MHTEST"
    output_root, purified_h5ad = _build_purified_h5ad(tmp_path, sample_id)

    sentinel = _run_stage(output_root, sample_id)
    assert sentinel == intermediate_path(
        output_root, sample_id, RUN_ID, "postprocess_sentinel",
    )
    assert sentinel.exists()

    import pandas as pd
    summary_csv = intermediate_path(
        output_root, sample_id, RUN_ID, "postprocess_summary",
    )
    assert summary_csv.exists()
    summary = pd.read_csv(summary_csv)
    assert summary.shape == (1, 10)
    assert summary.iloc[0]["sample_id"] == sample_id
    assert summary.iloc[0]["min_counts"] == 3
    assert bool(summary.iloc[0]["positive_x"]) is True

    adata = ad.read_h5ad(purified_h5ad)
    assert adata.n_obs == 40
    assert "counts" in adata.layers
    assert "X_clipped" in adata.layers
    assert "X_clipped_unpurified" in adata.layers
    assert adata.layers["X_clipped"].shape == adata.shape
    assert adata.layers["X_clipped_unpurified"].shape == adata.shape

    import numpy as np
    assert np.all(np.asarray(adata.layers["X_clipped"]) >= 0)
    assert np.all(np.asarray(adata.layers["X_clipped_unpurified"]) >= 0)

    for key in (
        "X_pca", "X_pca_unpurified",
        "X_umap_clipped_norm", "X_umap_clipped_norm_unpurified",
        "knn_indices", "knn_dists",
        "knn_indices_unpurified", "knn_dists_unpurified",
        "spatial",
    ):
        assert key in adata.obsm

    assert adata.obsm["X_umap_clipped_norm"].shape == (adata.n_obs, 2)
    assert adata.obsm["X_umap_clipped_norm_unpurified"].shape == (adata.n_obs, 2)

    assert "leiden_0.3" in adata.obs.columns
    assert str(adata.obs["leiden_0.3"].dtype) == "category"

    for col in ("first_type", "nCount_Proseg", "purification_status"):
        assert col in adata.obs.columns


def test_postprocess_sentinel_skips_second_run(tmp_path):
    sample_id = "MHTEST"
    output_root, purified_h5ad = _build_purified_h5ad(tmp_path, sample_id)

    _run_stage(output_root, sample_id)
    mtime_after_first = purified_h5ad.stat().st_mtime_ns

    _run_stage(output_root, sample_id)
    assert purified_h5ad.stat().st_mtime_ns == mtime_after_first


def test_postprocess_force_rerun_overwrites(tmp_path):
    sample_id = "MHTEST"
    output_root, purified_h5ad = _build_purified_h5ad(tmp_path, sample_id)

    _run_stage(output_root, sample_id)
    mtime_after_first = purified_h5ad.stat().st_mtime_ns

    _run_stage(output_root, sample_id, force_rerun=True)
    assert purified_h5ad.stat().st_mtime_ns != mtime_after_first


def test_postprocess_positive_x_false_allows_negatives(tmp_path):
    import anndata as ad
    import numpy as np

    sample_id = "MHTEST"
    output_root, purified_h5ad = _build_purified_h5ad(tmp_path, sample_id)

    _run_stage(output_root, sample_id, positive_x=False)
    adata = ad.read_h5ad(purified_h5ad)
    xc = np.asarray(adata.layers["X_clipped"])
    assert xc.max() <= 0.0 + 1e-6
    assert xc.min() < 0


def test_postprocess_process_unpurified_layer_false(tmp_path):
    import anndata as ad

    sample_id = "MHTEST"
    output_root, purified_h5ad = _build_purified_h5ad(tmp_path, sample_id)

    _run_stage(output_root, sample_id, process_unpurified_layer=False)
    adata = ad.read_h5ad(purified_h5ad)

    assert "X_clipped" in adata.layers
    assert "X_clipped_unpurified" not in adata.layers
    for k in (
        "X_pca_unpurified", "X_umap_clipped_norm_unpurified",
        "knn_indices_unpurified", "knn_dists_unpurified",
    ):
        assert k not in adata.obsm


def test_postprocess_missing_unpurified_layer_raises(tmp_path):
    import anndata as ad

    sample_id = "MHTEST"
    output_root, purified_h5ad = _build_purified_h5ad(tmp_path, sample_id)

    adata = ad.read_h5ad(purified_h5ad)
    del adata.layers["unpurified_counts"]
    adata.write_h5ad(purified_h5ad, compression="gzip")

    with pytest.raises(SystemExit) as exc:
        _run_stage(output_root, sample_id)
    assert "unpurified_counts" in str(exc.value)


def test_postprocess_missing_purified_h5ad_raises(tmp_path):
    from rctd_split._internal.layout import spatial_adata_dir

    sample_id = "MHTEST"
    output_root = tmp_path / "runs"
    spatial_adata_dir(output_root, sample_id, RUN_ID).mkdir(parents=True)

    with pytest.raises(SystemExit) as exc:
        _run_stage(output_root, sample_id)
    assert "mtx_to_h5ad" in str(exc.value)
