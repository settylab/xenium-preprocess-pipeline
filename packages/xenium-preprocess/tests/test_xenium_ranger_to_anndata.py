"""`xenium_ranger_to_anndata` stage: mocked bundle → h5ad with .obsm['spatial']."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest


def _synth_bundle(tmp_path: Path, n_cells: int = 4, n_genes: int = 3) -> Path:
    """Emit a synthetic xenium-ranger bundle:
        - cell_feature_matrix.h5   (10x-format h5)
        - cells.csv.gz             (Xenium cells file)

    Skips if h5py / scipy / anndata aren't available.
    """
    try:
        import h5py  # noqa: F401
        from scipy.sparse import csr_matrix  # noqa: F401
    except ImportError as e:
        pytest.skip(f"missing dep: {e.name}")

    bundle = tmp_path / "xenium_ranger"
    bundle.mkdir()

    cell_ids = [f"cid_{i:04d}" for i in range(n_cells)]
    gene_ids = [f"GID{i}" for i in range(n_genes)]
    gene_names = [f"GENE{i}" for i in range(n_genes)]

    # Build a 10x-style h5. sc.read_10x_h5 expects a specific schema
    # (matrix/{data,indices,indptr,shape,barcodes,features/{id,name,feature_type}}).
    # Use anndata's write_h5ad + read_10x_h5 dance isn't straightforward,
    # so we build the h5 by hand.
    import h5py
    from scipy.sparse import csr_matrix

    # Genes × cells sparse int32 counts.
    counts = np.array([
        [1, 2, 3, 4],
        [0, 1, 2, 3],
        [1, 0, 5, 2],
    ], dtype=np.int32)[:n_genes, :n_cells]
    csr = csr_matrix(counts.T).T.tocsc()  # genes × cells in CSC form

    h5_path = bundle / "cell_feature_matrix.h5"
    with h5py.File(h5_path, "w") as f:
        m = f.create_group("matrix")
        m.create_dataset("barcodes", data=np.array(cell_ids, dtype="S"))
        m.create_dataset("data", data=csr.data)
        m.create_dataset("indices", data=csr.indices.astype(np.int64))
        m.create_dataset("indptr", data=csr.indptr.astype(np.int64))
        m.create_dataset("shape", data=np.array([n_genes, n_cells], dtype=np.int32))
        feat = m.create_group("features")
        feat.create_dataset("id", data=np.array(gene_ids, dtype="S"))
        feat.create_dataset("name", data=np.array(gene_names, dtype="S"))
        feat.create_dataset(
            "feature_type",
            data=np.array(["Gene Expression"] * n_genes, dtype="S"),
        )
        feat.create_dataset("genome", data=np.array(["synth"] * n_genes, dtype="S"))
        # 10x-ranger writes "target_sets" as an empty group; sc.read_10x_h5
        # doesn't require it but modern scanpy is happy with the base schema.

    # cells.csv.gz
    import pandas as pd
    cells = pd.DataFrame({
        "cell_id": cell_ids,
        "x_centroid": np.arange(n_cells, dtype=np.float64) + 100.0,
        "y_centroid": np.arange(n_cells, dtype=np.float64) + 200.0,
        "transcript_counts": np.arange(n_cells) + 10,
        "control_probe_counts": np.zeros(n_cells, dtype=np.int64),
        "total_counts": np.arange(n_cells) + 10,
        "cell_area": np.linspace(50.0, 200.0, n_cells),
        "nucleus_area": np.linspace(20.0, 80.0, n_cells),
        "nucleus_count": np.ones(n_cells, dtype=np.int64),
        "segmentation_method": ["v1"] * n_cells,
    })
    cells.to_csv(bundle / "cells.csv.gz", index=False, compression="gzip")
    return bundle


def test_xenium_ranger_to_anndata_shape_and_spatial(tmp_path: Path):
    try:
        import anndata  # noqa: F401
        import scanpy  # noqa: F401
    except ImportError as e:
        pytest.skip(f"missing dep: {e.name}")

    from xenium_preprocess.stages.xenium_ranger_to_anndata import (
        run_xenium_ranger_to_anndata,
    )
    import anndata

    bundle = _synth_bundle(tmp_path, n_cells=4, n_genes=3)
    out = tmp_path / "out" / "SYNTH_xenium_ranger.h5ad"

    run_xenium_ranger_to_anndata(
        sample_id="SYNTH",
        run_id="test_run_42",
        xenium_ranger_dir=bundle,
        out_h5ad=out,
        gex_only=True,
        force_rerun=False,
    )
    assert out.exists()
    a = anndata.read_h5ad(out)
    assert a.n_obs == 4 and a.n_vars == 3
    assert "spatial" in a.obsm
    assert a.obsm["spatial"].shape == (4, 2)
    # x_centroid values were 100, 101, 102, 103.
    assert a.obsm["spatial"][:, 0].tolist() == [100.0, 101.0, 102.0, 103.0]
    # QC cols from cells.csv landed on obs.
    for col in ("transcript_counts", "cell_area", "nucleus_area"):
        assert col in a.obs.columns


def test_xenium_ranger_to_anndata_missing_bundle_fails_loud(tmp_path: Path):
    from xenium_preprocess.stages.xenium_ranger_to_anndata import (
        run_xenium_ranger_to_anndata,
    )
    with pytest.raises(SystemExit, match="xenium_ranger_dir"):
        run_xenium_ranger_to_anndata(
            sample_id="SYNTH",
            run_id="test_run_42",
            xenium_ranger_dir=None,
            out_h5ad=tmp_path / "out.h5ad",
            gex_only=True,
            force_rerun=False,
        )


def test_xenium_ranger_to_anndata_writes_uns_identity(tmp_path: Path):
    """`.uns['sample_id']` + `.uns['run_id']` are load-bearing for H&E
    integration downstream (identity-discovery mechanism per
    reports/he-reg-xenium-integration_2026-08-11_021430_*.md).
    Skeptic-amend, TracyY123-nexus#26 comment 5251398151."""
    try:
        import anndata  # noqa: F401
        import scanpy  # noqa: F401
    except ImportError as e:
        pytest.skip(f"missing dep: {e.name}")

    from xenium_preprocess.stages.xenium_ranger_to_anndata import (
        run_xenium_ranger_to_anndata,
    )
    import anndata

    bundle = _synth_bundle(tmp_path, n_cells=3, n_genes=2)
    out = tmp_path / "out" / "MH1_xenium_ranger.h5ad"

    run_xenium_ranger_to_anndata(
        sample_id="MH1",
        run_id="my_experiment_v2",
        xenium_ranger_dir=bundle,
        out_h5ad=out,
        gex_only=True,
        force_rerun=False,
    )
    a = anndata.read_h5ad(out)
    assert a.uns.get("sample_id") == "MH1"
    assert a.uns.get("run_id") == "my_experiment_v2"
