"""Round-trip test for the mtx_to_h5ad stage.

Synthetic mtx bundle (100 cells × 50 genes, integer counts, spatial
coords + metadata) → mtx_to_h5ad → AnnData with the expected n_obs /
n_vars / .obs columns / .obsm['spatial'] shape.

No R involved. Skipped if anndata + scanpy + scipy aren't importable.
"""
from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pytest


pytest.importorskip("anndata")
pytest.importorskip("scanpy")
pytest.importorskip("scipy")


RUN_ID = "42"


def _write_bundle(tmp_dir: Path, sample_id: str, variant: str,
                  n_cells: int, n_genes: int, gz: bool = True):
    """Create a synthetic 10X-style mtx bundle under `tmp_dir`."""
    import numpy as np
    from scipy.sparse import csr_matrix
    from scipy.io import mmwrite
    import pandas as pd

    rng = np.random.default_rng(42)
    counts = rng.poisson(lam=2.0, size=(n_genes, n_cells)).astype("int32")
    counts_csr = csr_matrix(counts)

    stem = tmp_dir / f"{sample_id}_{variant}"

    if gz:
        with gzip.open(f"{stem}_counts.mtx.gz", "wb") as f:
            mmwrite(f, counts_csr)
    else:
        mmwrite(str(f"{stem}_counts.mtx"), counts_csr)

    gene_names = [f"gene{i}" for i in range(n_genes)]
    if gz:
        with gzip.open(f"{stem}_features.tsv.gz", "wt") as f:
            w = csv.writer(f, delimiter="\t")
            for g in gene_names:
                w.writerow([g])
    else:
        with open(f"{stem}_features.tsv", "w") as f:
            w = csv.writer(f, delimiter="\t")
            for g in gene_names:
                w.writerow([g])

    barcodes = [f"cell{i}" for i in range(n_cells)]
    if gz:
        with gzip.open(f"{stem}_barcodes.tsv.gz", "wt") as f:
            w = csv.writer(f, delimiter="\t")
            for b in barcodes:
                w.writerow([b])
    else:
        with open(f"{stem}_barcodes.tsv", "w") as f:
            w = csv.writer(f, delimiter="\t")
            for b in barcodes:
                w.writerow([b])

    meta = pd.DataFrame(
        {
            "first_type": rng.choice(["tumor", "Fibroblast", "endothelial"], size=n_cells),
            "purification_status": rng.choice(
                ["singlet", "doublet_certain", "purified"], size=n_cells
            ),
            "w1_larger_w2": rng.choice([True, False], size=n_cells),
            "same_class": rng.choice([True, False], size=n_cells),
            "nCount_Proseg": counts.sum(axis=0),
            "x": rng.uniform(0, 1000, size=n_cells),
            "y": rng.uniform(0, 1000, size=n_cells),
            "cell_area": rng.uniform(10, 500, size=n_cells),
            "original_cell_id": [f"orig_{i}" for i in range(n_cells)],
        },
        index=barcodes,
    )
    meta.to_csv(f"{stem}_metadata.csv", index=True)

    coords = pd.DataFrame(
        {"x": meta["x"].to_numpy(), "y": meta["y"].to_numpy()},
        index=barcodes,
    )
    if gz:
        coords.to_csv(f"{stem}_spatial_coords.csv.gz", index=True, compression="gzip")
    else:
        coords.to_csv(f"{stem}_spatial_coords.csv", index=True)

    return tmp_dir


def _prepare_bundle_dirs(output_root: Path, sample_id: str):
    """Create the intermediate mtx bundle dirs the pipeline expects."""
    from rctd_split._internal.layout import mtx_bundle_dir

    unp_dir = mtx_bundle_dir(output_root, sample_id, RUN_ID, "unpurified")
    pur_dir = mtx_bundle_dir(output_root, sample_id, RUN_ID, "purified")
    unp_dir.mkdir(parents=True)
    pur_dir.mkdir(parents=True)
    return unp_dir, pur_dir


def test_mtx_to_h5ad_roundtrip(tmp_path):
    from rctd_split.stages.mtx_to_h5ad import run_mtx_to_h5ad

    sample_id = "MHTEST"
    n_cells, n_genes = 100, 50
    output_root = tmp_path / "runs"
    unp_dir, pur_dir = _prepare_bundle_dirs(output_root, sample_id)

    _write_bundle(unp_dir, sample_id, "unpurified", n_cells, n_genes, gz=True)
    _write_bundle(pur_dir, sample_id, "purified", n_cells, n_genes, gz=True)

    unp_h5ad, pur_h5ad = run_mtx_to_h5ad(
        sample_id=sample_id,
        run_id=RUN_ID,
        unpurified_mtx_dir=unp_dir,
        purified_mtx_dir=pur_dir,
        output_root=output_root,
        h5ad_compression="gzip",
        force_rerun=False,
    )

    assert unp_h5ad.exists()
    assert pur_h5ad.exists()
    # Purified now lives under spatial_adata/ with the new basename.
    assert pur_h5ad.name == f"{sample_id}_proseg_purified.h5ad"
    assert pur_h5ad.parent.name == "spatial_adata"
    # Unpurified is intermediate.
    assert unp_h5ad.parent.name == "adata"
    assert "intermediate" in unp_h5ad.parts

    import anndata as ad
    # Unpurified: everything from the metadata.csv lands on obs
    # (unpurified is not stripped — only purified is).
    unp_adata = ad.read_h5ad(unp_h5ad)
    assert unp_adata.n_obs == n_cells
    assert unp_adata.n_vars == n_genes
    for col in (
        "first_type", "purification_status", "w1_larger_w2",
        "same_class", "nCount_Proseg",
        # Foreign columns preserved on unpurified — they get stripped
        # from purified only.
        "cell_area", "original_cell_id",
    ):
        assert col in unp_adata.obs.columns
    assert "spatial" in unp_adata.obsm
    assert unp_adata.obsm["spatial"].shape == (n_cells, 2)

    # Purified: SPLIT-native columns only. mtx_to_h5ad applies the
    # allow-list from _internal/split_native.py before writing.
    pur_adata = ad.read_h5ad(pur_h5ad)
    assert pur_adata.n_obs == n_cells
    assert pur_adata.n_vars == n_genes
    for col in (
        "first_type", "purification_status", "w1_larger_w2",
        "same_class", "nCount_Proseg",
    ):
        assert col in pur_adata.obs.columns
    for foreign in ("cell_area", "original_cell_id", "x", "y"):
        assert foreign not in pur_adata.obs.columns, (
            f"purified.obs should not carry {foreign!r} — it's not "
            "SPLIT-native and must be dropped by mtx_to_h5ad."
        )
    assert "spatial" in pur_adata.obsm
    assert pur_adata.obsm["spatial"].shape == (n_cells, 2)


def _write_bundle_with_nan_obs(tmp_dir: Path, sample_id: str, variant: str,
                                n_cells: int, n_genes: int):
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix
    from scipy.io import mmwrite

    rng = np.random.default_rng(0)
    counts = rng.poisson(lam=2.0, size=(n_genes, n_cells)).astype("int32")
    counts_csr = csr_matrix(counts)

    stem = tmp_dir / f"{sample_id}_{variant}"
    with gzip.open(f"{stem}_counts.mtx.gz", "wb") as f:
        mmwrite(f, counts_csr)

    gene_names = [f"gene{i}" for i in range(n_genes)]
    with gzip.open(f"{stem}_features.tsv.gz", "wt") as f:
        w = csv.writer(f, delimiter="\t")
        for g in gene_names:
            w.writerow([g])

    barcodes = [f"cell{i}" for i in range(n_cells)]
    with gzip.open(f"{stem}_barcodes.tsv.gz", "wt") as f:
        w = csv.writer(f, delimiter="\t")
        for b in barcodes:
            w.writerow([b])

    half = n_cells // 2
    first_type = np.array(
        ["tumor"] * half + [np.nan] * (n_cells - half), dtype=object
    )
    first_class = np.array(
        [False] * half + [np.nan] * (n_cells - half), dtype=object
    )
    spot_class = np.array(
        ["singlet"] * half + [np.nan] * (n_cells - half), dtype=object
    )
    meta = pd.DataFrame(
        {
            "first_type": first_type,
            "first_class": first_class,
            "spot_class": spot_class,
            "x": rng.uniform(0, 1000, size=n_cells),
            "y": rng.uniform(0, 1000, size=n_cells),
        },
        index=barcodes,
    )
    meta.to_csv(f"{stem}_metadata.csv", index=True)

    coords = pd.DataFrame(
        {"x": meta["x"].to_numpy(), "y": meta["y"].to_numpy()},
        index=barcodes,
    )
    coords.to_csv(f"{stem}_spatial_coords.csv.gz", index=True, compression="gzip")


def test_mtx_to_h5ad_handles_nan_in_object_obs_cols(tmp_path):
    """Regression for TracyY123-nexus#21."""
    from rctd_split.stages.mtx_to_h5ad import run_mtx_to_h5ad

    sample_id = "MHNAN"
    n_cells, n_genes = 20, 10
    output_root = tmp_path / "runs"
    unp_dir, pur_dir = _prepare_bundle_dirs(output_root, sample_id)

    _write_bundle_with_nan_obs(unp_dir, sample_id, "unpurified", n_cells, n_genes)
    _write_bundle_with_nan_obs(pur_dir, sample_id, "purified", n_cells, n_genes)

    unp_h5ad, _pur_h5ad = run_mtx_to_h5ad(
        sample_id=sample_id,
        run_id=RUN_ID,
        unpurified_mtx_dir=unp_dir,
        purified_mtx_dir=pur_dir,
        output_root=output_root,
        h5ad_compression="gzip",
        force_rerun=False,
    )

    import anndata as ad
    adata = ad.read_h5ad(unp_h5ad)
    half = n_cells // 2
    first_type = adata.obs["first_type"].astype(str).tolist()
    assert first_type[:half] == ["tumor"] * half
    assert first_type[half:] == [""] * (n_cells - half)
    first_class = adata.obs["first_class"].astype(str).tolist()
    assert first_class[:half] == ["False"] * half
    assert first_class[half:] == [""] * (n_cells - half)


def _write_split_bundles(
    output_root: Path, sample_id: str,
    n_cells_unp: int, n_cells_pur: int, n_genes: int,
) -> tuple[Path, Path]:
    import numpy as np
    import pandas as pd
    from scipy.io import mmwrite
    from scipy.sparse import csr_matrix

    unp_dir, pur_dir = _prepare_bundle_dirs(output_root, sample_id)

    gene_names = [f"gene{i}" for i in range(n_genes)]
    unp_barcodes = [f"cell{i}" for i in range(n_cells_unp)]
    pur_barcodes = unp_barcodes[:n_cells_pur]

    def _write(bundle_dir: Path, variant: str, barcodes: list[str], seed: int):
        rng = np.random.default_rng(seed)
        counts = rng.poisson(
            lam=2.0, size=(n_genes, len(barcodes)),
        ).astype("int32")
        stem = bundle_dir / f"{sample_id}_{variant}"

        with gzip.open(f"{stem}_counts.mtx.gz", "wb") as f:
            mmwrite(f, csr_matrix(counts))
        with gzip.open(f"{stem}_features.tsv.gz", "wt") as f:
            w = csv.writer(f, delimiter="\t")
            for g in gene_names:
                w.writerow([g])
        with gzip.open(f"{stem}_barcodes.tsv.gz", "wt") as f:
            w = csv.writer(f, delimiter="\t")
            for b in barcodes:
                w.writerow([b])
        meta = pd.DataFrame(
            {
                "first_type": rng.choice(["tumor", "Fibroblast"], size=len(barcodes)),
                "x": rng.uniform(0, 1000, size=len(barcodes)),
                "y": rng.uniform(0, 1000, size=len(barcodes)),
            },
            index=barcodes,
        )
        meta.to_csv(f"{stem}_metadata.csv", index=True)
        coords = pd.DataFrame(
            {"x": meta["x"].to_numpy(), "y": meta["y"].to_numpy()},
            index=barcodes,
        )
        coords.to_csv(f"{stem}_spatial_coords.csv.gz",
                      index=True, compression="gzip")

    _write(unp_dir, "unpurified", unp_barcodes, seed=0)
    _write(pur_dir, "purified",   pur_barcodes, seed=1)
    return unp_dir, pur_dir


def test_purified_adata_carries_unpurified_counts_layer(tmp_path):
    from rctd_split.stages.mtx_to_h5ad import run_mtx_to_h5ad

    sample_id = "MHLAYER"
    n_cells_unp, n_cells_pur, n_genes = 30, 12, 20
    output_root = tmp_path / "runs"
    unp_dir, pur_dir = _write_split_bundles(
        output_root, sample_id, n_cells_unp, n_cells_pur, n_genes,
    )

    _, pur_h5ad = run_mtx_to_h5ad(
        sample_id=sample_id,
        run_id=RUN_ID,
        unpurified_mtx_dir=unp_dir,
        purified_mtx_dir=pur_dir,
        output_root=output_root,
        h5ad_compression="gzip",
        force_rerun=False,
    )

    import numpy as np
    import anndata as ad
    adata = ad.read_h5ad(pur_h5ad)
    assert adata.n_obs == n_cells_pur
    assert adata.n_vars == n_genes
    assert "unpurified_counts" in adata.layers
    layer = adata.layers["unpurified_counts"]
    x = adata.X
    layer_shape = layer.shape if hasattr(layer, "shape") else np.asarray(layer).shape
    assert layer_shape == x.shape
    x_arr = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
    layer_arr = (
        layer.toarray() if hasattr(layer, "toarray") else np.asarray(layer)
    )
    assert layer_arr.sum() > 0
    assert not np.array_equal(x_arr, layer_arr)


def test_unpurified_adata_has_no_layer(tmp_path):
    from rctd_split.stages.mtx_to_h5ad import run_mtx_to_h5ad

    sample_id = "MHLAYER2"
    output_root = tmp_path / "runs"
    unp_dir, pur_dir = _write_split_bundles(output_root, sample_id, 20, 8, 15)

    unp_h5ad, _ = run_mtx_to_h5ad(
        sample_id=sample_id,
        run_id=RUN_ID,
        unpurified_mtx_dir=unp_dir,
        purified_mtx_dir=pur_dir,
        output_root=output_root,
        h5ad_compression="gzip",
        force_rerun=False,
    )

    import anndata as ad
    adata = ad.read_h5ad(unp_h5ad)
    assert "unpurified_counts" not in adata.layers


def test_mtx_to_h5ad_sentinel_skip(tmp_path):
    from rctd_split.stages.mtx_to_h5ad import run_mtx_to_h5ad

    sample_id = "MHSKIP"
    output_root = tmp_path / "runs"
    unp_dir, pur_dir = _prepare_bundle_dirs(output_root, sample_id)
    _write_bundle(unp_dir, sample_id, "unpurified", 10, 5)
    _write_bundle(pur_dir, sample_id, "purified", 10, 5)

    unp1, pur1 = run_mtx_to_h5ad(
        sample_id=sample_id,
        run_id=RUN_ID,
        unpurified_mtx_dir=unp_dir,
        purified_mtx_dir=pur_dir,
        output_root=output_root,
        h5ad_compression="gzip",
        force_rerun=False,
    )
    m1_unp = unp1.stat().st_mtime_ns
    m1_pur = pur1.stat().st_mtime_ns

    unp2, pur2 = run_mtx_to_h5ad(
        sample_id=sample_id,
        run_id=RUN_ID,
        unpurified_mtx_dir=unp_dir,
        purified_mtx_dir=pur_dir,
        output_root=output_root,
        h5ad_compression="gzip",
        force_rerun=False,
    )
    assert unp2.stat().st_mtime_ns == m1_unp
    assert pur2.stat().st_mtime_ns == m1_pur
