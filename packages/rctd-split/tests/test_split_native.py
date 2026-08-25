"""Tests for the SPLIT-native obs allow-list.

Two layers:

  * ``rctd_split._internal.split_native`` — the allow-list itself
    (exact + prefix rules) and :func:`filter_purified_obs`.
  * ``rctd-split strip-purified-obs`` CLI subcommand — the retroactive
    cleanup path for Tracy's existing MH8_2 purified.h5ad
    (settylab/TracyY123-nexus#26 comment 5274767383).
"""
from __future__ import annotations

from pathlib import Path

import pytest


pytest.importorskip("anndata")


def test_is_split_native_covers_seurat_managed_prefixes():
    from rctd_split._internal.split_native import is_split_native_obs_col

    assert is_split_native_obs_col("orig.ident")
    assert is_split_native_obs_col("nCount_Proseg")
    assert is_split_native_obs_col("nCount_Xenium")
    assert is_split_native_obs_col("nFeature_Proseg")


def test_is_split_native_covers_split_result_cols():
    from rctd_split._internal.split_native import is_split_native_obs_col

    for col in (
        "first_type", "second_type", "spot_class",
        "purification_status", "w1_larger_w2", "same_class",
    ):
        assert is_split_native_obs_col(col), col


def test_is_split_native_covers_postprocess_additions():
    from rctd_split._internal.split_native import is_split_native_obs_col

    # sc.pp.filter_cells side-effect + leiden output.
    assert is_split_native_obs_col("n_counts")
    assert is_split_native_obs_col("leiden_0.5")
    assert is_split_native_obs_col("leiden_1")


def test_is_split_native_rejects_raw_side_leaks():
    from rctd_split._internal.split_native import is_split_native_obs_col

    # Tracy's motivating case: raw-side per-cell total counts leaked
    # onto purified via preserve_meta_from_unpurified.
    for col in (
        "total_counts", "cell_area",
        "original_cell_id", "x", "y", "qc_filtered",
        "xenium_id_match", "xenium_cell_id_nn",
    ):
        assert not is_split_native_obs_col(col), col


def test_is_split_native_keeps_centroid_cols_for_celltype_writeback():
    """centroid_x / centroid_y are xenium-preprocess-derived (not SPLIT-native
    proper), but the Stage-D pipeline's celltype_writeback reads them
    off purified.obs. Kept on the allow-list so a fresh pipeline run
    doesn't fail at Stage 8."""
    from rctd_split._internal.split_native import is_split_native_obs_col

    assert is_split_native_obs_col("centroid_x")
    assert is_split_native_obs_col("centroid_y")


def test_filter_purified_obs_drops_only_foreign(tmp_path: Path):
    import anndata as ad
    import numpy as np
    import pandas as pd

    from rctd_split._internal.split_native import filter_purified_obs

    obs = pd.DataFrame({
        "first_type": pd.Categorical(["tumor", "stroma"]),
        "purification_status": pd.Categorical(["purified", "purified"]),
        "nCount_Proseg": [10.0, 20.0],
        "leiden_0.5": pd.Categorical(["0", "1"]),
        # Foreign — must be dropped.
        "total_counts": [11.0, 22.0],
        "cell_area": [100.0, 200.0],
    }, index=["C0", "C1"])
    adata = ad.AnnData(
        X=np.ones((2, 3), dtype="float32"),
        obs=obs,
        var=pd.DataFrame(index=["g0", "g1", "g2"]),
    )
    dropped = filter_purified_obs(adata)

    assert dropped == ["cell_area", "total_counts"]
    assert set(adata.obs.columns) == {
        "first_type", "purification_status", "nCount_Proseg", "leiden_0.5",
    }


def test_strip_purified_obs_subcommand_backs_up_and_strips(tmp_path: Path):
    """`rctd-split strip-purified-obs --purified-h5ad ...` writes a
    timestamped backup then rewrites the file with only SPLIT-native
    obs cols."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    from rctd_split._internal.layout import atomic_write_h5ad
    from rctd_split.cli import main

    obs = pd.DataFrame({
        "first_type": pd.Categorical(["tumor", "stroma", "tumor"]),
        "purification_status": pd.Categorical(
            ["purified", "purified", "flagged"],
        ),
        "nCount_Proseg": [10.0, 20.0, 30.0],
        # Foreign columns Tracy called out.
        "total_counts": [11.0, 22.0, 33.0],
        "cell_area": [100.0, 200.0, 300.0],
    }, index=["C0", "C1", "C2"])
    adata = ad.AnnData(
        X=np.ones((3, 3), dtype="float32"),
        obs=obs,
        var=pd.DataFrame(index=["g0", "g1", "g2"]),
    )
    p = tmp_path / "MHTEST_proseg_purified.h5ad"
    atomic_write_h5ad(adata, p, compression="gzip")

    rc = main(["strip-purified-obs", "--purified-h5ad", str(p)])
    assert rc == 0

    # Backup created next to the file.
    backups = sorted(tmp_path.glob("MHTEST_proseg_purified.h5ad.bak-*"))
    assert len(backups) == 1

    # Rewritten in place with SPLIT-native cols only.
    stripped = ad.read_h5ad(p)
    assert set(stripped.obs.columns) == {
        "first_type", "purification_status", "nCount_Proseg",
    }
    assert "total_counts" not in stripped.obs.columns
    assert "cell_area" not in stripped.obs.columns

    # Backup still has the foreign columns.
    original = ad.read_h5ad(backups[0])
    assert "total_counts" in original.obs.columns
    assert "cell_area" in original.obs.columns


def test_strip_purified_obs_subcommand_dry_run_makes_no_changes(tmp_path: Path):
    import anndata as ad
    import numpy as np
    import pandas as pd

    from rctd_split._internal.layout import atomic_write_h5ad
    from rctd_split.cli import main

    obs = pd.DataFrame({
        "first_type": pd.Categorical(["tumor"]),
        "total_counts": [42.0],
    }, index=["C0"])
    adata = ad.AnnData(
        X=np.ones((1, 2), dtype="float32"),
        obs=obs,
        var=pd.DataFrame(index=["g0", "g1"]),
    )
    p = tmp_path / "MHTEST_proseg_purified.h5ad"
    atomic_write_h5ad(adata, p, compression="gzip")
    mtime_before = p.stat().st_mtime_ns

    rc = main([
        "strip-purified-obs", "--purified-h5ad", str(p), "--dry-run",
    ])
    assert rc == 0
    assert p.stat().st_mtime_ns == mtime_before
    # Backup NOT written in dry-run mode.
    assert not list(tmp_path.glob("*.bak-*"))
    # File still carries the foreign column.
    still = ad.read_h5ad(p)
    assert "total_counts" in still.obs.columns


def test_strip_purified_obs_subcommand_idempotent_on_clean_file(tmp_path: Path):
    """Running on an already-clean file is a no-op — no backup and no
    rewrite."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    from rctd_split._internal.layout import atomic_write_h5ad
    from rctd_split.cli import main

    obs = pd.DataFrame({
        "first_type": pd.Categorical(["tumor"]),
        "nCount_Proseg": [10.0],
    }, index=["C0"])
    adata = ad.AnnData(
        X=np.ones((1, 2), dtype="float32"),
        obs=obs,
        var=pd.DataFrame(index=["g0", "g1"]),
    )
    p = tmp_path / "MHCLEAN_proseg_purified.h5ad"
    atomic_write_h5ad(adata, p, compression="gzip")
    mtime_before = p.stat().st_mtime_ns

    rc = main(["strip-purified-obs", "--purified-h5ad", str(p)])
    assert rc == 0
    # No backup, no rewrite.
    assert not list(tmp_path.glob("*.bak-*"))
    assert p.stat().st_mtime_ns == mtime_before
