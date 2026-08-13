"""Unit tests for the filter_status stage.

Fabricate tiny intermediate unpurified.h5ad (10 cells, 6 RCTD-called +
4 with first_type=NaN) and a matching purified.h5ad in spatial_adata/
(5 cells, all a subset of the 6 RCTD-called); run the stage; verify:

  1. Sidecar CSV has 10 rows with correct boolean columns.
  2. Augmented unpurified.h5ad has the three columns AND matches the
     sidecar exactly.
  3. write_inplace=false leaves unpurified.h5ad unmodified.
  4. Idempotent — force_rerun=True re-writes without corrupting data.

Skipped cleanly if anndata / scanpy / scipy aren't importable.
"""
from __future__ import annotations

from pathlib import Path

import pytest


pytest.importorskip("anndata")
pytest.importorskip("scipy")
pytest.importorskip("pandas")


RUN_ID = "42"


def _write_unpurified(path: Path, rctd_called_ids, rctd_rejected_ids):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    cell_ids = rctd_called_ids + rctd_rejected_ids
    n = len(cell_ids)
    n_genes = 5

    rng = np.random.default_rng(0)
    X = csr_matrix(rng.poisson(2.0, size=(n, n_genes)).astype("float32"))

    first_type = (
        ["tumor"] * len(rctd_called_ids) + [np.nan] * len(rctd_rejected_ids)
    )
    obs = pd.DataFrame(
        {
            "first_type": first_type,
            "nCount_Proseg": rng.integers(10, 1000, size=n),
            "original_cell_id": [f"orig_{c}" for c in cell_ids],
        },
        index=pd.Index(cell_ids, name="cell_id"),
    )
    var = pd.DataFrame(index=[f"gene{i}" for i in range(n_genes)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path, compression="gzip")


def _write_purified(path: Path, purified_ids):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    n = len(purified_ids)
    n_genes = 5
    rng = np.random.default_rng(1)
    X = csr_matrix(rng.poisson(2.0, size=(n, n_genes)).astype("float32"))
    obs = pd.DataFrame(
        {"first_type": ["tumor"] * n},
        index=pd.Index(purified_ids, name="cell_id"),
    )
    var = pd.DataFrame(index=[f"gene{i}" for i in range(n_genes)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path, compression="gzip")


def _setup_scenario(tmp_path: Path, sample_id: str = "MHTEST"):
    from rctd_split._internal.layout import intermediate_path, spatial_adata_path

    output_root = tmp_path / "runs"

    called = [f"C{i}" for i in range(6)]
    rejected = [f"C{i}" for i in range(6, 10)]
    purified = [f"C{i}" for i in range(5)]

    unpurified = intermediate_path(
        output_root, sample_id, RUN_ID, "unpurified_h5ad",
    )
    purified_p = spatial_adata_path(
        output_root, sample_id, RUN_ID, "proseg_purified",
    )
    _write_unpurified(unpurified, called, rejected)
    _write_purified(purified_p, purified)

    sidecar = intermediate_path(
        output_root, sample_id, RUN_ID, "filter_status_csv",
    )

    expected = {}
    for i in range(10):
        cid = f"C{i}"
        p_rctd = i < 6
        p_split = i < 5
        p_filt = p_rctd and not p_split
        expected[cid] = (p_rctd, p_split, p_filt)
    return output_root, sidecar, unpurified, expected


def test_filter_status_writes_sidecar_and_augments_h5ad(tmp_path):
    from rctd_split.stages.filter_status import (
        DEFAULT_INPLACE_COLUMNS, run_filter_status,
    )

    sample_id = "MHTEST"
    output_root, sidecar, unpurified_h5ad, expected = _setup_scenario(
        tmp_path, sample_id
    )

    returned = run_filter_status(
        sample_id=sample_id,
        run_id=RUN_ID,
        output_root=output_root,
        write_inplace=True,
        inplace_columns=DEFAULT_INPLACE_COLUMNS,
        first_type_column="first_type",
        h5ad_compression="gzip",
        force_rerun=False,
    )
    assert returned == sidecar
    assert sidecar.exists()

    import pandas as pd
    df = pd.read_csv(sidecar, index_col=0)
    assert df.shape == (10, 3)
    assert list(df.columns) == list(DEFAULT_INPLACE_COLUMNS)
    for cid, (p_rctd, p_split, p_filt) in expected.items():
        row = df.loc[cid]
        assert bool(row["passed_rctd"]) is p_rctd, cid
        assert bool(row["passed_split_purify"]) is p_split, cid
        assert bool(row["filtered_by_purification"]) is p_filt, cid

    import anndata as ad
    adata = ad.read_h5ad(unpurified_h5ad)
    for col in DEFAULT_INPLACE_COLUMNS:
        assert col in adata.obs.columns
    for cid, (p_rctd, p_split, p_filt) in expected.items():
        obs_row = adata.obs.loc[cid]
        assert bool(obs_row["passed_rctd"]) is p_rctd
        assert bool(obs_row["passed_split_purify"]) is p_split
        assert bool(obs_row["filtered_by_purification"]) is p_filt

    for col in DEFAULT_INPLACE_COLUMNS:
        for cid in adata.obs_names:
            assert bool(adata.obs.loc[cid, col]) == bool(df.loc[cid, col])

    for col in ("first_type", "nCount_Proseg", "original_cell_id"):
        assert col in adata.obs.columns


def test_filter_status_write_inplace_false_leaves_h5ad_unchanged(tmp_path):
    from rctd_split.stages.filter_status import (
        DEFAULT_INPLACE_COLUMNS, run_filter_status,
    )

    sample_id = "MHTEST"
    output_root, sidecar, unpurified_h5ad, _ = _setup_scenario(
        tmp_path, sample_id
    )

    mtime_before = unpurified_h5ad.stat().st_mtime_ns
    size_before = unpurified_h5ad.stat().st_size

    run_filter_status(
        sample_id=sample_id,
        run_id=RUN_ID,
        output_root=output_root,
        write_inplace=False,
        inplace_columns=DEFAULT_INPLACE_COLUMNS,
        first_type_column="first_type",
        h5ad_compression="gzip",
        force_rerun=False,
    )
    assert sidecar.exists()

    assert unpurified_h5ad.stat().st_mtime_ns == mtime_before
    assert unpurified_h5ad.stat().st_size == size_before

    import anndata as ad
    adata = ad.read_h5ad(unpurified_h5ad)
    for col in DEFAULT_INPLACE_COLUMNS:
        assert col not in adata.obs.columns


def test_filter_status_force_rerun_overwrites(tmp_path):
    from rctd_split.stages.filter_status import (
        DEFAULT_INPLACE_COLUMNS, run_filter_status,
    )

    sample_id = "MHTEST"
    output_root, sidecar, unpurified_h5ad, expected = _setup_scenario(
        tmp_path, sample_id
    )

    run_filter_status(
        sample_id=sample_id,
        run_id=RUN_ID,
        output_root=output_root,
        write_inplace=True,
        inplace_columns=DEFAULT_INPLACE_COLUMNS,
        first_type_column="first_type",
        h5ad_compression="gzip",
        force_rerun=False,
    )
    run_filter_status(
        sample_id=sample_id,
        run_id=RUN_ID,
        output_root=output_root,
        write_inplace=True,
        inplace_columns=DEFAULT_INPLACE_COLUMNS,
        first_type_column="first_type",
        h5ad_compression="gzip",
        force_rerun=True,
    )
    import anndata as ad
    adata = ad.read_h5ad(unpurified_h5ad)
    for cid, (p_rctd, p_split, p_filt) in expected.items():
        obs_row = adata.obs.loc[cid]
        assert bool(obs_row["passed_rctd"]) is p_rctd
        assert bool(obs_row["passed_split_purify"]) is p_split
        assert bool(obs_row["filtered_by_purification"]) is p_filt


def test_filter_status_missing_first_type_raises(tmp_path):
    from rctd_split.stages.filter_status import (
        DEFAULT_INPLACE_COLUMNS, run_filter_status,
    )
    import anndata as ad

    sample_id = "MHTEST"
    output_root, _, unpurified_h5ad, _ = _setup_scenario(tmp_path, sample_id)

    adata = ad.read_h5ad(unpurified_h5ad)
    del adata.obs["first_type"]
    adata.write_h5ad(unpurified_h5ad, compression="gzip")

    with pytest.raises(SystemExit) as exc:
        run_filter_status(
            sample_id=sample_id,
            run_id=RUN_ID,
            output_root=output_root,
            write_inplace=True,
            inplace_columns=DEFAULT_INPLACE_COLUMNS,
            first_type_column="first_type",
            h5ad_compression="gzip",
            force_rerun=False,
        )
    assert "first_type" in str(exc.value)
