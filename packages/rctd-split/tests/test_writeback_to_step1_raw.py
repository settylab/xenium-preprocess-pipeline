"""Tests for the writeback_to_step1_raw stage.

Fabricate:
  * a 5-cell step-1 raw h5ad with .obs['qc_filtered'] (3 True, 2 False)
  * a filter_status.csv with rows for exactly the 3 qc_filtered=True cells
  * a step-4 unpurified.h5ad also carrying those 3 cells with rich obs

Verify:
  1. filter-status cols land on all 5 raw cells (False on the 2 missing).
  2. exclude_obs_cols is honoured; folded columns arrive on raw.
  3. Idempotence: rerunning does not double-write or grow the schema.
  4. Fail-loud invariant: mark one qc_filtered=True cell absent from
     filter_status.csv and the stage raises + emits a summary CSV.
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("anndata")
pytest.importorskip("scipy")


def _write_raw(path: Path, cell_ids, qc_flags, extra_obs=None):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    n = len(cell_ids)
    obs_data = {"qc_filtered": [bool(f) for f in qc_flags]}
    if extra_obs:
        obs_data.update(extra_obs)
    obs = pd.DataFrame(obs_data, index=pd.Index(cell_ids, name="cell_id"))
    var = pd.DataFrame(index=["g0", "g1"])
    X = csr_matrix(np.random.default_rng(0).poisson(2.0, size=(n, 2)).astype("float32"))
    path.parent.mkdir(parents=True, exist_ok=True)
    from rctd_split._internal.layout import atomic_write_h5ad
    atomic_write_h5ad(ad.AnnData(X=X, obs=obs, var=var), path, compression="gzip")


def _write_unpurified(path: Path, cell_ids, extra_obs):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    n = len(cell_ids)
    obs = pd.DataFrame(extra_obs, index=pd.Index(cell_ids, name="cell_id"))
    var = pd.DataFrame(index=["g0", "g1"])
    X = csr_matrix(np.random.default_rng(1).poisson(2.0, size=(n, 2)).astype("float32"))
    path.parent.mkdir(parents=True, exist_ok=True)
    from rctd_split._internal.layout import atomic_write_h5ad
    atomic_write_h5ad(ad.AnnData(X=X, obs=obs, var=var), path, compression="gzip")


def _write_purified(path: Path, cell_ids):
    """Fabricate a post-postprocess ``proseg_purified.h5ad``: only the
    cells that survived ``sc.pp.filter_cells(min_counts=...)``. Minimal
    schema — only ``.obs_names`` is read by the writeback stage."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    n = len(cell_ids)
    obs = pd.DataFrame(index=pd.Index(cell_ids, name="cell_id"))
    var = pd.DataFrame(index=["g0", "g1"])
    X = csr_matrix(np.random.default_rng(2).poisson(5.0, size=(n, 2)).astype("float32"))
    path.parent.mkdir(parents=True, exist_ok=True)
    from rctd_split._internal.layout import atomic_write_h5ad
    atomic_write_h5ad(ad.AnnData(X=X, obs=obs, var=var), path, compression="gzip")


def _write_filter_status(path: Path, cell_ids, cols):
    import pandas as pd

    df = pd.DataFrame(cols, index=pd.Index(cell_ids, name="cell_id"))
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=True)


def _setup(tmp_path: Path):
    from rctd_split._internal.layout import (
        intermediate_path,
        spatial_adata_path,
    )

    sample_id = "MHTEST"
    run_id = "42"
    output_root = tmp_path / "runs"

    # 5 raw cells: C0..C2 passed qc (in filter_status.csv), C3..C4
    # failed qc (absent from filter_status.csv).
    raw_ids = ["C0", "C1", "C2", "C3", "C4"]
    qc_flags = [True, True, True, False, False]
    raw_p = spatial_adata_path(output_root, sample_id, run_id, "proseg_raw")
    _write_raw(raw_p, raw_ids, qc_flags,
               extra_obs={"n_counts": [50, 60, 70, 5, 3]})

    # filter_status: only the qc-passed cells. C0,C1 passed everything;
    # C2 was filtered by purification.
    fs_ids = ["C0", "C1", "C2"]
    fs_cols = {
        "passed_rctd": [True, True, True],
        "passed_split_purify": [True, True, False],
        "filtered_by_purification": [False, False, True],
    }
    fs_p = intermediate_path(output_root, sample_id, run_id, "filter_status_csv")
    _write_filter_status(fs_p, fs_ids, fs_cols)

    # unpurified.h5ad: same 3 cells, with RCTD-derived cols.
    unp_p = intermediate_path(output_root, sample_id, run_id, "unpurified_h5ad")
    unp_p.parent.mkdir(parents=True, exist_ok=True)
    _write_unpurified(unp_p, fs_ids, extra_obs={
        "spot_class": ["singlet", "singlet", "doublet_certain"],
        "first_type": ["tumor", "stroma", "immune"],
        "second_type": ["", "", "tumor"],
        "weight_first_type": [0.9, 0.8, 0.6],
    })

    # purified.h5ad: cells that survived SPLIT::purify AND the postprocess
    # min_counts filter. C0, C1 survive both; C2 was filtered_by_purification
    # (and would fail min_counts anyway); C3, C4 never made it in.
    pur_p = spatial_adata_path(output_root, sample_id, run_id, "proseg_purified")
    _write_purified(pur_p, ["C0", "C1"])

    return output_root, sample_id, run_id, raw_p, fs_p, unp_p, pur_p


def test_writeback_folds_filter_status_and_obs(tmp_path):
    from rctd_split.stages.writeback_to_step1_raw import run_writeback_to_step1_raw

    output_root, sample_id, run_id, raw_p, _, _, _ = _setup(tmp_path)
    run_writeback_to_step1_raw(
        sample_id=sample_id,
        run_id=run_id,
        output_root=output_root,
        raw_h5ad=None,
        qc_filtered_col="qc_filtered",
        exclude_obs_cols=[],
        h5ad_compression="gzip",
        force_rerun=False,
    )

    import anndata as ad
    raw = ad.read_h5ad(raw_p)
    # Filter-status cols present on all 5 raw cells.
    for col in ("passed_rctd", "passed_split_purify", "filtered_by_purification"):
        assert col in raw.obs.columns
    # C0,C1 all-True; C2 mixed; C3,C4 all-False (missing from fs).
    expected = {
        "C0": (True, True, False),
        "C1": (True, True, False),
        "C2": (True, False, True),
        "C3": (False, False, False),
        "C4": (False, False, False),
    }
    for cid, (pr, ps, pf) in expected.items():
        assert bool(raw.obs.loc[cid, "passed_rctd"]) is pr, cid
        assert bool(raw.obs.loc[cid, "passed_split_purify"]) is ps, cid
        assert bool(raw.obs.loc[cid, "filtered_by_purification"]) is pf, cid
    # Unpurified obs cols folded onto raw (exclude_obs_cols=[] here).
    for col in ("spot_class", "first_type", "second_type", "weight_first_type"):
        assert col in raw.obs.columns
    # Present on the 3 unpurified cells, empty/NaN on C3/C4.
    assert raw.obs.loc["C0", "first_type"] == "tumor"
    assert raw.obs.loc["C2", "first_type"] == "immune"
    assert raw.obs.loc["C3", "first_type"] == ""
    # passed_purification: True iff raw cell is in purified.obs_names
    # (C0, C1). C2 was filtered_by_purification; C3, C4 never made it in.
    assert "passed_purification" in raw.obs.columns
    assert bool(raw.obs.loc["C0", "passed_purification"]) is True
    assert bool(raw.obs.loc["C1", "passed_purification"]) is True
    assert bool(raw.obs.loc["C2", "passed_purification"]) is False
    assert bool(raw.obs.loc["C3", "passed_purification"]) is False
    assert bool(raw.obs.loc["C4", "passed_purification"]) is False


def test_writeback_respects_exclude_obs_cols(tmp_path):
    from rctd_split.stages.writeback_to_step1_raw import run_writeback_to_step1_raw

    output_root, sample_id, run_id, raw_p, _, _, _ = _setup(tmp_path)
    run_writeback_to_step1_raw(
        sample_id=sample_id,
        run_id=run_id,
        output_root=output_root,
        raw_h5ad=None,
        qc_filtered_col="qc_filtered",
        exclude_obs_cols=["spot_class", "second_type"],
        h5ad_compression="gzip",
        force_rerun=False,
    )
    import anndata as ad
    raw = ad.read_h5ad(raw_p)
    assert "spot_class" not in raw.obs.columns
    assert "second_type" not in raw.obs.columns
    assert "first_type" in raw.obs.columns
    assert "weight_first_type" in raw.obs.columns


def test_writeback_is_idempotent(tmp_path):
    """Running the stage twice must not grow raw.obs schema or dims,
    and must produce byte-identical column values."""
    from rctd_split.stages.writeback_to_step1_raw import run_writeback_to_step1_raw

    output_root, sample_id, run_id, raw_p, _, _, _ = _setup(tmp_path)
    kwargs = dict(
        sample_id=sample_id, run_id=run_id, output_root=output_root,
        raw_h5ad=None, qc_filtered_col="qc_filtered",
        exclude_obs_cols=[], h5ad_compression="gzip",
    )
    run_writeback_to_step1_raw(force_rerun=False, **kwargs)
    import anndata as ad
    r1 = ad.read_h5ad(raw_p)
    n_obs = r1.n_obs
    cols_1 = list(r1.obs.columns)
    ft_1 = list(r1.obs["first_type"])
    pr_1 = [bool(x) for x in r1.obs["passed_rctd"]]

    run_writeback_to_step1_raw(force_rerun=True, **kwargs)
    r2 = ad.read_h5ad(raw_p)
    assert r2.n_obs == n_obs
    assert list(r2.obs.columns) == cols_1
    assert list(r2.obs["first_type"]) == ft_1
    assert [bool(x) for x in r2.obs["passed_rctd"]] == pr_1


def test_writeback_preserves_source_h5ad_shape_and_x(tmp_path):
    """Tracy's "source files otherwise unchanged" invariant
    (`settylab/TracyY123-nexus#26` comment 5260916505): the writeback
    only ADDS obs columns — it must not touch `.X`, `.var`, `.obs_names`,
    or the pre-existing obs columns.
    """
    import anndata as ad
    import numpy as np
    from rctd_split.stages.writeback_to_step1_raw import run_writeback_to_step1_raw

    output_root, sample_id, run_id, raw_p, _, _, _ = _setup(tmp_path)

    before = ad.read_h5ad(raw_p)
    x_before = before.X.toarray().copy()
    var_before = list(before.var.index)
    obs_names_before = list(before.obs_names)
    n_counts_before = list(before.obs["n_counts"])
    qc_before = [bool(v) for v in before.obs["qc_filtered"]]
    pre_existing_cols = set(before.obs.columns)

    run_writeback_to_step1_raw(
        sample_id=sample_id, run_id=run_id, output_root=output_root,
        raw_h5ad=None, qc_filtered_col="qc_filtered",
        exclude_obs_cols=[], h5ad_compression="gzip",
        force_rerun=False,
    )

    after = ad.read_h5ad(raw_p)
    assert after.n_obs == before.n_obs
    assert after.n_vars == before.n_vars
    assert list(after.obs_names) == obs_names_before
    assert list(after.var.index) == var_before
    np.testing.assert_array_equal(after.X.toarray(), x_before)
    assert list(after.obs["n_counts"]) == n_counts_before
    assert [bool(v) for v in after.obs["qc_filtered"]] == qc_before
    # Every pre-existing column still there; new columns are additions.
    assert pre_existing_cols.issubset(set(after.obs.columns))


def test_writeback_fails_loud_on_unaccounted_cells(tmp_path):
    """A raw cell that is qc_filtered=True but absent from
    filter_status.csv MUST fail loud + emit a summary CSV."""
    from rctd_split._internal.layout import intermediate_path
    from rctd_split.stages.writeback_to_step1_raw import run_writeback_to_step1_raw

    output_root, sample_id, run_id, raw_p, fs_p, _, _ = _setup(tmp_path)
    # Mark C3 as qc_filtered=True — but it's NOT in filter_status.csv.
    import anndata as ad
    raw = ad.read_h5ad(raw_p)
    raw.obs.loc["C3", "qc_filtered"] = True
    from rctd_split._internal.layout import atomic_write_h5ad
    atomic_write_h5ad(raw, raw_p, compression="gzip")

    with pytest.raises(SystemExit) as exc:
        run_writeback_to_step1_raw(
            sample_id=sample_id,
            run_id=run_id,
            output_root=output_root,
            raw_h5ad=None,
            qc_filtered_col="qc_filtered",
            exclude_obs_cols=[],
            h5ad_compression="gzip",
            force_rerun=False,
        )
    assert "unaccounted" in str(exc.value).lower() or "not in filter_status" in str(exc.value).lower()

    # The summary CSV must exist and mention C3.
    summary = intermediate_path(
        output_root, sample_id, run_id, "writeback_unaccounted",
    )
    assert summary.exists()
    body = summary.read_text()
    assert "C3" in body


def test_writeback_fails_loud_when_qc_col_missing(tmp_path):
    from rctd_split.stages.writeback_to_step1_raw import run_writeback_to_step1_raw

    output_root, sample_id, run_id, raw_p, _, _, _ = _setup(tmp_path)
    import anndata as ad
    raw = ad.read_h5ad(raw_p)
    del raw.obs["qc_filtered"]
    from rctd_split._internal.layout import atomic_write_h5ad
    atomic_write_h5ad(raw, raw_p, compression="gzip")

    with pytest.raises(SystemExit) as exc:
        run_writeback_to_step1_raw(
            sample_id=sample_id,
            run_id=run_id,
            output_root=output_root,
            raw_h5ad=None,
            qc_filtered_col="qc_filtered",
            exclude_obs_cols=[],
            h5ad_compression="gzip",
            force_rerun=False,
        )
    assert "qc_filter" in str(exc.value)


def test_writeback_bool_col_with_missing_cells_is_writeable(tmp_path):
    """Regression: a bool column on unpurified.obs, folded onto raw.obs
    where some cells are missing, must land as a proper bool dtype
    (NA → False) — not an object array of `True` / `False` / `pd.NA`
    that anndata's writer routes through `write_vlen_string_array` and
    fails with "Can't implicitly convert non-string objects to strings"
    (settylab/TracyY123-nexus#26 comment 5271431123)."""
    import anndata as ad
    from rctd_split.stages.writeback_to_step1_raw import run_writeback_to_step1_raw

    output_root, sample_id, run_id, raw_p, _, unp_p, _ = _setup(tmp_path)

    # Add a bool col on unpurified (present only on C0,C1,C2 — C3,C4
    # are absent from unpurified, so reindex will introduce NAs).
    unp = ad.read_h5ad(unp_p)
    unp.obs["some_bool_flag"] = [True, True, False]
    from rctd_split._internal.layout import atomic_write_h5ad
    atomic_write_h5ad(unp, unp_p, compression="gzip")

    # Must not raise. Prior to the fix this raised
    # "TypeError: Can't implicitly convert non-string objects to strings".
    run_writeback_to_step1_raw(
        sample_id=sample_id,
        run_id=run_id,
        output_root=output_root,
        raw_h5ad=None,
        qc_filtered_col="qc_filtered",
        exclude_obs_cols=[],
        h5ad_compression="gzip",
        force_rerun=False,
    )

    raw = ad.read_h5ad(raw_p)
    assert "some_bool_flag" in raw.obs.columns
    # Plain bool dtype (or nullable Boolean after read) — not object.
    assert raw.obs["some_bool_flag"].dtype != object, (
        f"expected bool/BooleanDtype, got {raw.obs['some_bool_flag'].dtype!r}"
    )
    # C0,C1 True, C2 False; C3,C4 default False (absent from unpurified).
    assert bool(raw.obs.loc["C0", "some_bool_flag"]) is True
    assert bool(raw.obs.loc["C1", "some_bool_flag"]) is True
    assert bool(raw.obs.loc["C2", "some_bool_flag"]) is False
    assert bool(raw.obs.loc["C3", "some_bool_flag"]) is False
    assert bool(raw.obs.loc["C4", "some_bool_flag"]) is False


def test_passed_purification_sum_equals_purified_n_obs(tmp_path):
    """Invariant (settylab/TracyY123-nexus#26 comment 5272585455):
    ``raw.obs['passed_purification'].sum() == purified.n_obs``. Every
    cell in ``proseg_purified.h5ad`` must map to exactly one raw cell —
    if the invariant fails, either purified has cells not in raw
    (impossible: purified is a filtered subset of raw's ids) or the
    membership test is broken."""
    import anndata as ad
    from rctd_split.stages.writeback_to_step1_raw import run_writeback_to_step1_raw

    output_root, sample_id, run_id, raw_p, _, _, pur_p = _setup(tmp_path)
    run_writeback_to_step1_raw(
        sample_id=sample_id,
        run_id=run_id,
        output_root=output_root,
        raw_h5ad=None,
        qc_filtered_col="qc_filtered",
        exclude_obs_cols=[],
        h5ad_compression="gzip",
        force_rerun=False,
    )

    raw = ad.read_h5ad(raw_p)
    purified = ad.read_h5ad(pur_p)
    assert "passed_purification" in raw.obs.columns
    assert raw.obs["passed_purification"].dtype != object
    assert int(raw.obs["passed_purification"].sum()) == purified.n_obs
