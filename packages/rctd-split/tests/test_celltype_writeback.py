"""Tests for the celltype_writeback stage.

Fabricate:
  * a 4-cell reference h5ad simulating <S>_proseg_purified.h5ad with
    .obs[first_type, centroid_x, centroid_y]
  * a 6-cell query h5ad simulating <S>_xenium_ranger.h5ad with
    .obsm['spatial']

Verify:
  1. NN direction: labels flow reference → query (query cells receive
     labels of their nearest reference cell).
  2. Every query cell gets a label under the default 'nearest_label'
     policy (no Unassigned).
  3. Idempotence: re-running overwrites, no dimensional growth.
  4. Reference col missing → fail loud.
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("anndata")
pytest.importorskip("sklearn")
pytest.importorskip("scipy")


def _write_reference(path: Path, sample_id: str = "MHTEST"):
    """Write a fake proseg_purified.h5ad with (x,y, first_type)."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    n_ref = 4
    xs = [0.0, 100.0, 0.0, 100.0]
    ys = [0.0, 0.0, 100.0, 100.0]
    types = ["tumor", "stroma", "immune", "endo"]
    obs = pd.DataFrame(
        {"first_type": types, "centroid_x": xs, "centroid_y": ys},
        index=pd.Index([f"REF{i}" for i in range(n_ref)], name="cell_id"),
    )
    var = pd.DataFrame(index=["g0", "g1"])
    X = csr_matrix(np.random.default_rng(0).poisson(2.0, size=(n_ref, 2)).astype("float32"))
    ad.AnnData(X=X, obs=obs, var=var).write_h5ad(path, compression="gzip")


def _write_query(path: Path, coords):
    """Write a fake xenium_ranger.h5ad with .obsm['spatial']."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    n = len(coords)
    obs = pd.DataFrame(index=pd.Index([f"Q{i}" for i in range(n)], name="cell_id"))
    var = pd.DataFrame(index=["g0", "g1"])
    X = csr_matrix(np.random.default_rng(1).poisson(2.0, size=(n, 2)).astype("float32"))
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm["spatial"] = np.asarray(coords, dtype=float)
    adata.write_h5ad(path, compression="gzip")


def _setup_run(tmp_path: Path):
    from rctd_split._internal.layout import spatial_adata_path

    sample_id = "MHTEST"
    run_id = "42"
    output_root = tmp_path / "runs"

    ref_p = spatial_adata_path(output_root, sample_id, run_id, "proseg_purified")
    query_p = spatial_adata_path(output_root, sample_id, run_id, "xenium_ranger")
    ref_p.parent.mkdir(parents=True, exist_ok=True)
    _write_reference(ref_p, sample_id)

    query_coords = [
        [0.1, 0.1],    # near REF0 (tumor)
        [99.9, 0.0],   # near REF1 (stroma)
        [0.0, 99.9],   # near REF2 (immune)
        [100.0, 100.0],# REF3 (endo)
        [50.0, 0.0],   # equidistant to REF0/REF1; sklearn picks the first-index tie
        [50.0, 50.0],  # equidistant to all four; sklearn picks lowest index
    ]
    _write_query(query_p, query_coords)
    return output_root, sample_id, run_id, ref_p, query_p


def test_celltype_writeback_maps_labels(tmp_path):
    from rctd_split.stages.celltype_writeback import run_celltype_writeback

    output_root, sample_id, run_id, _, query_p = _setup_run(tmp_path)

    run_celltype_writeback(
        sample_id=sample_id,
        run_id=run_id,
        output_root=output_root,
        xenium_ranger_h5ad=None,
        reference_label_col="first_type",
        reference_x_col="centroid_x",
        reference_y_col="centroid_y",
        k=1,
        algorithm="auto",
        metric="euclidean",
        tiebreak="min_dist",
        distance_threshold=None,
        unmatched_policy="nearest_label",
        celltype_col="celltype",
        distance_col="celltype_source_distance",
        h5ad_compression="gzip",
        force_rerun=False,
    )

    import anndata as ad
    query = ad.read_h5ad(query_p)
    assert "celltype" in query.obs.columns
    assert "celltype_source_distance" in query.obs.columns
    # Q0..Q3 should get tumor / stroma / immune / endo respectively.
    labels = list(query.obs["celltype"])
    assert labels[0] == "tumor"
    assert labels[1] == "stroma"
    assert labels[2] == "immune"
    assert labels[3] == "endo"
    # Q0 is right on REF0.
    d = query.obs["celltype_source_distance"].to_numpy()
    assert d[3] == 0.0  # exact match
    # Default policy: NO "Unassigned" labels.
    assert not any(lbl == "Unassigned" for lbl in labels)


def test_celltype_writeback_is_idempotent(tmp_path):
    from rctd_split.stages.celltype_writeback import run_celltype_writeback

    output_root, sample_id, run_id, _, query_p = _setup_run(tmp_path)

    common = dict(
        sample_id=sample_id, run_id=run_id, output_root=output_root,
        xenium_ranger_h5ad=None, reference_label_col="first_type",
        reference_x_col="centroid_x", reference_y_col="centroid_y",
        k=1, algorithm="auto", metric="euclidean", tiebreak="min_dist",
        distance_threshold=None, unmatched_policy="nearest_label",
        celltype_col="celltype", distance_col="celltype_source_distance",
        h5ad_compression="gzip",
    )
    run_celltype_writeback(force_rerun=False, **common)

    import anndata as ad
    q1 = ad.read_h5ad(query_p)
    n_obs = q1.n_obs
    cols_1 = list(q1.obs.columns)

    # Force re-run — must not double-write or grow dims.
    run_celltype_writeback(force_rerun=True, **common)
    q2 = ad.read_h5ad(query_p)
    assert q2.n_obs == n_obs
    assert list(q2.obs.columns) == cols_1
    # Same labels, byte-for-byte.
    assert list(q1.obs["celltype"]) == list(q2.obs["celltype"])


def test_celltype_writeback_preserves_source_h5ad_shape_and_x(tmp_path):
    """Tracy's "source files otherwise unchanged" invariant
    (`settylab/TracyY123-nexus#26` comment 5260916505): the writeback
    only ADDS `celltype` + `celltype_source_distance` to the query
    xenium_ranger h5ad — it must not touch `.X`, `.var`, `.obs_names`,
    `.obsm['spatial']`, or the reference (source) purified h5ad.
    """
    import anndata as ad
    import numpy as np
    from rctd_split.stages.celltype_writeback import run_celltype_writeback

    output_root, sample_id, run_id, ref_p, query_p = _setup_run(tmp_path)

    q_before = ad.read_h5ad(query_p)
    r_before = ad.read_h5ad(ref_p)
    qx_before = q_before.X.toarray().copy()
    qvar_before = list(q_before.var.index)
    qnames_before = list(q_before.obs_names)
    qspatial_before = np.asarray(q_before.obsm["spatial"]).copy()
    rx_before = r_before.X.toarray().copy()
    r_obs_before = r_before.obs.copy()

    run_celltype_writeback(
        sample_id=sample_id, run_id=run_id, output_root=output_root,
        xenium_ranger_h5ad=None, reference_label_col="first_type",
        reference_x_col="centroid_x", reference_y_col="centroid_y",
        k=1, algorithm="auto", metric="euclidean", tiebreak="min_dist",
        distance_threshold=None, unmatched_policy="nearest_label",
        celltype_col="celltype", distance_col="celltype_source_distance",
        h5ad_compression="gzip", force_rerun=False,
    )

    # Query (target): shape / X / obs_names / spatial / var unchanged;
    # only new obs cols added.
    q_after = ad.read_h5ad(query_p)
    assert q_after.n_obs == q_before.n_obs
    assert q_after.n_vars == q_before.n_vars
    assert list(q_after.obs_names) == qnames_before
    assert list(q_after.var.index) == qvar_before
    np.testing.assert_array_equal(q_after.X.toarray(), qx_before)
    np.testing.assert_array_equal(
        np.asarray(q_after.obsm["spatial"]), qspatial_before
    )
    added = set(q_after.obs.columns) - set(q_before.obs.columns)
    assert added == {"celltype", "celltype_source_distance"}

    # Reference (source): completely untouched.
    r_after = ad.read_h5ad(ref_p)
    np.testing.assert_array_equal(r_after.X.toarray(), rx_before)
    assert list(r_after.obs.columns) == list(r_obs_before.columns)
    assert r_after.n_obs == r_before.n_obs


def test_celltype_writeback_fails_on_missing_reference_col(tmp_path):
    from rctd_split.stages.celltype_writeback import run_celltype_writeback

    output_root, sample_id, run_id, ref_p, _ = _setup_run(tmp_path)

    # Strip reference_label_col.
    import anndata as ad
    ref = ad.read_h5ad(ref_p)
    del ref.obs["first_type"]
    ref.write_h5ad(ref_p, compression="gzip")

    with pytest.raises(SystemExit) as exc:
        run_celltype_writeback(
            sample_id=sample_id, run_id=run_id, output_root=output_root,
            xenium_ranger_h5ad=None, reference_label_col="first_type",
            reference_x_col="centroid_x", reference_y_col="centroid_y",
            k=1, algorithm="auto", metric="euclidean", tiebreak="min_dist",
            distance_threshold=None, unmatched_policy="nearest_label",
            celltype_col="celltype", distance_col="celltype_source_distance",
            h5ad_compression="gzip", force_rerun=False,
        )
    assert "first_type" in str(exc.value)


# ---------------------------------------------------------------------
# v0.2.9: RCTD results + filter provenance writeback on the
# xenium_ranger-only path (settylab/TracyY123-nexus#40).
# ---------------------------------------------------------------------

def test_celltype_writeback_v0_2_9_folds_rctd_results(tmp_path, monkeypatch):
    """When ``rctd_results_rds`` is supplied, celltype_writeback also
    writes the 12 rctd_* obs columns + ``in_rctd_input`` +
    ``in_purified`` booleans onto xenium_ranger.h5ad. Monkeypatches the
    R-side extractor since sandbox CI doesn't have spacexr/SPLIT; the
    schema the fake extractor returns is the exact shape the real one
    emits (see ``packages/rctd-split/src/rctd_split/r/rctd_results_to_csv.R``).
    """
    import anndata as ad
    import pandas as pd
    import numpy as np
    from rctd_split.stages.celltype_writeback import run_celltype_writeback

    output_root, sample_id, run_id, _, query_p = _setup_run(tmp_path)

    # Fake rctd_results covers 4 of the 6 query cells; 2 cells (Q4, Q5)
    # are NOT in the RCTD input set — they should get fillers ("" / NaN)
    # and ``in_rctd_input = False``.
    synth = pd.DataFrame(
        {
            "spot_class":          ["singlet", "doublet_certain",
                                    "doublet_uncertain", "reject"],
            "first_type":          ["tumor", "stroma", "immune", "endo"],
            "second_type":         ["", "immune", "", ""],
            "min_score":           [0.9, 0.85, 0.6, 0.3],
            "singlet_score":       [0.95, 0.4, 0.5, 0.1],
            "purification_status": ["purified", "raw", "raw", "rejected"],
        },
        index=pd.Index([f"Q{i}" for i in range(4)], name="cell_id"),
    )

    def _fake_extract(*, rctd_results_rds, rscript_bin, r_lib_paths,
                      tmp_csv=None):
        return synth

    # `_fold_rctd_results_onto_query` imports `extract_rctd_results_df`
    # LAZILY at function-body scope via
    # `from rctd_split._internal.rctd_results import extract_rctd_results_df`.
    # The lazy import resolves the symbol out of the
    # `rctd_split._internal.rctd_results` module at CALL time, so patching
    # THAT module replaces the function the lazy import will grab.
    # Prior attempt to also patch
    # `rctd_split.stages.celltype_writeback.extract_rctd_results_df`
    # errored at setup because the name isn't bound at the stage-module
    # level (monkeypatch.setattr raises AttributeError by default when the
    # attribute is missing, taking down the whole test before the real
    # patch ran). Skeptic pass on branch v0.2.9 (depth 1) caught this.
    import rctd_split._internal.rctd_results as _rr
    monkeypatch.setattr(_rr, "extract_rctd_results_df", _fake_extract)

    run_celltype_writeback(
        sample_id=sample_id, run_id=run_id, output_root=output_root,
        xenium_ranger_h5ad=None, reference_label_col="first_type",
        reference_x_col="centroid_x", reference_y_col="centroid_y",
        k=1, algorithm="auto", metric="euclidean", tiebreak="min_dist",
        distance_threshold=None, unmatched_policy="nearest_label",
        celltype_col="celltype", distance_col="celltype_source_distance",
        h5ad_compression="gzip", force_rerun=True,
        # v0.2.9 new args:
        rctd_results_rds=tmp_path / "fake_rctd_results.rds",  # bypassed
        rctd_obs_prefix="rctd_",
    )

    q = ad.read_h5ad(query_p)

    # All 6 rctd_* cols land, with the right fill behavior on Q4/Q5.
    for col in ("rctd_spot_class", "rctd_first_type", "rctd_second_type",
                "rctd_min_score", "rctd_singlet_score",
                "rctd_purification_status"):
        assert col in q.obs.columns, col

    # Cells Q0..Q3 carry real values; Q4/Q5 carry the fillers.
    assert list(q.obs["rctd_spot_class"].iloc[:4]) == [
        "singlet", "doublet_certain", "doublet_uncertain", "reject",
    ]
    assert list(q.obs["rctd_spot_class"].iloc[4:]) == ["", ""]
    assert list(q.obs["rctd_first_type"].iloc[:4]) == [
        "tumor", "stroma", "immune", "endo",
    ]
    # Numeric fill is NaN.
    assert np.isnan(q.obs["rctd_min_score"].iloc[4:]).all()

    # Provenance flags.
    # v0.2.9 Part 1 (#40 c6755625444): in_rctd_input semantics fix.
    # This test does NOT pass a barcodes path, so the code falls back
    # to `in_rctd_input == has_rctd_result` (= cells in rctd_df,
    # [T,T,T,T,F,F]) with a WARN log line. See the v0.2.9 Part 1 tests
    # below for the correctly-separated semantics with a barcodes path.
    assert list(q.obs["in_rctd_input"]) == [True, True, True, True, False, False]
    assert list(q.obs["has_rctd_result"]) == [True, True, True, True, False, False]
    # in_purified: the reference h5ad has REF0..REF3; none of the Q*
    # cell_ids match, so in_purified is False for every query cell.
    assert list(q.obs["in_purified"]) == [False, False, False, False, False, False]


# ---------------------------------------------------------------------
# v0.2.9 Part 1 — in_rctd_input semantic correction + Part-1 schema
# (settylab/TracyY123-nexus#40 c6755625444).
#
# The Part-1 fix splits what v0.2.9-prior called `in_rctd_input` into
# two distinct columns:
#   * `in_rctd_input`   — cell was SUBMITTED to RCTD (from split_prep
#                         barcodes.tsv.gz).
#   * `has_rctd_result` — cell has a row in rctd_results_df.
# And sources `rctd_purification_status` from the purified h5ad
# (populated by SPLIT::purify), NOT from the rctd_results_df (which
# doesn't carry it after only `SPLIT::run_post_process_RCTD`).
# ---------------------------------------------------------------------


def _write_barcodes_tsv_gz(path, barcodes):
    """Helper: write a 10x-style barcodes.tsv.gz (one cell id per line)."""
    import gzip
    with gzip.open(path, "wt") as f:
        for b in barcodes:
            f.write(f"{b}\n")


def _write_reference_with_purification_status(
    path, sample_id="MHTEST", sel_ids=("Q0", "Q1"),
    pstatus=("singlet", "doublet_certain"),
):
    """Like `_write_reference`, but rows are Q0/Q1 (so the Q-named
    query cells actually match) + carries `purification_status`.
    Simulates what a mtx_to_h5ad + postprocess-filtered purified h5ad
    looks like on an xenium_ranger-source run."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    n = len(sel_ids)
    obs = pd.DataFrame(
        {
            "first_type":           ["tumor"] * n,
            "centroid_x":           [0.0] * n,
            "centroid_y":           [0.0] * n,
            "purification_status":  list(pstatus),
        },
        index=pd.Index(list(sel_ids), name="cell_id"),
    )
    var = pd.DataFrame(index=["g0", "g1"])
    X = csr_matrix(np.ones((n, 2), dtype="float32"))
    ad.AnnData(X=X, obs=obs, var=var).write_h5ad(path, compression="gzip")


def test_celltype_writeback_v0_2_9_part1_barcodes_semantic_split(
    tmp_path, monkeypatch,
):
    """With a `rctd_input_barcodes_path` supplied, `in_rctd_input`
    sources from the barcodes file (= cells SUBMITTED to RCTD), and
    `has_rctd_result` sources from the rctd_df (= cells with a result).
    The two counts can differ (they do, in production: 299,585 vs
    255,661 on MH8_2 v0.3).
    """
    import anndata as ad
    import pandas as pd
    from rctd_split.stages.celltype_writeback import run_celltype_writeback
    from rctd_split._internal.layout import (
        spatial_adata_path, run_dir,
    )

    sample_id = "MHTEST"
    run_id = "42"
    output_root = tmp_path / "runs"

    # Query cells: Q0..Q5 (6 cells).
    ref_p = spatial_adata_path(output_root, sample_id, run_id, "proseg_purified")
    query_p = spatial_adata_path(output_root, sample_id, run_id, "xenium_ranger")
    ref_p.parent.mkdir(parents=True, exist_ok=True)
    # Purified ref has Q0 + Q1 only (= 2 cells "in purified").
    _write_reference_with_purification_status(
        ref_p, sample_id,
        sel_ids=("Q0", "Q1"),
        pstatus=("singlet", "doublet_certain"),
    )
    _write_query(query_p, [
        [0.0, 0.0], [1.0, 0.0], [2.0, 0.0],
        [3.0, 0.0], [4.0, 0.0], [5.0, 0.0],
    ])

    # Barcodes file: Q0..Q3 were submitted to RCTD (= 4 cells).
    split_prep = run_dir(output_root, sample_id, run_id) / "split_prep"
    split_prep.mkdir(parents=True)
    barc_p = split_prep / f"{sample_id}_test_barcodes.tsv.gz"
    _write_barcodes_tsv_gz(barc_p, [f"Q{i}" for i in range(4)])

    # RCTD results: Q0..Q2 got a result (= 3 cells). Q3 was submitted
    # but spacexr dropped it (mimicking the 43,924-drop behavior on the
    # real run).
    synth = pd.DataFrame(
        {
            "spot_class":   ["singlet", "doublet_certain", "reject"],
            "first_type":   ["tumor",  "stroma",           "immune"],
            "min_score":    [0.9,      0.8,                0.3],
        },
        index=pd.Index([f"Q{i}" for i in range(3)], name="cell_id"),
    )

    def _fake_extract(*, rctd_results_rds, rscript_bin, r_lib_paths,
                      tmp_csv=None):
        return synth

    import rctd_split._internal.rctd_results as _rr
    monkeypatch.setattr(_rr, "extract_rctd_results_df", _fake_extract)

    run_celltype_writeback(
        sample_id=sample_id, run_id=run_id, output_root=output_root,
        xenium_ranger_h5ad=None, reference_label_col="first_type",
        reference_x_col="centroid_x", reference_y_col="centroid_y",
        k=1, algorithm="auto", metric="euclidean", tiebreak="min_dist",
        distance_threshold=None, unmatched_policy="nearest_label",
        celltype_col="celltype", distance_col="celltype_source_distance",
        h5ad_compression="gzip", force_rerun=True,
        rctd_results_rds=tmp_path / "fake.rds",
        rctd_obs_prefix="rctd_",
        rctd_input_barcodes_path=barc_p,
    )

    q = ad.read_h5ad(query_p)

    # The split: in_rctd_input=4, has_rctd_result=3, in_purified=2.
    assert list(q.obs["in_rctd_input"]) == [
        True,  True,  True,  True,  False, False,
    ]
    assert list(q.obs["has_rctd_result"]) == [
        True,  True,  True,  False, False, False,
    ]
    assert list(q.obs["in_purified"]) == [
        True,  True,  False, False, False, False,
    ]

    # Invariant: has_rctd_result ⊆ in_rctd_input ⊇ in_purified.
    import numpy as np
    in_t = q.obs["in_rctd_input"].to_numpy(dtype=bool)
    has_t = q.obs["has_rctd_result"].to_numpy(dtype=bool)
    pur_t = q.obs["in_purified"].to_numpy(dtype=bool)
    assert (has_t <= in_t).all()
    assert (pur_t <= in_t).all()

    # rctd_purification_status sourced from purified h5ad (= Q0/Q1 only).
    pstatus = list(q.obs["rctd_purification_status"])
    assert pstatus[0] == "singlet"
    assert pstatus[1] == "doublet_certain"
    assert pstatus[2] == "" and pstatus[3] == "" and pstatus[4] == "" and pstatus[5] == ""

    # RCTD payload cols land for the 3 cells with results; others get fillers.
    assert list(q.obs["rctd_spot_class"]) == [
        "singlet", "doublet_certain", "reject", "", "", "",
    ]
    assert np.isnan(q.obs["rctd_min_score"].iloc[3:]).all()


def test_celltype_writeback_v0_2_9_part1_pstatus_overrides_rctd_df(
    tmp_path, monkeypatch,
):
    """If the rctd_df happens to include a `purification_status` column
    (older RCTD runs could), the Part-1 code OVERWRITES it with the
    purified-h5ad-sourced value — because the purified h5ad is the
    authoritative source (populated by SPLIT::purify, not just
    SPLIT::run_post_process_RCTD).
    """
    import anndata as ad
    import pandas as pd
    from rctd_split.stages.celltype_writeback import run_celltype_writeback
    from rctd_split._internal.layout import (
        spatial_adata_path, run_dir,
    )

    sample_id = "MHTEST"
    run_id = "42"
    output_root = tmp_path / "runs"

    ref_p = spatial_adata_path(output_root, sample_id, run_id, "proseg_purified")
    query_p = spatial_adata_path(output_root, sample_id, run_id, "xenium_ranger")
    ref_p.parent.mkdir(parents=True, exist_ok=True)
    _write_reference_with_purification_status(
        ref_p, sample_id,
        sel_ids=("Q0",),
        pstatus=("AUTHORITATIVE_FROM_PURIFIED",),
    )
    _write_query(query_p, [[0.0, 0.0]])

    split_prep = run_dir(output_root, sample_id, run_id) / "split_prep"
    split_prep.mkdir(parents=True)
    barc_p = split_prep / f"{sample_id}_test_barcodes.tsv.gz"
    _write_barcodes_tsv_gz(barc_p, ["Q0"])

    # rctd_df claims a bogus pstatus for Q0. Must be overwritten.
    synth = pd.DataFrame(
        {
            "spot_class":          ["singlet"],
            "first_type":          ["tumor"],
            "purification_status": ["BOGUS_FROM_RCTD_DF"],  # this is the one to OVERRIDE
        },
        index=pd.Index(["Q0"], name="cell_id"),
    )

    def _fake_extract(**_kwargs):
        return synth

    import rctd_split._internal.rctd_results as _rr
    monkeypatch.setattr(_rr, "extract_rctd_results_df", _fake_extract)

    run_celltype_writeback(
        sample_id=sample_id, run_id=run_id, output_root=output_root,
        xenium_ranger_h5ad=None, reference_label_col="first_type",
        reference_x_col="centroid_x", reference_y_col="centroid_y",
        k=1, algorithm="auto", metric="euclidean", tiebreak="min_dist",
        distance_threshold=None, unmatched_policy="nearest_label",
        celltype_col="celltype", distance_col="celltype_source_distance",
        h5ad_compression="gzip", force_rerun=True,
        rctd_results_rds=tmp_path / "fake.rds",
        rctd_obs_prefix="rctd_",
        rctd_input_barcodes_path=barc_p,
    )

    q = ad.read_h5ad(query_p)
    assert q.obs["rctd_purification_status"].iloc[0] == "AUTHORITATIVE_FROM_PURIFIED"


def test_celltype_writeback_v0_2_9_part1_fallback_without_barcodes(
    tmp_path, monkeypatch,
):
    """When `rctd_input_barcodes_path=None`, `in_rctd_input` falls
    back to `has_rctd_result` (= the v0.2.9-prior behavior that
    OVERSTATES on xenium_ranger-source runs). Just confirms the
    fallback path is live and doesn't crash.
    """
    import anndata as ad
    import pandas as pd
    from rctd_split.stages.celltype_writeback import run_celltype_writeback

    output_root, sample_id, run_id, _, query_p = _setup_run(tmp_path)

    synth = pd.DataFrame(
        {"spot_class": ["singlet", "doublet_certain", "reject"]},
        index=pd.Index([f"Q{i}" for i in range(3)], name="cell_id"),
    )

    def _fake_extract(**_kwargs):
        return synth

    import rctd_split._internal.rctd_results as _rr
    monkeypatch.setattr(_rr, "extract_rctd_results_df", _fake_extract)

    run_celltype_writeback(
        sample_id=sample_id, run_id=run_id, output_root=output_root,
        xenium_ranger_h5ad=None, reference_label_col="first_type",
        reference_x_col="centroid_x", reference_y_col="centroid_y",
        k=1, algorithm="auto", metric="euclidean", tiebreak="min_dist",
        distance_threshold=None, unmatched_policy="nearest_label",
        celltype_col="celltype", distance_col="celltype_source_distance",
        h5ad_compression="gzip", force_rerun=True,
        rctd_results_rds=tmp_path / "fake.rds",
        # NOTE: no rctd_input_barcodes_path kwarg - triggers the fallback.
    )

    q = ad.read_h5ad(query_p)
    assert list(q.obs["in_rctd_input"]) == list(q.obs["has_rctd_result"])
    assert list(q.obs["in_rctd_input"]) == [True, True, True, False, False, False]


def test_celltype_writeback_v0_2_9_back_compat_without_rctd(tmp_path):
    """When ``rctd_results_rds=None`` (the default — xenium-preprocess
    ran without a rctd_run pass), behavior is identical to v0.2.8 and
    none of the new columns appear.
    """
    import anndata as ad
    from rctd_split.stages.celltype_writeback import run_celltype_writeback

    output_root, sample_id, run_id, _, query_p = _setup_run(tmp_path)

    run_celltype_writeback(
        sample_id=sample_id, run_id=run_id, output_root=output_root,
        xenium_ranger_h5ad=None, reference_label_col="first_type",
        reference_x_col="centroid_x", reference_y_col="centroid_y",
        k=1, algorithm="auto", metric="euclidean", tiebreak="min_dist",
        distance_threshold=None, unmatched_policy="nearest_label",
        celltype_col="celltype", distance_col="celltype_source_distance",
        h5ad_compression="gzip", force_rerun=False,
        rctd_results_rds=None,  # explicit — same as not passing
    )

    q = ad.read_h5ad(query_p)
    added = set(q.obs.columns)
    # Only the v0.2.8 cols land; no rctd_* or in_* appear.
    assert any(c.startswith("celltype") for c in added)
    assert not any(c.startswith("rctd_") for c in added)
    assert "in_rctd_input" not in added
    assert "in_purified" not in added
