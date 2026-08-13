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
