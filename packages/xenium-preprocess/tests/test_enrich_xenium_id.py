"""Synthetic smoke test for the enrich_xenium_id stage.

Round-trip: build a tiny fake proseg h5ad (source of proseg IDs +
centroids) + a tiny fake xenium-ranger h5ad (the target we write onto).
Wire xenium centroids so a known subset lands exactly on proseg cells
and the rest is nudged / far off, then run the stage and check the
three new columns on the xenium adata.

Direction (reversed 2026-08-11, TracyY123-nexus#26 comment 5260213932):
FOR EACH xenium cell → find nearest proseg cell → write proseg_cell_id
onto the xenium h5ad's obs. The proseg h5ad is READ-ONLY.

Deps: anndata + pandas + numpy + scikit-learn. Any missing → skip
(same pattern as test_stage_imports and test_dual_matrix).
"""
from __future__ import annotations

import pytest


def _skip_if_missing(*mods):
    import importlib
    for m in mods:
        try:
            importlib.import_module(m)
        except ImportError:
            pytest.skip(f"missing dependency: {m}")


def _build_fixture(tmp_path):
    """Return (xenium_h5ad_path, proseg_h5ad_path, expected_nn_ids).

    Layout:
      - 3 proseg cells with ids ``proseg_0``, ``proseg_1``, ``proseg_2``
        at coordinates (0,0), (10,10), (20,20).
      - 5 xenium cells at coordinates chosen so each has a known
        nearest proseg cell:
          xen_0  at (0, 0)     → nearest proseg_0 (distance 0)
          xen_1  at (10, 10)   → nearest proseg_1 (distance 0)
          xen_2  at (20, 20)   → nearest proseg_2 (distance 0)
          xen_3  at (0.5, 0.5) → nearest proseg_0 (distance ~0.707)
          xen_4  at (11, 11)   → nearest proseg_1 (distance ~1.414)
    """
    import anndata
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    # --- Proseg h5ad (source of IDs) ---
    n_proseg = 3
    proseg_ids = [f"proseg_{i}" for i in range(n_proseg)]
    proseg_xy = np.array([[0.0, 0.0], [10.0, 10.0], [20.0, 20.0]])

    proseg_obs = pd.DataFrame(index=pd.Index(proseg_ids, name="cell_d"))
    proseg_var = pd.DataFrame(index=[f"gene_{k}" for k in range(3)])
    proseg_X = csr_matrix(np.eye(n_proseg, 3))
    proseg_adata = anndata.AnnData(
        X=proseg_X,
        obs=proseg_obs,
        var=proseg_var,
        obsm={"spatial": proseg_xy},
    )
    proseg_path = tmp_path / "proseg_raw.h5ad"
    proseg_adata.write_h5ad(proseg_path)

    # --- Xenium h5ad (target we write onto) ---
    n_xenium = 5
    xen_ids = [f"xen_{i}" for i in range(n_xenium)]
    xen_xy = np.array(
        [
            [0.0, 0.0],       # exact overlap with proseg_0
            [10.0, 10.0],     # exact overlap with proseg_1
            [20.0, 20.0],     # exact overlap with proseg_2
            [0.5, 0.5],       # nearest proseg_0, dist sqrt(0.5)
            [11.0, 11.0],     # nearest proseg_1, dist sqrt(2)
        ]
    )
    xen_obs = pd.DataFrame(
        {
            "transcript_counts": np.arange(n_xenium, dtype=int) + 1,
            "cell_area": np.arange(n_xenium, dtype=float) + 1.0,
        },
        index=pd.Index(xen_ids, name="cell_id"),
    )
    xen_var = pd.DataFrame(index=[f"gene_{k}" for k in range(4)])
    xen_X = csr_matrix(np.eye(n_xenium, 4))
    xen_adata = anndata.AnnData(
        X=xen_X, obs=xen_obs, var=xen_var, obsm={"spatial": xen_xy},
    )
    xen_path = tmp_path / "xenium_ranger.h5ad"
    xen_adata.write_h5ad(xen_path)

    expected_nn_ids = [
        "proseg_0",
        "proseg_1",
        "proseg_2",
        "proseg_0",
        "proseg_1",
    ]
    return xen_path, proseg_path, expected_nn_ids


def test_enrich_xenium_id_smoke(tmp_path):
    _skip_if_missing("anndata", "sklearn", "scipy")

    import anndata
    from xenium_preprocess.stages.enrich_xenium_id import run_enrich_xenium_id

    xen_path, proseg_path, expected_ids = _build_fixture(tmp_path)

    out_h5ad = run_enrich_xenium_id(
        sample_id="synth",
        run_id="test_run",
        xenium_ranger_h5ad=xen_path,
        proseg_h5ad=proseg_path,
        output_root=tmp_path,
        nn_k=1,
        distance_threshold=None,
        nn_id_col="proseg_cell_id_nn",
        distance_col="proseg_id_nn_distance",
        note_col="proseg_id_nn_note",
        write_back_h5ad=True,
        force_rerun=False,
    )
    # The write target is the xenium h5ad, NOT proseg.
    assert out_h5ad == xen_path

    xen = anndata.read_h5ad(out_h5ad)
    obs = xen.obs

    # New columns land on the xenium h5ad.
    for c in ("proseg_cell_id_nn", "proseg_id_nn_distance",
              "proseg_id_nn_note"):
        assert c in obs.columns, f"missing new column on xenium h5ad: {c}"
    # Pre-existing xenium obs cols survive untouched.
    for c in ("transcript_counts", "cell_area"):
        assert c in obs.columns, f"pre-existing column {c!r} was dropped"

    # NN picked the expected proseg ids.
    nn_ids = obs["proseg_cell_id_nn"].astype(object).tolist()
    assert nn_ids == expected_ids, (nn_ids, expected_ids)

    # Distances: rows 0–2 are exact matches (dist 0); rows 3–4 > 0.
    dist = obs["proseg_id_nn_distance"].to_numpy()
    assert dist.dtype.kind == "f"
    assert (dist >= 0).all()
    for i in range(3):
        assert dist[i] == 0.0, (i, dist[i])
    assert dist[3] > 0.0
    assert dist[4] > 0.0

    # Notes: all `nn_match` when no threshold is set.
    notes = obs["proseg_id_nn_note"].astype(object).tolist()
    for i, n in enumerate(notes):
        assert n == "nn_match", (i, n)

    # The proseg h5ad is READ-ONLY: no new columns and no sentinel there.
    proseg = anndata.read_h5ad(proseg_path)
    for c in ("proseg_cell_id_nn", "proseg_id_nn_distance",
              "proseg_id_nn_note",
              # Old-direction columns must NOT reappear on proseg.
              "xenium_cell_id_nn", "xenium_id_nn_distance",
              "xenium_id_nn_note", "xenium_id_match"):
        assert c not in proseg.obs.columns, (
            f"proseg h5ad unexpectedly carries {c!r} — direction reversed"
        )


def test_enrich_xenium_id_distance_threshold(tmp_path):
    """Cells beyond the distance threshold get null nn id + threshold note."""
    _skip_if_missing("anndata", "sklearn", "scipy")

    import anndata
    from xenium_preprocess.stages.enrich_xenium_id import run_enrich_xenium_id

    xen_path, proseg_path, _ = _build_fixture(tmp_path)

    # Row 4's distance is sqrt(2) ≈ 1.414; threshold at 1.0 catches
    # only rows 4 (and NOT row 3 whose dist ≈ 0.707).
    run_enrich_xenium_id(
        sample_id="synth",
        run_id="test_run",
        xenium_ranger_h5ad=xen_path,
        proseg_h5ad=proseg_path,
        output_root=tmp_path,
        nn_k=1,
        distance_threshold=1.0,
        nn_id_col="proseg_cell_id_nn",
        distance_col="proseg_id_nn_distance",
        note_col="proseg_id_nn_note",
        write_back_h5ad=True,
        force_rerun=False,
    )

    xen = anndata.read_h5ad(xen_path)
    obs = xen.obs
    nn_ids = obs["proseg_cell_id_nn"].astype(object).tolist()
    notes = obs["proseg_id_nn_note"].astype(object).tolist()

    # Row 4 got unassigned by threshold.
    assert notes[4] == "unassigned:distance_over_threshold", notes[4]
    assert nn_ids[4] in (None,) or nn_ids[4] != nn_ids[4]  # None or NaN

    # Rows 0..3 stay nn_match.
    for i in range(4):
        assert notes[i] == "nn_match", (i, notes[i])
    # Row 3 within threshold.
    assert nn_ids[3] == "proseg_0", nn_ids[3]


def test_enrich_xenium_id_sentinel_resume(tmp_path):
    """A second call with force_rerun=False is a no-op (sentinel hit)."""
    _skip_if_missing("anndata", "sklearn", "scipy")

    import anndata
    from xenium_preprocess.stages.enrich_xenium_id import run_enrich_xenium_id

    xen_path, proseg_path, _ = _build_fixture(tmp_path)

    run_enrich_xenium_id(
        sample_id="synth",
        run_id="test_run",
        xenium_ranger_h5ad=xen_path,
        proseg_h5ad=proseg_path,
        output_root=tmp_path,
        nn_k=1,
        distance_threshold=None,
        nn_id_col="proseg_cell_id_nn",
        distance_col="proseg_id_nn_distance",
        note_col="proseg_id_nn_note",
        write_back_h5ad=True,
        force_rerun=False,
    )

    # Manually corrupt an output column — sentinel should cause the
    # second call to short-circuit BEFORE reading either h5ad, so the
    # corruption survives (i.e. the stage didn't touch it).
    xen = anndata.read_h5ad(xen_path)
    xen.obs["proseg_cell_id_nn"] = "SENTINEL_CANARY"
    xen.write_h5ad(xen_path)

    run_enrich_xenium_id(
        sample_id="synth",
        run_id="test_run",
        xenium_ranger_h5ad=xen_path,
        proseg_h5ad=proseg_path,
        output_root=tmp_path,
        nn_k=1,
        distance_threshold=None,
        nn_id_col="proseg_cell_id_nn",
        distance_col="proseg_id_nn_distance",
        note_col="proseg_id_nn_note",
        write_back_h5ad=True,
        force_rerun=False,
    )

    xen = anndata.read_h5ad(xen_path)
    assert (xen.obs["proseg_cell_id_nn"] == "SENTINEL_CANARY").all(), (
        "sentinel resume broke — second call re-ran when it should have skipped."
    )


def test_enrich_xenium_id_refuses_nn_k_gt_1(tmp_path):
    _skip_if_missing("anndata", "sklearn", "scipy")

    from xenium_preprocess.stages.enrich_xenium_id import run_enrich_xenium_id

    xen_path, proseg_path, _ = _build_fixture(tmp_path)

    with pytest.raises(SystemExit):
        run_enrich_xenium_id(
            sample_id="synth",
            run_id="test_run",
            xenium_ranger_h5ad=xen_path,
            proseg_h5ad=proseg_path,
            output_root=tmp_path,
            nn_k=3,
            distance_threshold=None,
            nn_id_col="proseg_cell_id_nn",
            distance_col="proseg_id_nn_distance",
            note_col="proseg_id_nn_note",
            write_back_h5ad=True,
            force_rerun=False,
        )


def test_enrich_xenium_id_missing_xenium_h5ad(tmp_path):
    """Fails loud when the xenium-ranger h5ad hasn't been produced yet."""
    _skip_if_missing("anndata", "sklearn", "scipy")

    from xenium_preprocess.stages.enrich_xenium_id import run_enrich_xenium_id

    _, proseg_path, _ = _build_fixture(tmp_path)
    nonexistent_xen = tmp_path / "does_not_exist.h5ad"

    with pytest.raises(SystemExit, match="xenium-ranger h5ad not found"):
        run_enrich_xenium_id(
            sample_id="synth",
            run_id="test_run",
            xenium_ranger_h5ad=nonexistent_xen,
            proseg_h5ad=proseg_path,
            output_root=tmp_path,
            nn_k=1,
            distance_threshold=None,
            nn_id_col="proseg_cell_id_nn",
            distance_col="proseg_id_nn_distance",
            note_col="proseg_id_nn_note",
            write_back_h5ad=True,
            force_rerun=False,
        )
