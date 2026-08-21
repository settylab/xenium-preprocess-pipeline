"""End-to-end unit tests for Tracy's hybrid Rule 2 branch selection
(locked spec on `settylab/TracyY123-nexus#21` comment `5137958844`).

Covers both the census-side decisions and the assemble-side sampler
behavior for the two Rule 2 branches:

  Rule 2 HIGH — primary_count >= cell_min_instance
    * decision: `balanced`
    * assemble: per-donor cap = primary_count (existing sampler)

  Rule 2 LOW — primary_count <  cell_min_instance
    * decision: `hybrid_borrow`
    * assemble: primary + waterfall-balanced donor supplement with
      rebalancing, capped at donor_borrow_cap total

Three LOW-branch sub-cases from the spec:

  Case 1 — both donors have plenty (each >= per_donor_target)
  Case 2 — one donor short, need is rebalanced across the rest
  Case 3 — all donors combined can't reach the cap → take everything
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ad = pytest.importorskip("anndata")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")


def _synthetic_concat(counts_by_sample: dict[str, dict[str, int]],
                      celltype_col: str) -> "ad.AnnData":
    rows = []
    idx = 0
    for sid, per_ct in counts_by_sample.items():
        for ct, n in per_ct.items():
            for _ in range(n):
                rows.append({"cell_id": f"{sid}_c{idx:06d}", "sample_ID": sid, celltype_col: ct})
                idx += 1
    df = pd.DataFrame(rows)
    for col in ("cell_id", "sample_ID", celltype_col):
        df[col] = df[col].astype(object)
    obs = df.set_index("cell_id")
    obs.index = obs.index.astype(object)

    n_cells = len(obs)
    n_genes = 4
    X = np.ones((n_cells, n_genes), dtype=np.int32)
    var = pd.DataFrame(index=pd.Index([f"g{i}" for i in range(n_genes)], dtype=object))

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["counts"] = X.copy()
    return adata


def _write_marker_json(tmp_path: Path, celltypes: list[str]) -> Path:
    p = tmp_path / "markers.json"
    with open(p, "w") as f:
        json.dump({f"{ct}_marker": [f"g{i}" for i in range(3)]
                   for ct in celltypes},
                  f)
    return p


def _write_concat(tmp_path: Path, sample_id: str, concat: "ad.AnnData") -> Path:
    (tmp_path / sample_id / "loaded").mkdir(parents=True, exist_ok=True)
    p = tmp_path / sample_id / "loaded" / "concat.h5ad"
    concat.write_h5ad(p)
    return p


def _read_census(output_root: Path, sample_id: str) -> pd.DataFrame:
    return pd.read_csv(output_root / sample_id / f"{sample_id}_test" / "census" / "census.csv")


def _run_pair(tmp_path, counts_by_sample, expected_celltypes,
              donor_borrow_cap=100, cell_min_instance=20,
              primary_only_celltypes=("tumor", "liver")) -> Path:
    """Helper: run census + assemble on a fresh tmp; return the
    assembled-reference h5ad path."""
    from ref_build.stages.assemble import run_assemble
    from ref_build.stages.census import run_census

    concat = _synthetic_concat(
        counts_by_sample=counts_by_sample,
        celltype_col="Final_level1_celltype_annotation",
    )
    concat_path = _write_concat(tmp_path, "MHTEST", concat)
    marker_json = _write_marker_json(tmp_path, list(expected_celltypes))

    run_census(
        sample_id="MHTEST",
        concat_h5ad=concat_path,
        output_root=tmp_path,
        celltype_marker_json=marker_json,
        celltype_col="Final_level1_celltype_annotation",
        donor_borrow_cap=donor_borrow_cap,
        cell_min_instance=cell_min_instance,
        primary_only_celltypes=list(primary_only_celltypes),
        force_rerun=True,
        run_id="test",
    )
    run_assemble(
        sample_id="MHTEST",
        concat_h5ad=concat_path,
        census_csv=tmp_path / "MHTEST" / "MHTEST_test" / "census" / "census.csv",
        output_root=tmp_path,
        celltype_col="Final_level1_celltype_annotation",
        donor_borrow_cap=donor_borrow_cap,
        random_state=42,
        force_rerun=True,
        run_id="test",
    )
    return tmp_path / "MHTEST" / "MHTEST_test" / "rctd" / "MHTEST_reference_post_rules.h5ad"


# ------------------------------------------------------------------
# Rule 2 HIGH — primary >= cell_min_instance → per-donor cap = primary
# ------------------------------------------------------------------

def test_rule2_high_decision_when_primary_at_or_above_min_instance(tmp_path):
    """primary=50 (>= 20), donors have plenty → decision `balanced`
    (Rule 2 HIGH branch). Per-donor cap = primary_count = 50."""
    from ref_build.stages.census import run_census

    concat = _synthetic_concat(
        counts_by_sample={
            "MHTEST": {"fibroblast": 50},
            "MHDONOR1": {"fibroblast": 200},
            "MHDONOR2": {"fibroblast": 200},
        },
        celltype_col="Final_level1_celltype_annotation",
    )
    concat_path = _write_concat(tmp_path, "MHTEST", concat)
    marker_json = _write_marker_json(tmp_path, ["fibroblast"])

    run_census(
        sample_id="MHTEST",
        concat_h5ad=concat_path,
        output_root=tmp_path,
        celltype_marker_json=marker_json,
        celltype_col="Final_level1_celltype_annotation",
        donor_borrow_cap=100,
        cell_min_instance=20,
        primary_only_celltypes=["tumor", "liver"],
        force_rerun=True,
        run_id="test",
    )
    df = _read_census(tmp_path, "MHTEST")
    fib = df.loc[df["celltype"] == "fibroblast"].iloc[0]
    assert fib["decision"] == "balanced", (
        f"expected `balanced` for primary=50 (>=20); got {fib['decision']!r}"
    )
    assert "rule 2 HIGH" in fib["note"]
    assert "per-donor cap = primary_count = 50" in fib["note"]


def test_rule2_high_sampler_caps_each_donor_at_primary_count(tmp_path):
    """Assemble stage: per-donor cap = primary_count (=60). Each donor
    with 100 cells gets sampled down to 60. Total = 60 + 60 + 60 = 180."""
    ref_path = _run_pair(
        tmp_path,
        counts_by_sample={
            "MHTEST": {"fibroblast": 60},
            "MHDONOR1": {"fibroblast": 100},
            "MHDONOR2": {"fibroblast": 100},
        },
        expected_celltypes=["fibroblast"],
    )
    ref = ad.read_h5ad(ref_path)
    from_primary = int((ref.obs["sample_ID"] == "MHTEST").sum())
    from_d1 = int((ref.obs["sample_ID"] == "MHDONOR1").sum())
    from_d2 = int((ref.obs["sample_ID"] == "MHDONOR2").sum())
    assert from_primary == 60
    assert from_d1 == 60, f"MHDONOR1 capped at primary_count=60; got {from_d1}"
    assert from_d2 == 60, f"MHDONOR2 capped at primary_count=60; got {from_d2}"


# ------------------------------------------------------------------
# Rule 2 LOW — primary < cell_min_instance → hybrid_borrow
# ------------------------------------------------------------------

def test_rule2_low_decision_below_min_instance(tmp_path):
    """primary=10 (< 20), donors have plenty → decision `hybrid_borrow`."""
    from ref_build.stages.census import run_census

    concat = _synthetic_concat(
        counts_by_sample={
            "MHTEST": {"fibroblast": 10},
            "MHDONOR1": {"fibroblast": 200},
            "MHDONOR2": {"fibroblast": 200},
        },
        celltype_col="Final_level1_celltype_annotation",
    )
    concat_path = _write_concat(tmp_path, "MHTEST", concat)
    marker_json = _write_marker_json(tmp_path, ["fibroblast"])

    run_census(
        sample_id="MHTEST",
        concat_h5ad=concat_path,
        output_root=tmp_path,
        celltype_marker_json=marker_json,
        celltype_col="Final_level1_celltype_annotation",
        donor_borrow_cap=100,
        cell_min_instance=20,
        primary_only_celltypes=["tumor", "liver"],
        force_rerun=True,
        run_id="test",
    )
    df = _read_census(tmp_path, "MHTEST")
    fib = df.loc[df["celltype"] == "fibroblast"].iloc[0]
    assert fib["decision"] == "hybrid_borrow"
    assert "rule 2 LOW" in fib["note"]


def test_rule2_low_case1_both_donors_have_plenty(tmp_path):
    """LOW Case 1: primary=10, donors=(200, 200), cap=100.
    need = 100-10 = 90. per_donor_target = 90//2 = 45.
    Both donors have >= 45 → each contributes exactly 45.
    Total: 10 primary + 45 + 45 = 100."""
    ref_path = _run_pair(
        tmp_path,
        counts_by_sample={
            "MHTEST": {"fibroblast": 10},
            "MHDONOR1": {"fibroblast": 200},
            "MHDONOR2": {"fibroblast": 200},
        },
        expected_celltypes=["fibroblast"],
        donor_borrow_cap=100,
        cell_min_instance=20,
    )
    ref = ad.read_h5ad(ref_path)
    from_primary = int((ref.obs["sample_ID"] == "MHTEST").sum())
    from_d1 = int((ref.obs["sample_ID"] == "MHDONOR1").sum())
    from_d2 = int((ref.obs["sample_ID"] == "MHDONOR2").sum())
    assert from_primary == 10
    assert from_d1 == 45, f"MHDONOR1 expected 45, got {from_d1}"
    assert from_d2 == 45, f"MHDONOR2 expected 45, got {from_d2}"
    assert from_primary + from_d1 + from_d2 == 100


def test_rule2_low_case2_one_donor_short_rebalances(tmp_path):
    """LOW Case 2: primary=10, donors=(3, 200, 200), cap=100.
    Iter 1: need=90, n=3, per_donor_target=30. MHDONOR1 has 3 < 30 → take all 3.
    Iter 2: n=2, need=87, per_donor_target=43. Both remaining have >= 43 →
            take 43 each.
    Total donor supplement: 3 + 43 + 43 = 89. Grand total: 10 + 89 = 99."""
    ref_path = _run_pair(
        tmp_path,
        counts_by_sample={
            "MHTEST": {"fibroblast": 10},
            "MHDONOR1": {"fibroblast": 3},
            "MHDONOR2": {"fibroblast": 200},
            "MHDONOR3": {"fibroblast": 200},
        },
        expected_celltypes=["fibroblast"],
        donor_borrow_cap=100,
        cell_min_instance=20,
    )
    ref = ad.read_h5ad(ref_path)
    from_primary = int((ref.obs["sample_ID"] == "MHTEST").sum())
    from_d1 = int((ref.obs["sample_ID"] == "MHDONOR1").sum())
    from_d2 = int((ref.obs["sample_ID"] == "MHDONOR2").sum())
    from_d3 = int((ref.obs["sample_ID"] == "MHDONOR3").sum())
    assert from_primary == 10
    assert from_d1 == 3, f"MHDONOR1 (3 avail) → take all 3; got {from_d1}"
    # After MHDONOR1 taken: need=87, split across 2 donors → 43 each
    # (integer div; last cell dropped because need <= per_donor_target * n).
    assert from_d2 == 43, f"MHDONOR2 rebalanced share = 43; got {from_d2}"
    assert from_d3 == 43, f"MHDONOR3 rebalanced share = 43; got {from_d3}"


def test_rule2_low_case3_all_donors_short_take_everything(tmp_path):
    """LOW Case 3: primary=5, donors=(3, 2), cap=100.
    need=95, n=2, per_donor_target=47. Both donors < 47 → take all 3 + 2.
    Rule 2 LOW total = 5+3+2 = 10, which is still < cell_min_instance,
    but no fallback → decision stays `hybrid_borrow` (under-represented
    warning in note)."""
    ref_path = _run_pair(
        tmp_path,
        # Primary needs another celltype to keep the ref non-empty.
        counts_by_sample={
            "MHTEST": {"fibroblast": 5, "myeloid": 100},
            "MHDONOR1": {"fibroblast": 3},
            "MHDONOR2": {"fibroblast": 2},
        },
        expected_celltypes=["fibroblast", "myeloid"],
        donor_borrow_cap=100,
        cell_min_instance=20,
    )
    ref = ad.read_h5ad(ref_path)
    fib_cells = ref.obs[
        ref.obs["Final_level1_celltype_annotation"] == "fibroblast"
    ]
    assert len(fib_cells) == 10, (
        f"expected 5 primary + 3 + 2 = 10 fibroblast cells; got {len(fib_cells)}"
    )
    from_d1 = int((fib_cells["sample_ID"] == "MHDONOR1").sum())
    from_d2 = int((fib_cells["sample_ID"] == "MHDONOR2").sum())
    assert from_d1 == 3, f"take all MHDONOR1 (3); got {from_d1}"
    assert from_d2 == 2, f"take all MHDONOR2 (2); got {from_d2}"

    # Census-side check: under-represented warning in the note.
    df = _read_census(tmp_path, "MHTEST")
    fib_row = df.loc[df["celltype"] == "fibroblast"].iloc[0]
    assert fib_row["decision"] == "hybrid_borrow"
    assert "still < cell_min_instance" in fib_row["note"], (
        f"expected under-represented warning; got note={fib_row['note']!r}"
    )


def test_rule2_low_primary_zero_still_hybrid_borrow(tmp_path):
    """primary=0 falls into the LOW branch (0 < 20). Donors supply the
    full cap. Decision should still be `hybrid_borrow`, not any legacy
    `borrowed` decision."""
    from ref_build.stages.census import run_census

    concat = _synthetic_concat(
        counts_by_sample={
            "MHTEST": {"tumor": 50, "fibroblast": 0},
            "MHDONOR": {"tumor": 10, "fibroblast": 200},
        },
        celltype_col="Final_level1_celltype_annotation",
    )
    concat_path = _write_concat(tmp_path, "MHTEST", concat)
    marker_json = _write_marker_json(tmp_path, ["fibroblast"])

    run_census(
        sample_id="MHTEST",
        concat_h5ad=concat_path,
        output_root=tmp_path,
        celltype_marker_json=marker_json,
        celltype_col="Final_level1_celltype_annotation",
        donor_borrow_cap=100,
        cell_min_instance=20,
        primary_only_celltypes=["tumor", "liver"],
        force_rerun=True,
        run_id="test",
    )
    df = _read_census(tmp_path, "MHTEST")
    fib = df.loc[df["celltype"] == "fibroblast"].iloc[0]
    assert fib["decision"] == "hybrid_borrow", (
        f"primary=0 with a single donor → hybrid_borrow; got {fib['decision']!r}"
    )
    per_donor = json.loads(fib["per_donor_counts"])
    assert per_donor.get("MHDONOR") == 200


def test_rule2_low_primary_zero_single_donor_takes_full_cap(tmp_path):
    """primary=0, single donor with 200 cells, cap=100. need=100, n=1,
    per_donor_target=100 → sample 100."""
    ref_path = _run_pair(
        tmp_path,
        counts_by_sample={
            "MHTEST": {"tumor": 200, "fibroblast": 0},
            "MHDONOR": {"fibroblast": 500},
        },
        expected_celltypes=["fibroblast", "tumor"],
        donor_borrow_cap=100,
        cell_min_instance=20,
    )
    ref = ad.read_h5ad(ref_path)
    fib_cells = ref.obs[
        ref.obs["Final_level1_celltype_annotation"] == "fibroblast"
    ]
    assert len(fib_cells) == 100, (
        f"expected 100 fibroblast cells (donor supplement to cap); got "
        f"{len(fib_cells)}"
    )
    assert set(fib_cells["sample_ID"]) == {"MHDONOR"}


def test_hybrid_deterministic_under_random_seed(tmp_path):
    """Two identical runs pick the same donor cells (random_seed=42)."""
    def _run_once(tmp: Path) -> list[str]:
        ref_path = _run_pair(
            tmp,
            counts_by_sample={
                "MHTEST": {"fibroblast": 10},
                "MHDONOR1": {"fibroblast": 200},
                "MHDONOR2": {"fibroblast": 200},
            },
            expected_celltypes=["fibroblast"],
        )
        r = ad.read_h5ad(ref_path)
        return sorted(r.obs_names.tolist())

    a = _run_once(tmp_path / "a")
    b = _run_once(tmp_path / "b")
    assert a == b, "hybrid sampler must be deterministic under fixed random_seed"
