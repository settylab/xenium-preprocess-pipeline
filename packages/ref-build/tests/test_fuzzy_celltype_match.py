"""Fuzzy celltype label matching (Tracy request 2026-07-10 comment
`4938815275` on `settylab/TracyY123-nexus#14`).

Two layers of assertion:

  1. **Pure unit tests** on `_internal.celltype_match` — Rule 1
     (exact) / Rule 2 (`unknown_maybe_`) / Rule 3 (composite delimited),
     including the `T/NK` slash-in-name edge case.
  2. **Stage-level end-to-end** covering the two spec asks:
     - census + assemble roll composite/unknown_maybe cells up under
       every expected celltype they name, and assemble dedupes them.
     - Rule 2 fires against `primary_count > 0` when the primary has
       zero exact-`T/NK` cells but 100 `unknown_maybe_T/NK` cells.

The stage tests skip (rather than fail) when anndata / scipy are
missing — same discipline as `test_stage_imports.py`.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest


# --- Pure unit tests on the helper -----------------------------------

def test_matcher_rule1_exact():
    """Rule 1: strict equality always fires."""
    from ref_build._internal.celltype_match import build_matcher

    m = build_matcher("T/NK")
    assert m("T/NK") is True


def test_matcher_rule2_unknown_maybe_prefix():
    """Rule 2: `unknown_maybe_{X}` matches X when `include_unknown_maybe`."""
    from ref_build._internal.celltype_match import build_matcher

    m = build_matcher("T/NK", fuzzy=False, include_unknown_maybe=True)
    assert m("unknown_maybe_T/NK") is True
    # But NOT a random unknown_maybe_ label.
    assert m("unknown_maybe_Fibroblast") is False

    # And when disabled AND fuzzy is off — strict equality only.
    m_off = build_matcher("T/NK", fuzzy=False, include_unknown_maybe=False)
    assert m_off("unknown_maybe_T/NK") is False


def test_matcher_rule3_composite_delimited():
    """Rule 3: token surrounded by `_`, `/`, or string boundary."""
    from ref_build._internal.celltype_match import build_matcher

    m = build_matcher("T/NK", fuzzy=True, include_unknown_maybe=False)
    # `T/NK` embedded between `_` markers
    assert m("B/Plasma_T/NK_rbc") is True
    # `B/Plasma` at start followed by `_`
    m_bp = build_matcher("B/Plasma", fuzzy=True, include_unknown_maybe=False)
    assert m_bp("B/Plasma_T/NK_rbc") is True
    # `rbc` at end preceded by `_`
    m_rbc = build_matcher("rbc", fuzzy=True, include_unknown_maybe=False)
    assert m_rbc("B/Plasma_T/NK_rbc") is True

    # Non-token-boundary substring must NOT match. `Fibroblast_only` doesn't
    # end with `Fibroblast` at a token boundary — no trailing `_` or `/`
    # after `Fibroblast`, but `_only` follows without a delimiter that
    # would isolate `Fibroblast` from the rest.
    m_fib = build_matcher("Fibroblast", fuzzy=True, include_unknown_maybe=False)
    # `Fibroblast_only` — starts with `Fibroblast_`, so `Fibroblast` IS
    # surrounded by `^` and `_`. That IS a token boundary match. Verify.
    assert m_fib("Fibroblast_only") is True
    # But a suffix-only substring should NOT match — `PreFibroblast` has
    # no left-hand delimiter.
    assert m_fib("PreFibroblast") is False


def test_matcher_backcompat_strict_equality():
    """`fuzzy=False, include_unknown_maybe=False` reduces to strict equality
    — the pre-2026-07-10 behavior. Backward-compat guarantee."""
    from ref_build._internal.celltype_match import build_matcher

    m = build_matcher("T/NK", fuzzy=False, include_unknown_maybe=False)
    assert m("T/NK") is True
    assert m("unknown_maybe_T/NK") is False
    assert m("B/Plasma_T/NK_rbc") is False


def test_matching_labels_dedupes_and_orders():
    """`matching_labels` returns first-seen-order, unique matches."""
    from ref_build._internal.celltype_match import matching_labels

    labels = [
        "Fibroblast",
        "unknown_maybe_Fibroblast",
        "T/NK",
        "B/Plasma_T/NK_rbc",
        "Fibroblast",  # duplicate
    ]
    got = matching_labels("T/NK", labels)
    assert got == ["T/NK", "B/Plasma_T/NK_rbc"]


# --- Stage-level end-to-end -----------------------------------------

# Skip stage tests if anndata / pandas / scipy are missing.
ad = pytest.importorskip("anndata")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")


def _synthetic_concat(counts_by_sample_and_label: dict[str, dict[str, int]],
                      celltype_col: str) -> "ad.AnnData":
    """Build a concat AnnData given `{sample_id: {actual_label: n_cells, ...}}`.

    Genes are a tiny synthetic panel with integer counts (assemble's
    integer check passes). Cell ids are unique across (sample, label).
    """
    rows = []
    idx = 0
    for sid, per_label in counts_by_sample_and_label.items():
        for lab, n in per_label.items():
            for _ in range(n):
                rows.append({
                    "cell_id": f"{sid}_c{idx:06d}",
                    "sample_ID": sid,
                    celltype_col: lab,
                })
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


def test_census_fuzzy_counts_composite_and_unknown_maybe(tmp_path):
    """Composite + unknown_maybe cells roll up under `T/NK` in the census.

    Synthetic primary:
        - 1 cell labelled `T/NK`
        - 1 cell labelled `unknown_maybe_T/NK`
        - 1 cell labelled `B/Plasma_T/NK_rbc`
        - 1 cell labelled `Fibroblast_only`

    Expected: primary_count for `T/NK` = 3; primary_count for `B/Plasma`
    = 1; primary_count for `Fibroblast` = 1 (the `_only` suffix is
    token-boundary-delimited); `matched_labels` for `T/NK` lists the
    three actual labels.
    """
    from ref_build.stages.census import run_census

    concat = _synthetic_concat(
        {
            "MHTEST": {
                "T/NK": 1,
                "unknown_maybe_T/NK": 1,
                "B/Plasma_T/NK_rbc": 1,
                "Fibroblast_only": 1,
                # A rule-1 keeper so we don't accidentally trigger the
                # rule-2 no-primary path for every celltype.
                "Myeloid": 200,
            },
        },
        celltype_col="Final_level1_celltype_annotation",
    )
    concat_dir = tmp_path / "MHTEST" / "loaded"
    concat_dir.mkdir(parents=True)
    concat.write_h5ad(concat_dir / "concat.h5ad")

    marker_json = _write_marker_json(
        tmp_path, ["T/NK", "B/Plasma", "Fibroblast", "Myeloid"],
    )

    run_census(
        sample_id="MHTEST",
        concat_h5ad=concat_dir / "concat.h5ad",
        output_root=tmp_path,
        celltype_marker_json=marker_json,
        celltype_col="Final_level1_celltype_annotation",
        donor_borrow_cap=100,
        cell_min_instance=20,
        primary_only_celltypes=["tumor", "liver"],
        force_rerun=True,
        fuzzy_matching=True,
        include_unknown_maybe=True,
    )
    df = pd.read_csv(tmp_path / "MHTEST" / "census" / "census.csv")

    # T/NK: 1 exact + 1 unknown_maybe + 1 composite = 3.
    tnk = df.loc[df["celltype"] == "T/NK"].iloc[0]
    assert int(tnk["primary_count"]) == 3, (
        f"T/NK should catch exact + unknown_maybe + composite = 3, "
        f"got {tnk['primary_count']}\nmatched_labels={tnk['matched_labels']!r}"
    )
    # matched_labels lists all three actual labels.
    matched = set(str(tnk["matched_labels"]).split(","))
    assert matched == {"T/NK", "unknown_maybe_T/NK", "B/Plasma_T/NK_rbc"}, (
        f"matched_labels for T/NK = {matched!r}"
    )

    # B/Plasma: 1 composite cell = 1.
    bp = df.loc[df["celltype"] == "B/Plasma"].iloc[0]
    assert int(bp["primary_count"]) == 1, (
        f"B/Plasma should catch composite `B/Plasma_T/NK_rbc` = 1, "
        f"got {bp['primary_count']}"
    )

    # Fibroblast: `Fibroblast_only` is a token-boundary match — 1.
    fib = df.loc[df["celltype"] == "Fibroblast"].iloc[0]
    assert int(fib["primary_count"]) == 1, (
        f"Fibroblast should token-match `Fibroblast_only` = 1, "
        f"got {fib['primary_count']}"
    )


def test_assemble_dedupes_composite_cells(tmp_path):
    """A composite-labelled cell appears exactly ONCE in the assembled
    reference, even though it's counted under multiple expected celltypes.
    Original label preserved (no relabeling)."""
    from ref_build.stages.assemble import run_assemble
    from ref_build.stages.census import run_census

    concat = _synthetic_concat(
        {
            "MHTEST": {
                # Composite cell that matches T/NK + B/Plasma + rbc.
                "B/Plasma_T/NK_rbc": 1,
                # Something for rule-1 to keep so the run doesn't fail.
                "Myeloid": 200,
            },
        },
        celltype_col="Final_level1_celltype_annotation",
    )
    concat_dir = tmp_path / "MHTEST" / "loaded"
    concat_dir.mkdir(parents=True)
    concat.write_h5ad(concat_dir / "concat.h5ad")

    marker_json = _write_marker_json(
        tmp_path, ["T/NK", "B/Plasma", "rbc", "Myeloid"],
    )

    run_census(
        sample_id="MHTEST",
        concat_h5ad=concat_dir / "concat.h5ad",
        output_root=tmp_path,
        celltype_marker_json=marker_json,
        celltype_col="Final_level1_celltype_annotation",
        donor_borrow_cap=100,
        cell_min_instance=20,
        primary_only_celltypes=["tumor", "liver"],
        force_rerun=True,
        fuzzy_matching=True,
        include_unknown_maybe=True,
    )

    run_assemble(
        sample_id="MHTEST",
        concat_h5ad=concat_dir / "concat.h5ad",
        census_csv=tmp_path / "MHTEST" / "census" / "census.csv",
        output_root=tmp_path,
        celltype_col="Final_level1_celltype_annotation",
        donor_borrow_cap=100,
        random_state=1,
        force_rerun=True,
        fuzzy_matching=True,
        include_unknown_maybe=True,
    )

    ref = ad.read_h5ad(tmp_path / "MHTEST" / "assembled" / "reference.h5ad")
    # Cell ids are unique.
    assert ref.n_obs == len(set(ref.obs_names)), (
        f"assembled reference has duplicate cell ids: n_obs={ref.n_obs}, "
        f"unique={len(set(ref.obs_names))}"
    )
    # The composite cell survives with its ORIGINAL label.
    composite_cells = ref.obs[
        ref.obs["Final_level1_celltype_annotation"] == "B/Plasma_T/NK_rbc"
    ]
    assert len(composite_cells) == 1, (
        f"expected exactly 1 composite-labelled cell in the reference, "
        f"got {len(composite_cells)}"
    )


def test_rule2_fuzzy_fires_with_primary_when_all_unknown_maybe(tmp_path):
    """Primary has 0 exact-`T/NK` cells but 100 `unknown_maybe_T/NK` cells.

    Under fuzzy matching, primary_count for T/NK = 100 → rule 2 does NOT
    fire (no borrow) — the primary is the source, decision is either
    `balanced` (Rule 2 HIGH; primary >= cell_min_instance) or
    `primary_only` (Rule 1; we set donor_borrow_cap=50 so 100 > 50
    forces primary_only).

    The key assertion: fuzzy matching resolves the missing_no_donor
    false-positive Tracy reported in comment 4938815275.
    """
    from ref_build.stages.census import run_census

    concat = _synthetic_concat(
        {
            "MHTEST": {"unknown_maybe_T/NK": 100},
        },
        celltype_col="Final_level1_celltype_annotation",
    )
    concat_dir = tmp_path / "MHTEST" / "loaded"
    concat_dir.mkdir(parents=True)
    concat.write_h5ad(concat_dir / "concat.h5ad")

    marker_json = _write_marker_json(tmp_path, ["T/NK"])

    run_census(
        sample_id="MHTEST",
        concat_h5ad=concat_dir / "concat.h5ad",
        output_root=tmp_path,
        celltype_marker_json=marker_json,
        celltype_col="Final_level1_celltype_annotation",
        donor_borrow_cap=50,
        cell_min_instance=20,
        primary_only_celltypes=["tumor", "liver"],
        force_rerun=True,
        fuzzy_matching=True,
        include_unknown_maybe=True,
    )
    df = pd.read_csv(tmp_path / "MHTEST" / "census" / "census.csv")
    tnk = df.loc[df["celltype"] == "T/NK"].iloc[0]

    # Not `missing_no_donor` — the strict-equality regression Tracy caught.
    assert tnk["decision"] != "missing_no_donor", (
        f"T/NK with 100 unknown_maybe_ cells must NOT be missing_no_donor "
        f"under fuzzy matching (Tracy request 2026-07-10 comment 4938815275). "
        f"Got decision={tnk['decision']!r}"
    )
    # 100 > donor_borrow_cap=50 → primary_only fires (Rule 1).
    assert tnk["decision"] == "primary_only"
    assert int(tnk["primary_count"]) == 100

    # Fuzzy off (strict equality) — the pre-fix regression: 0 primary
    # and no donor → missing_no_donor. Sanity re-check on the same data.
    run_census(
        sample_id="MHTEST",
        concat_h5ad=concat_dir / "concat.h5ad",
        output_root=tmp_path,
        celltype_marker_json=marker_json,
        celltype_col="Final_level1_celltype_annotation",
        donor_borrow_cap=50,
        cell_min_instance=20,
        primary_only_celltypes=["tumor", "liver"],
        force_rerun=True,
        fuzzy_matching=False,
        include_unknown_maybe=False,
    )
    df_strict = pd.read_csv(tmp_path / "MHTEST" / "census" / "census.csv")
    tnk_strict = df_strict.loc[df_strict["celltype"] == "T/NK"].iloc[0]
    assert tnk_strict["decision"] == "missing_no_donor", (
        f"strict-equality (fuzzy=False, include_unknown_maybe=False) must "
        f"reproduce the pre-fix behavior — got {tnk_strict['decision']!r}"
    )
    assert int(tnk_strict["primary_count"]) == 0


def test_roundtrip_realistic_4sample(tmp_path):
    """Round-trip sanity: 4-sample scenario mimicking Tracy's MH9 setup.

    Primary MH9 has zero exact `T/NK` / `B/Plasma` / `rbc` cells, but
    the concat contains composite `B/Plasma_T/NK_rbc` and
    `unknown_maybe_*` labels + donor cells. Assert:

      1. Census picks up T/NK across composite + unknown_maybe.
      2. Assemble dedupes so composite cells appear ONCE.
      3. Final reference AnnData has expected per-celltype counts.
    """
    from ref_build.stages.assemble import run_assemble
    from ref_build.stages.census import run_census

    # Numbers below chosen to exercise all three rules:
    #  - Myeloid > 100 primary → rule 1 (primary_only)
    #  - Fibroblast 30 primary + donor supply → balanced
    #  - T/NK: 0 EXACT primary but composite / unknown_maybe fuzz → 15 primary
    #  - B/Plasma: same composite → 10 primary (shared cells)
    concat = _synthetic_concat(
        {
            "MH9": {
                "Myeloid": 200,
                "Fibroblast": 30,
                "unknown_maybe_T/NK": 5,
                "B/Plasma_T/NK_rbc": 10,
                "liver": 50,
                "endothelial": 20,
            },
            "MH10": {
                "Fibroblast": 40,
                "T/NK": 20,
                "B/Plasma": 15,
                "Myeloid": 30,
            },
            "MH8": {
                "Fibroblast": 10,
                "T/NK": 5,
                "endothelial": 8,
            },
            "MH7": {
                "unknown_maybe_Myeloid": 6,
                "unknown_maybe_endothelial": 4,
            },
        },
        celltype_col="Final_level1_celltype_annotation",
    )
    concat_dir = tmp_path / "MH9" / "loaded"
    concat_dir.mkdir(parents=True)
    concat.write_h5ad(concat_dir / "concat.h5ad")

    marker_json = _write_marker_json(
        tmp_path,
        ["Myeloid", "Fibroblast", "T/NK", "B/Plasma", "endothelial"],
    )

    run_census(
        sample_id="MH9",
        concat_h5ad=concat_dir / "concat.h5ad",
        output_root=tmp_path,
        celltype_marker_json=marker_json,
        celltype_col="Final_level1_celltype_annotation",
        donor_borrow_cap=100,
        cell_min_instance=20,
        primary_only_celltypes=["tumor", "liver"],
        force_rerun=True,
        fuzzy_matching=True,
        include_unknown_maybe=True,
    )
    df = pd.read_csv(tmp_path / "MH9" / "census" / "census.csv")

    # T/NK in MH9: 5 unknown_maybe + 10 composite = 15 primary.
    tnk = df.loc[df["celltype"] == "T/NK"].iloc[0]
    assert int(tnk["primary_count"]) == 15
    # Decision: primary_count=15 < cell_min_instance (20) → Rule 2 LOW.
    assert tnk["decision"] == "hybrid_borrow"

    # B/Plasma in MH9: 10 composite = 10 primary; also donor MH10 with 15.
    bp = df.loc[df["celltype"] == "B/Plasma"].iloc[0]
    assert int(bp["primary_count"]) == 10
    # primary=10 < cell_min_instance (20) → Rule 2 LOW.
    assert bp["decision"] == "hybrid_borrow"

    # Myeloid in MH9: 200 primary → rule 1.
    my = df.loc[df["celltype"] == "Myeloid"].iloc[0]
    assert my["decision"] == "primary_only"
    assert int(my["primary_count"]) == 200

    # endothelial in MH9: 20 primary → Rule 2 HIGH (>= 20).
    en = df.loc[df["celltype"] == "endothelial"].iloc[0]
    assert int(en["primary_count"]) == 20
    assert en["decision"] == "balanced"

    run_assemble(
        sample_id="MH9",
        concat_h5ad=concat_dir / "concat.h5ad",
        census_csv=tmp_path / "MH9" / "census" / "census.csv",
        output_root=tmp_path,
        celltype_col="Final_level1_celltype_annotation",
        donor_borrow_cap=100,
        random_state=1,
        force_rerun=True,
        fuzzy_matching=True,
        include_unknown_maybe=True,
    )
    ref = ad.read_h5ad(tmp_path / "MH9" / "assembled" / "reference.h5ad")

    # Dedupe check: every cell in the reference is unique.
    assert ref.n_obs == len(set(ref.obs_names))

    # The 10 `B/Plasma_T/NK_rbc` primary cells appear ONCE each despite
    # matching T/NK, B/Plasma, and (if it were expected) rbc.
    composite = ref.obs[
        ref.obs["Final_level1_celltype_annotation"] == "B/Plasma_T/NK_rbc"
    ]
    assert len(composite) == 10, (
        f"composite `B/Plasma_T/NK_rbc` cells should dedupe to 10, "
        f"got {len(composite)}"
    )

    # Labels preserved — no relabeling.
    label_set = set(ref.obs["Final_level1_celltype_annotation"].astype(str).unique())
    assert "B/Plasma_T/NK_rbc" in label_set
    assert "unknown_maybe_T/NK" in label_set
