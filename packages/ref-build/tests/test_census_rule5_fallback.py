"""End-to-end unit tests for Tracy's Rule 5 under the locked hybrid rule
set (TracyY123-nexus#21 comment 5137958844):

    Rule 5 fires when the Rule 2 LOW total (primary + donor supplement
    from the hybrid sampler) is still below `cell_min_instance` AND a
    fallback donor has the celltype. Assemble tops the total up to
    `donor_borrow_cap` from the fallback donor(s).

Coverage:

  Target-list JSON handling
    - `_marker` / `_markers` suffix stripping (flat schema).
    - nested `{tissue: {celltype: markers}}` + `celltype_target_key`.
    - target-list JSON OVERRIDES marker JSON's keys.
    - default (no target-list) uses marker JSON's keys.

  Rule 5 behavior
    - Rule 2 LOW total >= cell_min_instance → fallback NOT consulted.
    - Rule 2 LOW total < cell_min_instance AND fallback has celltype →
      `fallback_borrow` decision.
    - Rule 2 LOW total < cell_min_instance AND fallback lacks celltype →
      `hybrid_borrow` with an under-represented warning in the note
      (kept the few cells we have; no missing_no_donor unless everything
      is zero).
    - Everything zero (primary + donor + fallback) → `missing_no_donor`.
    - Fallback donors are excluded from Rules 1-4 (per_donor_counts
      does not sum in fallback contribution).

  Assemble
    - `fallback_borrow` tops up to donor_borrow_cap.
    - Round-robin across fallbacks that have the celltype.

  Fallback column auto-detection
    - When the fallback h5ad's celltype column is one of the documented
      alternatives (`refined_celltype`, etc.), the load stage
      resolves it under the primary `celltype_col` name without
      operator intervention.

  CLI
    - `--fallback-donor-h5ad` (multi), `--celltype-target-list`,
      `--celltype-target-key`, `--cell-min-instance` all parse.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ad = pytest.importorskip("anndata")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

def _synthetic_h5ad(tmp_path: Path, name: str, counts_by_ct: dict[str, int],
                    celltype_col: str = "Final_level1_celltype_annotation") -> Path:
    """Write a per-sample h5ad with the given per-celltype cell counts."""
    rows = []
    idx = 0
    for ct, n in counts_by_ct.items():
        for _ in range(n):
            rows.append({"cell_id": f"{name}_c{idx:06d}", celltype_col: ct})
            idx += 1
    df = pd.DataFrame(rows)
    for col in ("cell_id", celltype_col):
        df[col] = df[col].astype(object)
    obs = df.set_index("cell_id")
    obs.index = obs.index.astype(object)

    n_cells = len(obs)
    n_genes = 4
    X = np.ones((n_cells, n_genes), dtype=np.int32)
    var = pd.DataFrame(index=pd.Index([f"g{i}" for i in range(n_genes)], dtype=object))
    a = ad.AnnData(X=X, obs=obs, var=var)
    a.layers["counts"] = X.copy()
    p = tmp_path / f"{name}_preprocessed_scRNA.h5ad"
    a.write_h5ad(p)
    return p


def _write_marker_json(tmp_path: Path, celltypes: list[str],
                       filename: str = "markers.json") -> Path:
    p = tmp_path / filename
    with open(p, "w") as f:
        json.dump({f"{ct}_marker": [f"g{i}" for i in range(3)]
                   for ct in celltypes},
                  f)
    return p


def _write_nested_target_json(tmp_path: Path,
                              tissues: dict[str, list[str]]) -> Path:
    p = tmp_path / "nested_targets.json"
    with open(p, "w") as f:
        json.dump({tissue: {f"{ct}_marker": ["g0", "g1"]
                            for ct in celltypes}
                   for tissue, celltypes in tissues.items()},
                  f)
    return p


def _read_census(output_root: Path, sample_id: str) -> pd.DataFrame:
    return pd.read_csv(output_root / sample_id / f"{sample_id}_test" / "census" / "census.csv")


def _run_load_and_census(
    tmp_path: Path,
    primary_h5ad: Path,
    donor_h5ads: list[Path],
    fallback_donor_h5ads: list[Path],
    marker_json: Path,
    sample_id: str = "MHTEST",
    celltype_target_list: Path | None = None,
    celltype_target_key: str | None = None,
    donor_borrow_cap: int = 100,
    cell_min_instance: int = 20,
    primary_only_celltypes: tuple[str, ...] = ("tumor", "liver"),
) -> Path:
    """Helper: run stage 1 (load) + stage 2 (census). Returns the concat path."""
    from ref_build.stages.census import run_census
    from ref_build.stages.load_primary_and_donors import run_load_primary_and_donors

    run_load_primary_and_donors(
        sample_id=sample_id,
        primary_h5ad=primary_h5ad,
        donor_h5ads=donor_h5ads,
        output_root=tmp_path,
        celltype_col="Final_level1_celltype_annotation",
        force_rerun=True,
        run_id="test",
        tumor_type=None,
        fallback_donor_h5ads=fallback_donor_h5ads,
    )
    concat_h5ad = tmp_path / sample_id / f"{sample_id}_test" / "loaded" / "concat.h5ad"
    run_census(
        sample_id=sample_id,
        concat_h5ad=concat_h5ad,
        output_root=tmp_path,
        celltype_marker_json=marker_json,
        celltype_col="Final_level1_celltype_annotation",
        donor_borrow_cap=donor_borrow_cap,
        cell_min_instance=cell_min_instance,
        primary_only_celltypes=list(primary_only_celltypes),
        force_rerun=True,
        celltype_target_list=celltype_target_list,
        celltype_target_key=celltype_target_key,
        run_id="test",
    )
    return concat_h5ad


# ------------------------------------------------------------------
# Target-list JSON handling
# ------------------------------------------------------------------

def test_target_list_json_keys_extracted_from_flat_schema(tmp_path):
    """`_load_expected_celltypes` strips `_marker`/`_markers` suffixes."""
    from ref_build.stages.census import _load_expected_celltypes

    p = tmp_path / "targets.json"
    with open(p, "w") as f:
        json.dump({
            "B/Plasma_marker": ["CD79A"],
            "T/NK_marker": ["CD3E"],
            "rbc_markers": ["HBB"],       # plural suffix
            "Fibroblast_marker": ["POSTN"],
        }, f)
    out = _load_expected_celltypes(marker_json=p)
    assert out == ["B/Plasma", "T/NK", "rbc", "Fibroblast"], (
        f"expected suffix-stripped names in insertion order; got {out}"
    )


def test_target_list_json_overrides_marker_json_keys(tmp_path):
    """When target-list is set, its keys drive the expected celltype set."""
    from ref_build.stages.census import _load_expected_celltypes

    marker = _write_marker_json(tmp_path, ["fibroblast", "myeloid"], "marker.json")
    target = _write_marker_json(tmp_path, ["T/NK", "B/Plasma", "rbc"], "target.json")

    out = _load_expected_celltypes(marker_json=marker, target_list_json=target)
    assert out == ["T/NK", "B/Plasma", "rbc"], (
        f"target-list must override marker JSON; got {out}"
    )


def test_target_list_json_nested_schema_with_key(tmp_path):
    from ref_build.stages.census import _load_expected_celltypes

    nested = _write_nested_target_json(tmp_path, {
        "Bladder": ["urothelial", "smooth_muscle"],
        "Breast": ["luminal", "basal", "adipocyte"],
    })
    marker = _write_marker_json(tmp_path, ["fibroblast"], "marker.json")

    out = _load_expected_celltypes(
        marker_json=marker,
        target_list_json=nested,
        target_list_key="Breast",
    )
    assert out == ["luminal", "basal", "adipocyte"], (
        f"nested-schema read failed; got {out}"
    )


def test_target_list_json_missing_key_raises(tmp_path):
    from ref_build.stages.census import _load_expected_celltypes

    nested = _write_nested_target_json(tmp_path, {"Bladder": ["a", "b"]})
    marker = _write_marker_json(tmp_path, ["fibroblast"], "marker.json")

    with pytest.raises(SystemExit) as exc:
        _load_expected_celltypes(
            marker_json=marker,
            target_list_json=nested,
            target_list_key="Prostate",
        )
    assert "Prostate" in str(exc.value)


def test_default_still_uses_marker_json_keys(tmp_path):
    from ref_build.stages.census import _load_expected_celltypes

    marker = _write_marker_json(tmp_path, ["A", "B", "C"], "marker.json")
    out = _load_expected_celltypes(marker_json=marker)
    assert out == ["A", "B", "C"]


# ------------------------------------------------------------------
# Rule 5 census-level behavior
# ------------------------------------------------------------------

def test_rule5_no_fallback_reads_when_rule2_total_meets_min_instance(tmp_path):
    """Rule 2 LOW/HIGH total >= cell_min_instance → fallback NOT
    consulted (no fallback_borrow decision), even though a fallback
    donor exists for the celltype."""
    # primary=60 (>= 20) → Rule 2 HIGH; total >> cell_min_instance.
    primary = _synthetic_h5ad(tmp_path, "MHTEST",
                              {"fibroblast": 60, "myeloid": 30})
    donor = _synthetic_h5ad(tmp_path, "MHDONOR",
                            {"fibroblast": 100, "myeloid": 200})
    fb = _synthetic_h5ad(tmp_path, "IMMUNEATLAS",
                         {"fibroblast": 500, "myeloid": 500})
    markers = _write_marker_json(tmp_path, ["fibroblast", "myeloid"])

    _run_load_and_census(
        tmp_path, primary, [donor], [fb], markers,
    )
    df = _read_census(tmp_path, "MHTEST")

    for ct in ("fibroblast", "myeloid"):
        row = df.loc[df["celltype"] == ct].iloc[0]
        assert row["decision"] != "fallback_borrow", (
            f"{ct} should not fire fallback_borrow when Rule 2 total already "
            f">= cell_min_instance; got {row['decision']!r}"
        )


def test_rule5_fires_when_rule2_low_total_below_min_instance(tmp_path):
    """primary=0, donor=0, fallback>0 → Rule 2 LOW total = 0 < 20 →
    `fallback_borrow`."""
    primary = _synthetic_h5ad(tmp_path, "MHTEST", {"myeloid": 200})
    donor = _synthetic_h5ad(tmp_path, "MHDONOR", {"myeloid": 100})
    fb = _synthetic_h5ad(tmp_path, "IMMUNEATLAS",
                         {"T/NK": 500, "B/Plasma": 400})
    targets = _write_marker_json(tmp_path,
                                 ["myeloid", "T/NK", "B/Plasma"],
                                 "targets.json")

    _run_load_and_census(
        tmp_path, primary, [donor], [fb], marker_json=targets,
        celltype_target_list=targets,
    )
    df = _read_census(tmp_path, "MHTEST")

    tnk = df.loc[df["celltype"] == "T/NK"].iloc[0]
    assert tnk["decision"] == "fallback_borrow", (
        f"T/NK missing from primary+donor but in fallback → expected "
        f"fallback_borrow; got {tnk['decision']!r}"
    )
    assert "fallback:IMMUNEATLAS" in tnk["source_samples"]
    per_fb = json.loads(tnk["per_fallback_counts"])
    assert per_fb.get("fallback:IMMUNEATLAS") == 500
    assert "rule 5" in tnk["note"]


def test_rule5_fires_when_primary_plus_donor_below_min_instance(tmp_path):
    """primary=5, donor=5, fallback=500 → Rule 2 LOW total = 10 < 20 →
    `fallback_borrow` (top-up to reach 100)."""
    primary = _synthetic_h5ad(tmp_path, "MHTEST",
                              {"myeloid": 200, "T/NK": 5})
    donor = _synthetic_h5ad(tmp_path, "MHDONOR", {"T/NK": 5})
    fb = _synthetic_h5ad(tmp_path, "IMMUNEATLAS", {"T/NK": 500})
    targets = _write_marker_json(tmp_path, ["myeloid", "T/NK"], "targets.json")

    _run_load_and_census(
        tmp_path, primary, [donor], [fb], marker_json=targets,
        celltype_target_list=targets,
    )
    df = _read_census(tmp_path, "MHTEST")

    tnk = df.loc[df["celltype"] == "T/NK"].iloc[0]
    assert tnk["decision"] == "fallback_borrow", (
        f"primary+donor total=10 < 20 → fallback_borrow; got {tnk['decision']!r}"
    )
    assert "rule 5" in tnk["note"]


def test_rule5_noop_when_fallback_lacks_celltype(tmp_path):
    """primary=5, donor=3, no fallback coverage for T/NK → decision
    is `hybrid_borrow` with an under-represented warning; NOT
    `missing_no_donor` because we still have SOMETHING."""
    primary = _synthetic_h5ad(tmp_path, "MHTEST",
                              {"myeloid": 200, "T/NK": 5})
    donor = _synthetic_h5ad(tmp_path, "MHDONOR", {"T/NK": 3})
    # Fallback has cells but for a different celltype.
    fb = _synthetic_h5ad(tmp_path, "IMMUNEATLAS", {"platelets": 500})
    targets = _write_marker_json(tmp_path, ["myeloid", "T/NK"], "targets.json")

    _run_load_and_census(
        tmp_path, primary, [donor], [fb], marker_json=targets,
        celltype_target_list=targets,
    )
    df = _read_census(tmp_path, "MHTEST")

    tnk = df.loc[df["celltype"] == "T/NK"].iloc[0]
    assert tnk["decision"] == "hybrid_borrow", (
        f"fallback lacks T/NK → keep what we have as hybrid_borrow; "
        f"got {tnk['decision']!r}"
    )
    assert "still < cell_min_instance" in tnk["note"], (
        f"expected under-represented warning; got {tnk['note']!r}"
    )


def test_terminal_missing_no_donor_when_absolutely_everything_empty(
    tmp_path,
):
    """primary+donor+fallback all zero for a target celltype →
    `missing_no_donor` with WARN log."""
    import io
    import re
    import sys

    primary = _synthetic_h5ad(tmp_path, "MHTEST", {"myeloid": 200})
    donor = _synthetic_h5ad(tmp_path, "MHDONOR", {"myeloid": 100})
    fb = _synthetic_h5ad(tmp_path, "IMMUNEATLAS", {"platelets": 500})
    targets = _write_marker_json(tmp_path, ["myeloid", "rbc"], "targets.json")

    buf = io.StringIO()
    old_stdout = sys.stdout
    sys.stdout = buf
    try:
        _run_load_and_census(
            tmp_path, primary, [donor], [fb], marker_json=targets,
            celltype_target_list=targets,
        )
    finally:
        sys.stdout = old_stdout
    logs = buf.getvalue()

    df = _read_census(tmp_path, "MHTEST")
    rbc = df.loc[df["celltype"] == "rbc"].iloc[0]
    assert rbc["decision"] == "missing_no_donor"
    assert re.search(r"WARN.*rbc", logs), (
        f"WARN log must fire; got:\n{logs}"
    )


def test_rule5_fallback_donor_excluded_from_regular_donor_counts(tmp_path):
    """Fallback cells must NOT be counted as regular donor cells for
    Rules 1-4 — otherwise a big immune atlas would trigger Rule 1
    or Rule 2 accidentally."""
    primary = _synthetic_h5ad(tmp_path, "MHTEST", {"T/NK": 30})
    donor = _synthetic_h5ad(tmp_path, "MHDONOR", {"T/NK": 20})
    fb = _synthetic_h5ad(tmp_path, "IMMUNEATLAS", {"T/NK": 5000})
    markers = _write_marker_json(tmp_path, ["T/NK"])

    _run_load_and_census(
        tmp_path, primary, [donor], [fb], markers,
    )
    df = _read_census(tmp_path, "MHTEST")
    row = df.loc[df["celltype"] == "T/NK"].iloc[0]
    per_donor = json.loads(row["per_donor_counts"])
    per_fb = json.loads(row["per_fallback_counts"])
    assert per_donor == {"MHDONOR": 20}, (
        f"regular donor counts must exclude fallback donors; got {per_donor}"
    )
    assert per_fb == {"fallback:IMMUNEATLAS": 5000}
    # primary=30 (>= 20) → Rule 2 HIGH; fallback should NOT fire.
    assert row["decision"] != "fallback_borrow"


def test_load_stage_writes_fallback_ids_sidecar(tmp_path):
    """The load stage always writes `fallback_ids.json` next to the
    concat, containing the fallback sample_IDs (or [])."""
    primary = _synthetic_h5ad(tmp_path, "MHTEST", {"myeloid": 20})
    fb1 = _synthetic_h5ad(tmp_path, "IMMUNEATLAS", {"T/NK": 100})
    fb2 = _synthetic_h5ad(tmp_path, "SECONDFALLBACK", {"B/Plasma": 50})

    from ref_build.stages.load_primary_and_donors import run_load_primary_and_donors
    run_load_primary_and_donors(
        sample_id="MHTEST",
        primary_h5ad=primary,
        donor_h5ads=[],
        output_root=tmp_path,
        celltype_col="Final_level1_celltype_annotation",
        force_rerun=True,
        run_id="test",
        fallback_donor_h5ads=[fb1, fb2],
    )
    sidecar = tmp_path / "MHTEST" / "MHTEST_test" / "loaded" / "fallback_ids.json"
    assert sidecar.exists(), f"sidecar missing at {sidecar}"
    with open(sidecar) as f:
        d = json.load(f)
    assert d == {
        "fallback_sample_ids": ["fallback:IMMUNEATLAS", "fallback:SECONDFALLBACK"],
    }


# ------------------------------------------------------------------
# Assemble-level: fallback_borrow tops up to cap
# ------------------------------------------------------------------

def test_assemble_fallback_borrow_tops_up_to_cap_from_fallback(tmp_path):
    """`fallback_borrow` → primary + Rule 2 LOW donor + fallback top-up
    lands at donor_borrow_cap total."""
    primary = _synthetic_h5ad(tmp_path, "MHTEST", {"myeloid": 200})
    fb1 = _synthetic_h5ad(tmp_path, "IMMUNEATLAS", {"T/NK": 400})
    fb2 = _synthetic_h5ad(tmp_path, "SECONDFALLBACK", {"T/NK": 300})
    targets = _write_marker_json(tmp_path, ["myeloid", "T/NK"], "targets.json")

    concat_path = _run_load_and_census(
        tmp_path, primary, [], [fb1, fb2],
        marker_json=targets,
        celltype_target_list=targets,
    )
    from ref_build.stages.assemble import run_assemble
    run_assemble(
        sample_id="MHTEST",
        concat_h5ad=concat_path,
        census_csv=tmp_path / "MHTEST" / "MHTEST_test" / "census" / "census.csv",
        output_root=tmp_path,
        celltype_col="Final_level1_celltype_annotation",
        donor_borrow_cap=100,
        random_state=42,
        force_rerun=True,
        run_id="test",
    )
    ref_path = tmp_path / "MHTEST" / "MHTEST_test" / "rctd" / "MHTEST_reference_post_rules.h5ad"
    ref = ad.read_h5ad(ref_path)

    tnk = ref.obs[ref.obs["Final_level1_celltype_annotation"] == "T/NK"]
    # No primary or regular donor T/NK — the whole 100 come from
    # fallback top-up.
    assert len(tnk) == 100, (
        f"expected 100 T/NK cells (fallback top-up to cap); got {len(tnk)}"
    )
    tnk_sources = set(tnk["sample_ID"])
    assert tnk_sources <= {"fallback:IMMUNEATLAS", "fallback:SECONDFALLBACK"}, (
        f"T/NK must be fallback-only; got sample_IDs {tnk_sources}"
    )
    # Round-robin across both fallbacks (50/50 split, cap=100).
    assert len(tnk_sources) == 2, (
        f"round-robin should hit both fallbacks; got sources {tnk_sources}"
    )


def test_assemble_fallback_borrow_keeps_primary_and_donor(tmp_path):
    """primary=5, donor=5, fallback=500 → keep the 5+5 from
    primary+donor, then top up 90 from fallback = 100 total."""
    primary = _synthetic_h5ad(tmp_path, "MHTEST",
                              {"myeloid": 200, "T/NK": 5})
    donor = _synthetic_h5ad(tmp_path, "MHDONOR", {"T/NK": 5})
    fb = _synthetic_h5ad(tmp_path, "IMMUNEATLAS", {"T/NK": 500})
    targets = _write_marker_json(tmp_path, ["myeloid", "T/NK"], "targets.json")

    concat_path = _run_load_and_census(
        tmp_path, primary, [donor], [fb], marker_json=targets,
        celltype_target_list=targets,
    )
    from ref_build.stages.assemble import run_assemble
    run_assemble(
        sample_id="MHTEST",
        concat_h5ad=concat_path,
        census_csv=tmp_path / "MHTEST" / "MHTEST_test" / "census" / "census.csv",
        output_root=tmp_path,
        celltype_col="Final_level1_celltype_annotation",
        donor_borrow_cap=100,
        random_state=42,
        force_rerun=True,
        run_id="test",
    )
    ref_path = tmp_path / "MHTEST" / "MHTEST_test" / "rctd" / "MHTEST_reference_post_rules.h5ad"
    ref = ad.read_h5ad(ref_path)

    tnk = ref.obs[ref.obs["Final_level1_celltype_annotation"] == "T/NK"]
    assert len(tnk) == 100, (
        f"expected 5 primary + 5 donor + 90 fallback = 100; got {len(tnk)}"
    )
    assert (tnk["sample_ID"] == "MHTEST").sum() == 5
    assert (tnk["sample_ID"] == "MHDONOR").sum() == 5
    assert (tnk["sample_ID"] == "fallback:IMMUNEATLAS").sum() == 90


# ------------------------------------------------------------------
# Fallback celltype column auto-detection
# ------------------------------------------------------------------

def test_fallback_column_autodetect_refined_celltype(tmp_path):
    """A fallback h5ad whose ONLY celltype column is `refined_celltype`
    is transparently rewritten to `Final_level1_celltype_annotation`
    at load time — no operator intervention needed."""
    from ref_build.stages.load_primary_and_donors import run_load_primary_and_donors

    primary = _synthetic_h5ad(tmp_path, "MHTEST",
                              {"myeloid": 200, "T/NK": 0})
    # Fallback carries `refined_celltype` (Tracy's atlases ship this).
    fb = _synthetic_h5ad(tmp_path, "IMMUNEATLAS", {"T/NK": 500},
                         celltype_col="refined_celltype")

    run_load_primary_and_donors(
        sample_id="MHTEST",
        primary_h5ad=primary,
        donor_h5ads=[],
        output_root=tmp_path,
        celltype_col="Final_level1_celltype_annotation",
        force_rerun=True,
        run_id="test",
        fallback_donor_h5ads=[fb],
    )
    concat_h5ad = tmp_path / "MHTEST" / "MHTEST_test" / "loaded" / "concat.h5ad"
    concat = ad.read_h5ad(concat_h5ad)
    fb_cells = concat.obs[concat.obs["sample_ID"] == "fallback:IMMUNEATLAS"]
    # The fallback's cells should carry the standardised column name.
    assert (fb_cells["Final_level1_celltype_annotation"] == "T/NK").sum() == 500


def test_fallback_column_autodetect_typo_column(tmp_path):
    """The `refined_celltype_update_lymphocyes` typo (real column in
    Tracy's atlases) is included in the auto-detect shortlist."""
    from ref_build.stages.load_primary_and_donors import run_load_primary_and_donors

    primary = _synthetic_h5ad(tmp_path, "MHTEST", {"myeloid": 200})
    fb = _synthetic_h5ad(tmp_path, "TYPOATLAS", {"B/Plasma": 300},
                         celltype_col="refined_celltype_update_lymphocyes")

    run_load_primary_and_donors(
        sample_id="MHTEST",
        primary_h5ad=primary,
        donor_h5ads=[],
        output_root=tmp_path,
        celltype_col="Final_level1_celltype_annotation",
        force_rerun=True,
        run_id="test",
        fallback_donor_h5ads=[fb],
    )
    concat_h5ad = tmp_path / "MHTEST" / "MHTEST_test" / "loaded" / "concat.h5ad"
    concat = ad.read_h5ad(concat_h5ad)
    fb_cells = concat.obs[concat.obs["sample_ID"] == "fallback:TYPOATLAS"]
    assert (fb_cells["Final_level1_celltype_annotation"] == "B/Plasma").sum() == 300


def test_fallback_column_autodetect_fails_loudly_on_unknown(tmp_path):
    """A fallback h5ad with an entirely unrecognised celltype column
    fails LOUDLY (not silently skipped). Guards MH9-audit regression."""
    from ref_build.stages.load_primary_and_donors import run_load_primary_and_donors

    primary = _synthetic_h5ad(tmp_path, "MHTEST", {"myeloid": 200})
    fb = _synthetic_h5ad(tmp_path, "MYSTERYATLAS", {"B/Plasma": 300},
                         celltype_col="totally_custom_column_name")

    with pytest.raises(SystemExit) as exc:
        run_load_primary_and_donors(
            sample_id="MHTEST",
            primary_h5ad=primary,
            donor_h5ads=[],
            output_root=tmp_path,
            celltype_col="Final_level1_celltype_annotation",
            force_rerun=True,
        run_id="test",
            fallback_donor_h5ads=[fb],
        )
    msg = str(exc.value)
    assert "Final_level1_celltype_annotation" in msg, (
        f"error should name the missing standard column; got {msg!r}"
    )
    assert "fallback" in msg.lower(), (
        f"error should identify the file as a fallback; got {msg!r}"
    )


# ------------------------------------------------------------------
# CLI surface
# ------------------------------------------------------------------

def test_cli_parses_fallback_and_target_list_flags():
    """`--fallback-donor-h5ad` (multi), `--celltype-target-list`,
    `--celltype-target-key`, `--cell-min-instance`, `--donor-borrow-cap`
    all parse."""
    from ref_build.cli import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "run",
        "--sample-id", "MH7",
        "--primary-h5ad", "/tmp/primary.h5ad",
        "--celltype-marker-json", "/tmp/markers.json",
        "--output-root", "/tmp/out",
        "--fallback-donor-h5ad", "/tmp/fb1.h5ad",
        "--fallback-donor-h5ad", "/tmp/fb2.h5ad",
        "--donor-borrow-cap", "200",
        "--cell-min-instance", "25",
        "--celltype-target-list", "/tmp/targets.json",
        "--celltype-target-key", "Breast",
    ])
    assert len(args.fallback_donor_h5ads) == 2
    assert str(args.fallback_donor_h5ads[0]).endswith("fb1.h5ad")
    assert args.donor_borrow_cap == 200
    assert args.cell_min_instance == 25
    assert str(args.celltype_target_list).endswith("targets.json")
    assert args.celltype_target_key == "Breast"
