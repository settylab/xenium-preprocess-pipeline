"""End-to-end unit tests for Tracy's hybrid per-celltype migration rules
(locked spec on `settylab/TracyY123-nexus#21` comment `5137958844`).

Rule 1 + primary_only_celltypes safety-net coverage lives here. Rule 2
HIGH (per-donor cap), Rule 2 LOW (hybrid with rebalancing), and Rule 5
(fallback top-up when Rule 2 LOW total < cell_min_instance) each have
their own test module (`test_census_hybrid_rule2.py`,
`test_census_rule5_fallback.py`).

Each test constructs a synthetic concat AnnData, writes it + a marker
JSON to a tmp path, invokes `run_census`, and asserts on the resulting
`census.csv`. Skips (rather than fails) when anndata / scipy aren't
installed — same discipline as `test_stage_imports.py`.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

# Skip the whole module if any of the science deps are missing.
ad = pytest.importorskip("anndata")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")


def _synthetic_concat(counts_by_sample: dict[str, dict[str, int]],
                      celltype_col: str) -> "ad.AnnData":
    """Build a concat AnnData with the given per-sample per-celltype cell counts.

    `counts_by_sample` is `{sample_id: {celltype: n_cells, ...}, ...}`.
    """
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


def _read_census(output_root: Path, sample_id: str) -> pd.DataFrame:
    return pd.read_csv(output_root / sample_id / f"{sample_id}_test" / "census" / "census.csv")


def test_rule1_primary_only_when_count_exceeds_cap(tmp_path):
    """Rule 1: primary has 150 tumor + 150 fibroblast, no donors.

    Both counts > donor_borrow_cap (100) → decision is `primary_only`
    for both. And since there ARE no donors, source_samples must be
    `<primary>` only.
    """
    from ref_build.stages.census import run_census

    concat = _synthetic_concat(
        counts_by_sample={"MHTEST": {"tumor": 150, "fibroblast": 150}},
        celltype_col="Final_level1_celltype_annotation",
    )
    (tmp_path / "MHTEST" / "loaded").mkdir(parents=True)
    concat.write_h5ad(tmp_path / "MHTEST" / "loaded" / "concat.h5ad")

    marker_json = _write_marker_json(tmp_path, ["tumor", "fibroblast"])

    run_census(
        sample_id="MHTEST",
        concat_h5ad=tmp_path / "MHTEST" / "loaded" / "concat.h5ad",
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

    # tumor + fibroblast declared in marker JSON; `liver` also gets
    # injected because it's in primary_only_celltypes. All three
    # resolve to `primary_only`.
    decisions = dict(zip(df["celltype"], df["decision"]))
    assert decisions == {
        "tumor": "primary_only",
        "fibroblast": "primary_only",
        "liver": "primary_only",
    }

    fib = df.loc[df["celltype"] == "fibroblast"].iloc[0]
    assert int(fib["primary_count"]) == 150
    assert fib["source_samples"] == "<primary>"
    assert "rule 1 fires" in fib["note"]


def test_rule1_gate_at_exactly_cap_still_reaches_rule2(tmp_path):
    """Rule 1 fires on `primary_count > donor_borrow_cap` (strict).
    At the boundary (primary == cap), we fall through to Rule 2."""
    from ref_build.stages.census import run_census

    concat = _synthetic_concat(
        counts_by_sample={
            "MHTEST": {"fibroblast": 100},
            "MHDONOR": {"fibroblast": 200},
        },
        celltype_col="Final_level1_celltype_annotation",
    )
    (tmp_path / "MHTEST" / "loaded").mkdir(parents=True)
    concat.write_h5ad(tmp_path / "MHTEST" / "loaded" / "concat.h5ad")
    marker_json = _write_marker_json(tmp_path, ["fibroblast"])

    run_census(
        sample_id="MHTEST",
        concat_h5ad=tmp_path / "MHTEST" / "loaded" / "concat.h5ad",
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
    # primary=100 == cap → not primary_only; primary >= 20 → Rule 2 HIGH
    assert fib["decision"] == "balanced", (
        f"primary == cap → Rule 2 HIGH branch; got {fib['decision']!r}"
    )


def test_primary_only_celltypes_override_beats_borrowed(tmp_path):
    """Sanity: even if the primary has 0 tumor cells but a donor does,
    `primary_only_celltypes` override still forces `primary_only`.

    This guards ref-build-summary v3 Caveat §1 — the union-not-
    intersection recombine leak on tumor/liver.
    """
    from ref_build.stages.census import run_census

    concat = _synthetic_concat(
        counts_by_sample={
            "MHTEST": {"tumor": 0, "myeloid": 30},
            "MHDONOR": {"tumor": 500, "myeloid": 40},
        },
        celltype_col="Final_level1_celltype_annotation",
    )
    concat_dir = tmp_path / "MHTEST" / "loaded"
    concat_dir.mkdir(parents=True)
    concat.write_h5ad(concat_dir / "concat.h5ad")

    marker_json = _write_marker_json(tmp_path, ["tumor", "myeloid"])

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
        run_id="test",
    )
    df = _read_census(tmp_path, "MHTEST")
    tumor = df.loc[df["celltype"] == "tumor"].iloc[0]
    assert tumor["decision"] == "primary_only", (
        "primary_only_celltypes override must beat 'borrowed' — "
        "otherwise donor tumor leaks into the reference "
        "(ref-build-summary v3 Caveat §1)."
    )


def test_missing_from_all_including_fallback_stays_missing_no_donor(tmp_path):
    """Terminal case: primary + donor + fallback all zero for a target
    celltype → `missing_no_donor` with WARN log."""
    import io
    import re
    import sys

    from ref_build.stages.census import run_census

    concat = _synthetic_concat(
        counts_by_sample={
            "MHTEST": {"tumor": 10, "myeloid": 20},
            "MHDONOR": {"tumor": 30, "myeloid": 40},
        },
        celltype_col="Final_level1_celltype_annotation",
    )
    concat_dir = tmp_path / "MHTEST" / "loaded"
    concat_dir.mkdir(parents=True)
    concat.write_h5ad(concat_dir / "concat.h5ad")

    marker_json = _write_marker_json(tmp_path, ["fibroblast", "tumor"])

    buf = io.StringIO()
    old_stdout = sys.stdout
    sys.stdout = buf
    try:
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
        run_id="test",
        )
    finally:
        sys.stdout = old_stdout

    logs = buf.getvalue()
    df = _read_census(tmp_path, "MHTEST")

    fib = df.loc[df["celltype"] == "fibroblast"].iloc[0]
    assert fib["decision"] == "missing_no_donor", (
        f"expected missing_no_donor, got {fib['decision']!r}"
    )
    assert int(fib["primary_count"]) == 0
    assert str(fib["source_samples"]) in ("", "nan")
    assert re.search(r"WARN.*fibroblast", logs), (
        f"expected WARN log for the missing_no_donor edge case; got:\n{logs}"
    )
