"""Stage 2: per-celltype migration census + audit-trail CSV.

Applies the user's hybrid per-celltype migration rules (locked spec on
(internal issue review):

  - Rule 1 — fires if EITHER the celltype is in
    `primary_only_celltypes` (safety-net override) OR the primary's
    cell count for that celltype exceeds `donor_borrow_cap` (default
    100). Decision: `primary_only` — keep ALL primary cells, no
    donor borrow.
  - Rule 2 (hybrid) — for celltypes NOT caught by Rule 1:
      * `primary_count >= cell_min_instance` (default 20) —
        HIGH branch. Decision `balanced`: keep every primary cell
        plus a per-donor supplement capped at `primary_count`
        cells per donor (the classical Stage-B `donor_balanced_-
        sample_by_reference` recipe).
      * `primary_count <  cell_min_instance` — LOW branch. Decision
        `hybrid_borrow`: keep every primary cell plus a
        total-across-donors supplement capped at `donor_borrow_cap`
        (default 100) using balanced sampling with rebalancing —
        if a donor has fewer cells than its per-donor share, take
        all its cells and redistribute the residual across the
        remaining donors.
  - Rule 3 — fuzzy celltype-label matching. Composite delimited
    labels (`B/Plasma_T/NK_rbc` counts toward `T/NK`) and
    `unknown_maybe_X` counts toward `X`. Delegated to
    `ref_build._internal.celltype_match`. Unchanged.
  - Rule 5 — fires when the Rule 2 LOW total (primary + donor
    supplement) is still below `cell_min_instance` AND a fallback
    donor has the celltype. Decision `fallback_borrow`: also
    borrow from fallback donors to top the total up to
    `donor_borrow_cap`.

Terminal decision when even the fallback lacks the celltype AND
primary + donors have zero cells: `missing_no_donor` — the celltype
is skipped and RCTD won't detect it. A WARN log fires.

Two additional guards:

  - `primary_only_celltypes` (default `[tumor, liver]`) — always
    `primary_only` regardless of the census decision. Guards
    ref-build-summary Caveat §1 (union-not-intersection recombine
    leak).
  - Expected-celltype set is the marker JSON's keys with a trailing
    `_marker` / `_markers` stripped (matches xenium-preprocess's
    `--global-non-tumor-json` shape).

Writes `<output_root>/<sample_id>/census/census.csv` — one row per
celltype, columns:

  celltype, primary_count, per_donor_counts (JSON),
  per_fallback_counts (JSON), decision, source_samples, note,
  matched_labels

The audit trail is human-readable — operators inspect it to answer
"why did we borrow X from Y" without having to re-run the pipeline.
"""
from __future__ import annotations

import json
from pathlib import Path

from ref_build._internal.celltype_match import matching_labels
from ref_build._internal.compat import sentinel_exists
from ref_build._internal.logging import log


def _load_expected_celltypes(
    marker_json: Path,
    target_list_json: Path | None = None,
    target_list_key: str | None = None,
) -> list[str]:
    """Extract the expected-celltype set from a marker (or target-list) JSON.

    Two ways to declare the "required celltype set":

      * `marker_json` (always required) — the marker-gene JSON. Keys are
        `<CellType>_marker` (e.g. `Fibroblast_marker`), values are lists
        of marker genes. This is the default source; celltypes = keys
        with trailing `_marker` / `_markers` stripped.
      * `target_list_json` (optional; from (internal issue review) comment
        (internal issue review)) — a separate JSON whose keys enumerate the CANONICAL
        celltype list the user uses to label the flex + spatial data. When
        set, this OVERRIDES the marker JSON's key set as the target
        list; the marker JSON is still consumed elsewhere (RCTD panel
        selection, downstream markers), but Rules 1-5 iterate over the
        target-list's keys. This lets Rule 5 borrow immune cells (T/NK,
        B/Plasma, rbc) from a fallback donor even when the marker JSON
        doesn't otherwise mention them.

    When `target_list_key` is set, the top-level JSON must be a nested
    dict `{key: {celltype: markers}}` and the celltype list is drawn
    from `d[key]`. Otherwise the JSON is treated as flat.
    """
    src = target_list_json if target_list_json is not None else marker_json
    with open(src) as f:
        d = json.load(f)
    if not isinstance(d, dict):
        raise SystemExit(
            f"[census] celltype-list JSON must be a dict; "
            f"got {type(d).__name__} at {src}"
        )
    if target_list_key is not None:
        if target_list_key not in d:
            raise SystemExit(
                f"[census] --celltype-target-key {target_list_key!r} not present "
                f"in {src} (top-level keys: {list(d.keys())})"
            )
        d = d[target_list_key]
        if not isinstance(d, dict):
            raise SystemExit(
                f"[census] nested value at {target_list_key!r} in {src} "
                f"must be a dict of {{celltype: markers}}; got "
                f"{type(d).__name__}"
            )
    out = []
    for k in d.keys():
        # Strip both `_marker` (singular; most keys) AND `_markers`
        # (plural; some legacy notebook JSONs mix conventions, e.g.
        # `rbc_markers`, `epithelial_markers` in markers_NonTumor_level1.json).
        if k.endswith("_markers"):
            out.append(k[: -len("_markers")])
        elif k.endswith("_marker"):
            out.append(k[: -len("_marker")])
        else:
            # Tolerate keys without the suffix — same shape as some of
            # the older marker JSONs in the notebook family.
            out.append(k)
    return out


def _load_fallback_ids(loaded_dir: Path) -> list[str]:
    """Return the ordered fallback-donor sample_IDs, or [] if none.

    Written by the load stage as a sidecar next to the concat h5ad. The
    census + assemble stages use it to partition concat sample_IDs into
    `regular` vs `fallback` (regular = counted for Rules 1-4; fallback =
    only consulted by Rule 5 when the celltype is missing from
    primary + regular donors).
    """
    path = loaded_dir / "fallback_ids.json"
    if not path.exists():
        return []
    with open(path) as f:
        d = json.load(f)
    return list(d.get("fallback_sample_ids", []))


def _simulate_hybrid_low_total(
    per_donor_counts: dict[str, int],
    need: int,
) -> int:
    """Simulate the Rule 2 LOW allocation, returning the total donor
    cells that WOULD be kept (used by the census to decide whether
    Rule 5 needs to fire).

    Matches `_balanced_borrow_with_rebalance` in assemble.py: iterative
    waterfall — per_donor_target = need // n_donors_remaining; if any
    donor has fewer cells than per_donor_target, take all its cells,
    drop it from the pool, and rebalance the residual need across the
    remaining donors. If per_donor_target is 0 (need < n_donors), pick
    one cell from `need` donors round-robin.
    """
    if need <= 0:
        return 0
    remaining = {d: n for d, n in per_donor_counts.items() if n > 0}
    if not remaining:
        return 0

    kept_total = 0
    residual = need
    while remaining and residual > 0:
        n_donors = len(remaining)
        per_donor_target = residual // n_donors
        if per_donor_target == 0:
            # need < n_donors — pick one cell from `residual` donors.
            kept_total += residual
            residual = 0
            break
        under = [d for d, cnt in remaining.items() if cnt < per_donor_target]
        if under:
            for d in under:
                kept_total += remaining[d]
                residual -= remaining[d]
                del remaining[d]
            continue
        # All donors have >= per_donor_target; take per_donor_target
        # from each and stop.
        take = per_donor_target * n_donors
        kept_total += take
        residual -= take
        break
    return kept_total


def _decide(
    primary_count: int,
    per_donor_counts: dict[str, int],
    per_fallback_counts: dict[str, int],
    donor_borrow_cap: int,
    cell_min_instance: int,
    primary_only: bool,
) -> tuple[str, list[str], str]:
    """Return (decision, source_samples, note) for one celltype under the
    locked hybrid rule set (internal issue review).

    Decisions:
      - primary_only     — Rule 1 fires OR celltype in primary_only_celltypes.
      - balanced         — Rule 2 HIGH: primary >= cell_min_instance;
                           donors that have the celltype supplement per-donor
                           at primary_count cells.
      - hybrid_borrow    — Rule 2 LOW: primary < cell_min_instance; donors
                           supplement to total cap = donor_borrow_cap with
                           balanced rebalance.
      - fallback_borrow  — Rule 2 LOW + Rule 5: primary + donor total is
                           still below cell_min_instance and a fallback
                           donor has the celltype; top up from fallback to
                           reach donor_borrow_cap.
      - missing_no_donor — primary + donor total is 0 AND no fallback has
                           the celltype. WARN + skip.
    """
    # ---- Rule 1 ----
    if primary_only:
        return "primary_only", ["<primary>"], "primary_only_celltypes override"
    if primary_count > donor_borrow_cap:
        return "primary_only", ["<primary>"], (
            f"rule 1 fires (primary_count > donor_borrow_cap = {donor_borrow_cap})"
        )

    donors_with = [d for d, n in per_donor_counts.items() if n > 0]
    fallbacks_with = [f for f, n in per_fallback_counts.items() if n > 0]

    # ---- Rule 2 HIGH ----
    if primary_count >= cell_min_instance:
        sources = ["<primary>"] + donors_with
        note = (
            f"rule 2 HIGH: primary_count ({primary_count}) >= "
            f"cell_min_instance ({cell_min_instance}); per-donor cap = "
            f"primary_count = {primary_count}"
        )
        return "balanced", sources, note

    # ---- Rule 2 LOW ----
    # Simulate the donor allocation to see if Rule 5 should fire.
    need = donor_borrow_cap - primary_count
    donor_kept_total = _simulate_hybrid_low_total(per_donor_counts, need)
    total_after_r2 = primary_count + donor_kept_total

    # ---- Rule 5 top-up ----
    if total_after_r2 < cell_min_instance and fallbacks_with:
        sources = (
            (["<primary>"] if primary_count > 0 else [])
            + donors_with
            + [f"<fallback>{f}" for f in fallbacks_with]
        )
        note = (
            f"rule 5 — rule 2 LOW total ({total_after_r2}) < "
            f"cell_min_instance ({cell_min_instance}); borrowing from "
            f"fallback donor(s) {fallbacks_with} to reach cap = "
            f"{donor_borrow_cap}"
        )
        return "fallback_borrow", sources, note

    if total_after_r2 == 0:
        # primary + donor + fallback all empty for this celltype.
        return "missing_no_donor", [], (
            "no primary, no donor, no fallback — SKIPPED "
            "(RCTD cannot detect this celltype)"
        )

    # Rule 2 LOW stands (fallback either has the ct but total >= 20, or
    # fallback lacks it but donors gave us something).
    sources = (["<primary>"] if primary_count > 0 else []) + donors_with
    note = (
        f"rule 2 LOW: primary_count ({primary_count}) < "
        f"cell_min_instance ({cell_min_instance}); donors supplement "
        f"to total cap = {donor_borrow_cap} (balanced with rebalance); "
        f"projected total = {total_after_r2}"
    )
    if total_after_r2 < cell_min_instance:
        note += (
            f" — WARN: total ({total_after_r2}) still < cell_min_instance "
            f"({cell_min_instance}) but no fallback donor has this celltype"
        )
    return "hybrid_borrow", sources, note


def run_census(
    sample_id: str,
    concat_h5ad: Path,
    output_root: Path,
    celltype_marker_json: Path,
    celltype_col: str,
    donor_borrow_cap: int,
    cell_min_instance: int,
    primary_only_celltypes: list[str],
    force_rerun: bool,
    run_id: str,
    fuzzy_matching: bool = True,
    include_unknown_maybe: bool = True,
    celltype_target_list: Path | None = None,
    celltype_target_key: str | None = None,
) -> Path:
    """Build the per-celltype migration census. Writes
    `<run_dir>/census/census.csv` (intermediate — dropped after
    `rctd_reference_build` succeeds).

    Returns the path to the census CSV.
    """
    import anndata as ad
    import pandas as pd

    from ref_build._internal.layout import run_dir as _run_dir

    out_dir = _run_dir(output_root, sample_id, run_id) / "census"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "census.csv"

    if sentinel_exists(out_path, force_rerun):
        log(f"[census] sentinel exists: {out_path} — skipping "
            f"(pass --force-rerun to re-run).")
        return out_path

    log(f"[census] reading {concat_h5ad}")
    adata = ad.read_h5ad(concat_h5ad)
    log(f"[census]   concat shape: {adata.n_obs} cells × {adata.n_vars} genes")

    if celltype_col not in adata.obs.columns:
        raise SystemExit(
            f"[census] concat is missing the celltype column {celltype_col!r} "
            f"in .obs. Columns present: {list(adata.obs.columns)}"
        )
    if "sample_ID" not in adata.obs.columns:
        raise SystemExit(
            "[census] concat is missing the 'sample_ID' column in .obs — "
            "the load_primary_and_donors stage should have written it. "
            "Re-run stage 1 with --force-rerun."
        )

    expected = _load_expected_celltypes(
        marker_json=celltype_marker_json,
        target_list_json=celltype_target_list,
        target_list_key=celltype_target_key,
    )
    _target_src = (
        celltype_target_list.name if celltype_target_list is not None
        else celltype_marker_json.name
    )
    log(f"[census] expected celltypes ({len(expected)}) from {_target_src}"
        f"{f' [key={celltype_target_key!r}]' if celltype_target_key else ''}: "
        f"{expected}")
    log(f"[census] donor_borrow_cap={donor_borrow_cap} "
        f"cell_min_instance={cell_min_instance}")

    # Extend the iteration set with any `primary_only_celltypes` not
    # already declared by the marker JSON. This lets sample types the
    # JSON doesn't cover (e.g. `liver` for benign-liver primaries,
    # `tumor` when the JSON is non-tumor-only) still land in the
    # reference as primary-only cells, instead of being silently
    # dropped by the assemble stage. (user request 2026-07-09.)
    extra_primary_only = [ct for ct in primary_only_celltypes if ct not in expected]
    if extra_primary_only:
        expected = list(expected) + extra_primary_only
        log(f"[census]   + primary_only_celltypes not in marker JSON: "
            f"{extra_primary_only} (added to iteration as primary-only)")

    sample_ids_present = list(adata.obs["sample_ID"].unique())
    if sample_id not in sample_ids_present:
        raise SystemExit(
            f"[census] the primary sample_id={sample_id!r} is not among the "
            f"sample_ID values present in the concat: {sample_ids_present}. "
            "Check the load_primary_and_donors stage output."
        )
    # Partition sample_IDs: primary | regular donors | fallback donors.
    # Fallback IDs are declared by the load stage in a sidecar next to
    # the concat h5ad (see `_load_fallback_ids`); they are excluded from
    # Rules 1-4 and consulted only by Rule 5.
    fallback_ids = _load_fallback_ids(concat_h5ad.parent)
    donor_ids = [
        s for s in sample_ids_present
        if s != sample_id and s not in fallback_ids
    ]
    fallback_ids_present = [f for f in fallback_ids if f in sample_ids_present]
    log(f"[census]   primary sample_ID:   {sample_id!r}")
    log(f"[census]   donor sample_IDs:    {donor_ids}")
    log(f"[census]   fallback sample_IDs: {fallback_ids_present}")

    log(f"[census] fuzzy_matching={fuzzy_matching} "
        f"include_unknown_maybe={include_unknown_maybe}")
    obs = adata.obs[[celltype_col, "sample_ID"]].copy()
    tally = obs.groupby(["sample_ID", celltype_col], observed=True).size().unstack(fill_value=0)
    all_labels = list(tally.columns) if not tally.empty else []

    rows = []
    for ct in expected:
        matched = matching_labels(
            ct, all_labels,
            fuzzy=fuzzy_matching,
            include_unknown_maybe=include_unknown_maybe,
        )
        primary_count = int(
            sum(tally.loc[sample_id, L] for L in matched
                if sample_id in tally.index and L in tally.columns)
        )
        per_donor = {
            d: int(sum(tally.loc[d, L] for L in matched
                       if d in tally.index and L in tally.columns))
            for d in donor_ids
        }
        per_fallback = {
            f: int(sum(tally.loc[f, L] for L in matched
                       if f in tally.index and L in tally.columns))
            for f in fallback_ids_present
        }
        forced_primary_only = ct in set(primary_only_celltypes)
        decision, sources, note = _decide(
            primary_count=primary_count,
            per_donor_counts=per_donor,
            per_fallback_counts=per_fallback,
            donor_borrow_cap=donor_borrow_cap,
            cell_min_instance=cell_min_instance,
            primary_only=forced_primary_only,
        )

        if decision == "missing_no_donor":
            log(f"[census] WARN celltype {ct!r}: primary + all donors + all "
                f"fallbacks have 0 cells — RCTD won't detect this celltype "
                f"(hybrid rule set edge case, no primary, no donor, "
                f"no fallback)")
        elif decision == "fallback_borrow":
            log(f"[census] rule 5: celltype {ct!r}: rule 2 LOW total below "
                f"cell_min_instance; borrowing from fallback donor(s) "
                f"{[f for f, n in per_fallback.items() if n > 0]} to reach "
                f"cap = {donor_borrow_cap}")

        rows.append({
            "celltype": ct,
            "primary_count": primary_count,
            "per_donor_counts": json.dumps(per_donor, sort_keys=True),
            "per_fallback_counts": json.dumps(per_fallback, sort_keys=True),
            "decision": decision,
            "source_samples": ";".join(sources),
            "note": note,
            "matched_labels": ",".join(matched),
        })

    df = pd.DataFrame(rows, columns=[
        "celltype", "primary_count", "per_donor_counts", "per_fallback_counts",
        "decision", "source_samples", "note", "matched_labels",
    ])

    log(f"[census] decision breakdown:")
    for decision, n in df["decision"].value_counts().sort_index().items():
        log(f"[census]   {decision}: {int(n)}")

    df.to_csv(out_path, index=False)
    log(f"[census] wrote {out_path}")
    return out_path
