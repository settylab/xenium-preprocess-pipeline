"""Stage 3: apply the census → merged reference AnnData.

For each celltype in the census, apply the recorded decision
(internal issue review):

  - primary_only     — keep every primary cell of that celltype.
                       Rule 1.
  - balanced         — Rule 2 HIGH. Keep every primary cell; supplement
                       with donor cells via
                       `_donor_balanced_sample_by_reference`
                       (per-donor cap = primary's per-celltype count,
                       `random_state=random_state`).
  - hybrid_borrow    — Rule 2 LOW. Keep every primary cell; supplement
                       with donor cells via
                       `_balanced_borrow_with_rebalance` up to a total
                       cap of `donor_borrow_cap - primary_count` (so the
                       primary + donor total lands at
                       `donor_borrow_cap`). If a donor has fewer cells
                       than its per-donor share, take all its cells and
                       rebalance across remaining donors.
  - fallback_borrow  — Rule 5. Runs the same primary + Rule 2 LOW
                       donor allocation as `hybrid_borrow`, then tops
                       up from fallback donor(s) via the same
                       rebalancing sampler until the total reaches
                       `donor_borrow_cap`.
  - missing_no_donor — skip; the census stage already logged a warning.

Preserves raw integer counts. `.layers["counts"]` (or a fallback layer)
must be integer-valued or the Stage C R side's `require_int=TRUE` will
refuse the counts. This stage refuses to proceed on a non-integer
counts layer rather than silently truncating.

The `noGeneFilter` invariant is preserved — no `sc.pp.filter_genes`
call anywhere on this path.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ref_build._internal.celltype_match import build_matcher
from ref_build._internal.compat import sentinel_exists
from ref_build._internal.layout import atomic_write_h5ad, rctd_path
from ref_build._internal.logging import log


def _donor_balanced_sample_by_reference(
    donor_ids: list[str],
    per_donor_ids: dict[str, list[str]],
    per_donor_cap: int,
    rng: "np.random.Generator",
) -> list[str]:
    """Cap each donor's contribution to `per_donor_cap` cells.

    Faithful port of the `donor_balanced_sample_by_reference` recipe
    described in ref-build-summary-v3.md Stage B step 9: for each donor,
    if they have fewer cells than the cap, keep all; if more, randomly
    sample down to the cap. Returns the concatenated kept obs_names in
    donor order.

    Takes a `Generator` (not a seed) so the caller can hold ONE RNG
    across every celltype's sample. Re-seeding per-celltype would draw
    identical first samples for every celltype (see skeptic finding on
    per-celltype re-seed bias).
    """
    kept: list[str] = []
    for d in donor_ids:
        ids = per_donor_ids.get(d, [])
        if len(ids) <= per_donor_cap:
            kept.extend(ids)
        else:
            idx = np.array(sorted(ids))
            picked = rng.choice(idx, size=per_donor_cap, replace=False)
            kept.extend(picked.tolist())
    return kept


def _balanced_borrow_with_rebalance(
    donor_ids: list[str],
    per_donor_ids: dict[str, list[str]],
    need: int,
    rng: "np.random.Generator",
) -> list[str]:
    """Rule 2 LOW sampler — distribute `need` cells across donors evenly,
    rebalancing when a donor has fewer cells than its share.

    Locked spec (internal issue review):

      * per_donor_target = need // n_donors_remaining
      * If a donor has < per_donor_target → take all its cells, drop it
        from the pool, redistribute the residual need across the rest
      * If per_donor_target == 0 (need < n_donors), pick one cell from
        `need` donors round-robin so no donor gets a systematic
        low-lex advantage across celltypes
      * If all donors combined can't reach `need` → take everything

    Takes a `Generator` (shared across every celltype) for the same
    reasoning as `_donor_balanced_sample_by_reference`.
    """
    if need <= 0 or not donor_ids:
        return []
    # Snapshot available cells per donor; drop empty donors.
    remaining: dict[str, list[str]] = {
        d: sorted(per_donor_ids.get(d, []))
        for d in donor_ids
    }
    remaining = {d: ids for d, ids in remaining.items() if ids}
    if not remaining:
        return []

    kept: list[str] = []
    residual = need
    while remaining and residual > 0:
        n = len(remaining)
        per_donor_target = residual // n
        if per_donor_target == 0:
            # residual < n_donors — pick 1 cell from `residual` donors
            # (deterministic order to match the census's simulated total).
            picked_donors = list(remaining.keys())[:residual]
            for d in picked_donors:
                idx = np.array(remaining[d])
                pick = rng.choice(idx, size=1, replace=False)
                kept.extend(pick.tolist())
            residual = 0
            break
        under = [d for d, ids in remaining.items() if len(ids) < per_donor_target]
        if under:
            for d in under:
                kept.extend(remaining[d])
                residual -= len(remaining[d])
                del remaining[d]
            continue
        # All remaining donors have >= per_donor_target. Sample exactly
        # per_donor_target from each and stop.
        for d, ids in remaining.items():
            if len(ids) == per_donor_target:
                kept.extend(ids)
            else:
                idx = np.array(ids)
                picked = rng.choice(idx, size=per_donor_target, replace=False)
                kept.extend(picked.tolist())
        break
    return kept


def _is_integer(matrix) -> bool:
    """True if `matrix` values are (element-wise) equal to their int cast.

    Handles both sparse (via `scipy.sparse.issparse`) and dense. On
    matrices larger than 1M non-zeros, we check a stratified sample —
    the first, middle, and last thirds — instead of just the head, so a
    stray float in a late gene isn't missed (see skeptic finding on the
    head-only sampler).

    NB: cannot dispatch on `hasattr(matrix, "data")` — dense numpy
    ndarrays also expose `.data` as a memoryview, so that check gave
    a false-positive sparse branch and crashed on `vals.size` for
    dense inputs (fuzzy-match round-trip test flushed it out).
    """
    try:
        from scipy.sparse import issparse as _issparse
    except ImportError:  # pragma: no cover
        _issparse = lambda _m: False  # noqa: E731
    if _issparse(matrix):
        vals = matrix.data
    else:
        vals = np.asarray(matrix).ravel()
    if vals.size == 0:
        return True
    if vals.size <= 1_000_000:
        sample = vals
    else:
        third = 1_000_000 // 3
        mid_start = (vals.size - third) // 2
        sample = np.concatenate([
            vals[:third],
            vals[mid_start:mid_start + third],
            vals[-third:],
        ])
    return bool(np.allclose(sample, np.round(sample), rtol=0, atol=1e-9))


def _pick_counts(adata, layer_hint: str):
    """Return the counts matrix + the name of the layer used.

    Tries `layer_hint` first (default `counts`), then `raw_count`
    (the notebook family's other convention), then `.X`.
    """
    for name in (layer_hint, "raw_count"):
        if name and name in adata.layers:
            return adata.layers[name], name
    return adata.X, "X"


def run_assemble(
    sample_id: str,
    concat_h5ad: Path,
    census_csv: Path,
    output_root: Path,
    celltype_col: str,
    donor_borrow_cap: int,
    random_state: int,
    force_rerun: bool,
    run_id: str,
    tumor_type: str | None = None,
    fuzzy_matching: bool = True,
    include_unknown_maybe: bool = True,
) -> Path:
    """Apply the census. Writes the assembled reference h5ad DIRECTLY to
    the final locked layout path
    `<run_dir>/rctd/<sample_id>_reference_post_rules.h5ad`
    (internal issue review).

    Returns the path to the assembled reference h5ad.
    """
    import anndata as ad
    import pandas as pd

    out_path = rctd_path(output_root, sample_id, run_id, "reference_post_rules")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if sentinel_exists(out_path, force_rerun):
        log(f"[assemble] sentinel exists: {out_path} — skipping "
            f"(pass --force-rerun to re-run).")
        return out_path

    log(f"[assemble] reading {concat_h5ad}")
    adata = ad.read_h5ad(concat_h5ad)
    log(f"[assemble]   concat shape: {adata.n_obs} cells × {adata.n_vars} genes")

    if celltype_col not in adata.obs.columns:
        raise SystemExit(
            f"[assemble] concat is missing the celltype column {celltype_col!r} in .obs."
        )
    if "sample_ID" not in adata.obs.columns:
        raise SystemExit(
            "[assemble] concat is missing the 'sample_ID' column in .obs."
        )

    log(f"[assemble] reading {census_csv}")
    census = pd.read_csv(census_csv)
    log(f"[assemble]   census rows: {len(census)}")

    # ONE RNG for the whole assemble run — reused across every celltype's
    # sample. Per-celltype re-seeding would draw the same first samples
    # for every celltype and bias systematically toward the same
    # lower-lex obs_names.
    rng = np.random.default_rng(random_state)

    obs = adata.obs
    is_primary_mask = (obs["sample_ID"] == sample_id)

    log(f"[assemble] fuzzy_matching={fuzzy_matching} "
        f"include_unknown_maybe={include_unknown_maybe}")
    label_series = obs[celltype_col].astype(object)

    kept_ids: list[str] = []
    per_row_diag: list[dict] = []

    for _, row in census.iterrows():
        ct = row["celltype"]
        decision = row["decision"]
        primary_count = int(row["primary_count"])
        per_donor_counts = json.loads(row["per_donor_counts"])
        try:
            _pfc = row["per_fallback_counts"]
        except (KeyError, IndexError):
            _pfc = None
        if _pfc is None or (isinstance(_pfc, float) and pd.isna(_pfc)):
            per_fallback_counts: dict[str, int] = {}
        else:
            per_fallback_counts = json.loads(_pfc)

        matcher = build_matcher(
            ct, fuzzy=fuzzy_matching, include_unknown_maybe=include_unknown_maybe,
        )
        ct_mask = label_series.map(matcher).astype(bool)
        primary_ct_ids = obs.index[is_primary_mask & ct_mask].tolist()

        donor_ids_with = [d for d, n in per_donor_counts.items() if n > 0]
        per_donor_ids = {
            d: obs.index[(obs["sample_ID"] == d) & ct_mask].tolist()
            for d in donor_ids_with
        }
        fb_ids_with = [f for f, n in per_fallback_counts.items() if n > 0]
        per_fb_ids = {
            f: obs.index[(obs["sample_ID"] == f) & ct_mask].tolist()
            for f in fb_ids_with
        }

        if decision == "primary_only":
            picked = list(primary_ct_ids)
        elif decision == "balanced":
            picked = list(primary_ct_ids)
            supplement = _donor_balanced_sample_by_reference(
                donor_ids=donor_ids_with,
                per_donor_ids=per_donor_ids,
                per_donor_cap=primary_count,
                rng=rng,
            )
            picked.extend(supplement)
        elif decision == "hybrid_borrow":
            picked = list(primary_ct_ids)
            need = donor_borrow_cap - primary_count
            supplement = _balanced_borrow_with_rebalance(
                donor_ids=donor_ids_with,
                per_donor_ids=per_donor_ids,
                need=need,
                rng=rng,
            )
            picked.extend(supplement)
        elif decision == "fallback_borrow":
            # Rule 5: primary + Rule 2 LOW donor allocation, then top up
            # from fallback donors to reach donor_borrow_cap.
            picked = list(primary_ct_ids)
            donor_need = donor_borrow_cap - primary_count
            donor_supplement = _balanced_borrow_with_rebalance(
                donor_ids=donor_ids_with,
                per_donor_ids=per_donor_ids,
                need=donor_need,
                rng=rng,
            )
            picked.extend(donor_supplement)
            fallback_need = donor_borrow_cap - len(picked)
            if fallback_need > 0 and fb_ids_with:
                fallback_supplement = _balanced_borrow_with_rebalance(
                    donor_ids=fb_ids_with,
                    per_donor_ids=per_fb_ids,
                    need=fallback_need,
                    rng=rng,
                )
                picked.extend(fallback_supplement)
        elif decision == "missing_no_donor":
            picked = []
        else:
            raise SystemExit(
                f"[assemble] unknown census decision {decision!r} for celltype {ct!r}. "
                "Was the census.csv written by a compatible version?"
            )

        per_row_diag.append({
            "celltype": ct,
            "decision": decision,
            "kept": len(picked),
        })
        kept_ids.extend(picked)

    log(f"[assemble] decision → kept-cells breakdown:")
    for d in per_row_diag:
        log(f"[assemble]   {d['celltype']:24s} {d['decision']:18s} kept={d['kept']}")

    if not kept_ids:
        raise SystemExit(
            "[assemble] no cells survived the census — check the census.csv and the "
            "primary/donor h5ad files. Every expected celltype resolved to a "
            "missing_no_donor decision or an empty primary."
        )

    # Preserve insertion order without dupes.
    seen: set[str] = set()
    ordered_kept: list[str] = []
    for cid in kept_ids:
        if cid not in seen:
            ordered_kept.append(cid)
            seen.add(cid)

    ref = adata[ordered_kept].copy()
    log(f"[assemble] assembled reference: {ref.n_obs} cells × {ref.n_vars} genes")

    # Integrity check: the counts layer must be integer. The Stage-C R
    # side runs `require_int=TRUE`; casting a float layer to int32 in
    # export_mtx would silently truncate, so we refuse here instead.
    layer_hint = "counts"
    counts, used_layer = _pick_counts(ref, layer_hint)
    if not _is_integer(counts):
        raise SystemExit(
            f"[assemble] the counts layer {used_layer!r} contains non-integer values. "
            "The Stage C R side runs require_int=TRUE and casting to int32 would "
            "silently truncate. Fix the upstream preprocessing to preserve integer "
            "counts, or route a different --export-layer that IS integer."
        )
    log(f"[assemble]   counts source layer: {used_layer!r} (integer-check passed)")

    log(f"[assemble] writing {out_path}")
    atomic_write_h5ad(ref, out_path)
    return out_path
