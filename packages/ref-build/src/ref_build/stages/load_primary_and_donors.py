"""Stage 1: load primary + donor AnnData files → concatenated AnnData.

Adapted from the shared Stage-B backbone in
    the internal reference summary
(Stage B, steps 1-3). Reads each h5ad, concatenates with
`anndata.concat(objs, join="outer", fill_value=0, keys=[sample_ids])`
so the union of genes is preserved (mismatched panels 0-fill), and
attaches a `sample_ID` column indicating which sample each cell came
from.

Rule 5 (internal issue review) — an optional
`fallback_donor_h5ads` list is also concatenated in the same pass,
with a sidecar `fallback_ids.json` next to the concat that lists
those sample_IDs. The census stage reads that sidecar to distinguish
regular donors from fallback donors: fallbacks are excluded from
Rules 1-4 and are only consulted by Rule 5 when a target-list
celltype is missing from primary + regular donors.

The `noGeneFilter` invariant (per the pipeline directive) is preserved:
NO sc.pp.filter_genes(min_cells=20) is called here or anywhere else
in the pipeline. See ref-build-summary-v3.md line 30 and Stage B rctd-split.
"""
from __future__ import annotations

import json
from pathlib import Path

from ref_build._internal.compat import patch_legacy_null_encoding, sentinel_exists
from ref_build._internal.layout import atomic_write_h5ad, run_dir
from ref_build._internal.logging import log


# Fallback donors often come from external immune atlases whose obs
# column names don't match the primary's. Auto-detect: try the user's
# --celltype-col first, then walk this documented shortlist
# (internal issue review). If none match, fail LOUDLY —
# a silent zero-cell fallback would leave every Rule-5 celltype empty
# (see internal audit: the fallback h5ad only had `refined_celltype`, so
# the run silently skipped every Rule-5 borrow).
FALLBACK_CELLTYPE_COL_CANDIDATES = (
    "Final_level1_celltype_annotation",
    "refined_celltype",
    # Yes, "lymphocyes" — the typo is real, kept verbatim to match the
    # actual column name the user's atlases ship with.
    "refined_celltype_update_lymphocyes",
)


def _sample_id_from_path(p: Path) -> str:
    """Best-effort sample id extraction from a `_preprocessed_scRNA.h5ad` filename.

    The notebook family assumes files look like
    `<sample>_preprocessed_scRNA.h5ad`, so we strip the trailing suffix.
    Fall back to the stem if the suffix isn't present.
    """
    name = p.name
    for suffix in ("_preprocessed_scRNA.h5ad", "_scRNA_forSPLIT.h5ad"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    # Fall back to stem (drops one `.h5ad`).
    return p.stem


def _resolve_celltype_col(adata, primary_col: str, fallbacks: list[str], label: str):
    """Ensure ``adata.obs[primary_col]`` exists.

    If ``primary_col`` isn't in ``.obs``, walks ``fallbacks`` in order,
    finds the first candidate present, and COPIES it under ``primary_col``.
    Raises SystemExit if none match.
    """
    if primary_col in adata.obs.columns:
        return
    for fb in fallbacks:
        if fb in adata.obs.columns:
            adata.obs[primary_col] = adata.obs[fb]
            log(f"[load]   {label}: fallback celltype column {fb!r} -> {primary_col!r}")
            return
    raise SystemExit(
        f"[load] {label} is missing the celltype column {primary_col!r} in .obs "
        f"(also tried fallbacks {fallbacks}). Columns present: "
        f"{list(adata.obs.columns)}"
    )


def run_load_primary_and_donors(
    sample_id: str,
    primary_h5ad: Path,
    donor_h5ads: list[Path],
    output_root: Path,
    celltype_col: str,
    force_rerun: bool,
    run_id: str,
    tumor_type: str | None = None,
    celltype_col_fallbacks: list[str] | None = None,
    fallback_donor_h5ads: list[Path] | None = None,
) -> Path:
    """Read + concat primary + donors. Writes the concat h5ad under
    `<run_dir>/loaded/concat.h5ad` (intermediate — dropped after
    `rctd_reference_build` succeeds).

    Returns the path to the concat h5ad.
    """
    import anndata as ad

    patch_legacy_null_encoding()

    out_path = run_dir(output_root, sample_id, run_id) / "loaded" / "concat.h5ad"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if sentinel_exists(out_path, force_rerun):
        log(f"[load] sentinel exists: {out_path} — skipping "
            f"(pass --force-rerun to re-run).")
        return out_path

    log(f"[load] reading primary: {primary_h5ad}")
    primary = ad.read_h5ad(primary_h5ad)
    primary_id = sample_id  # by convention the primary is named after `sample_id`
    log(f"[load]   primary shape: {primary.n_obs} cells × {primary.n_vars} genes; id={primary_id!r}")

    fallbacks = celltype_col_fallbacks or []
    _resolve_celltype_col(primary, celltype_col, fallbacks, label="primary")

    donors = []
    donor_ids = []
    for dp in donor_h5ads:
        did = _sample_id_from_path(dp)
        log(f"[load] reading donor:   {dp}")
        d = ad.read_h5ad(dp)
        log(f"[load]   donor {did!r} shape: {d.n_obs} cells × {d.n_vars} genes")
        _resolve_celltype_col(d, celltype_col, fallbacks, label=f"donor {did!r}")
        donors.append(d)
        donor_ids.append(did)

    # Rule 5 fallback donors — read + concatenated with a distinct
    # sample_ID prefix `fallback:` so downstream stages can tell them
    # apart from regular donors without an obs schema change.
    #
    # Column resolution on fallback h5ad is EXTENDED (locked spec on
    # (internal issue review): try `celltype_col` first,
    # then walk the documented FALLBACK_CELLTYPE_COL_CANDIDATES
    # shortlist. The internal audit surfaced fallback h5ads whose only
    # celltype column was `refined_celltype`; without auto-detect we
    # silently skipped every Rule-5 borrow.
    fallback_paths = list(fallback_donor_h5ads or [])
    fallbacks_ad = []
    fallback_ids: list[str] = []
    fb_fallbacks = list(fallbacks) + [
        c for c in FALLBACK_CELLTYPE_COL_CANDIDATES if c not in fallbacks
    ]
    for fp in fallback_paths:
        did = f"fallback:{_sample_id_from_path(fp)}"
        # Disambiguate if the operator supplies two fallbacks whose
        # `_sample_id_from_path` collide.
        base_did = did
        n = 1
        while did in fallback_ids or did in donor_ids or did == primary_id:
            n += 1
            did = f"{base_did}#{n}"
        log(f"[load] reading fallback: {fp}")
        f = ad.read_h5ad(fp)
        log(f"[load]   fallback {did!r} shape: {f.n_obs} cells × {f.n_vars} genes")
        _resolve_celltype_col(f, celltype_col, fb_fallbacks, label=f"fallback {did!r}")
        fallbacks_ad.append(f)
        fallback_ids.append(did)

    objs = [primary] + donors + fallbacks_ad
    keys = [primary_id] + donor_ids + fallback_ids

    # anndata.concat with join="outer" preserves the union of genes;
    # mismatched panels get 0-filled. The `keys` list drives a MultiIndex
    # on the concatenated obs, which we then collapse into a `sample_ID`
    # obs column (matches Stage B step 2 in the summary).
    log(f"[load] concatenating {len(objs)} objects with join='outer', fill_value=0")
    concat = ad.concat(
        objs,
        join="outer",
        fill_value=0,
        keys=keys,
        label="sample_ID",
        index_unique="-",
    )
    log(f"[load]   concat shape: {concat.n_obs} cells × {concat.n_vars} genes")

    # Emit the per-sample cell tally so the log carries provenance.
    tally = concat.obs["sample_ID"].value_counts().sort_index()
    for sid, n in tally.items():
        log(f"[load]   sample_ID={sid!r}: {int(n)} cells")

    log(f"[load] writing {out_path}")
    atomic_write_h5ad(concat, out_path)

    # Sidecar for Rule 5 — census + assemble read this to partition
    # concat sample_IDs into regular vs fallback. Always written (empty
    # list when no fallbacks) so downstream stages don't have to guard
    # on file existence separately from list emptiness.
    sidecar = out_path.parent / "fallback_ids.json"
    with open(sidecar, "w") as _f:
        json.dump({"fallback_sample_ids": fallback_ids}, _f, indent=2)
    log(f"[load] wrote {sidecar} (fallback sample_IDs: {fallback_ids})")

    return out_path
