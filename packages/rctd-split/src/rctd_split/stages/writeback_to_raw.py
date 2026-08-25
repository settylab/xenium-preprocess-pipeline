"""Sub-stage 7: writeback filter status + unpurified obs onto
xenium-preprocess's raw h5ad.

Three joined operations, all additive-only obs joins on xenium-preprocess's
``<S>_proseg_raw.h5ad``:

**Part 1 — filter-status columns for every raw cell.**
``filter_status.csv`` has one row per unpurified cell; three booleans
(``passed_rctd``, ``passed_split_purify``, ``filtered_by_purification``)
are reindexed onto ``raw.obs_names``. Cells NOT in ``filter_status.csv``
MUST have ``xenium-preprocess's qc_filter`` say they were excluded
(``qc_filtered==False`` in xenium-preprocess's convention: cells that
FAILED the count threshold have ``qc_filtered=False``). Anything else is
an "unaccounted-for cells" fail-loud invariant — a CSV writeup is emitted
and the stage raises.

**Part 2 — fold unpurified.h5ad's .obs onto raw.**
``unpurified.h5ad`` is loaded, its .obs columns are reindexed onto
raw's obs axis (NA where missing) and folded onto ``raw.obs`` under the
configured ``exclude_obs_cols`` policy (default: exclude duplicates and
the three filter-status columns already written by Part 1).

**Part 3 — ``passed_purification`` from proseg_purified.h5ad.**
The postprocess sub-stage (sub-stage 6) runs
``sc.pp.filter_cells(min_counts=...)`` on ``proseg_purified.h5ad`` in
place, dropping cells that fall below the count threshold. Cells still
present after that filter are the "truly-purified" set; recording
``raw.obs_names.isin(purified.obs_names)`` on raw captures the effect
of BOTH ``SPLIT::purify`` AND the postprocess min-counts filter in one
bool column (internal issue review).

Idempotent: writes are overwrites, never appends. Rerunning does not
grow ``raw.obs`` schema or double-write. Sentinel is a marker file
under ``intermediate/adata/``.

**Ordering constraint**: this sub-stage MUST run AFTER xenium-preprocess's
``qc_filter`` sub-stage has written ``.obs['qc_filtered']`` on
``raw.h5ad``, and AFTER sub-stage 6 ``postprocess`` has filtered
``proseg_purified.h5ad`` in place.
"""
from __future__ import annotations

from pathlib import Path

from rctd_split._internal.compat import sentinel_exists
from rctd_split._internal.layout import (
    atomic_write_h5ad,
    intermediate_path,
    spatial_adata_path,
)
from rctd_split._internal.logging import log


def _emit_unaccounted_summary(
    raw_obs,
    unaccounted_mask,
    out_path: Path,
) -> None:
    """Write a CSV summary of raw cells that are neither in
    ``filter_status.csv`` nor xenium-preprocess-QC-dropped. Called only in
    the fail-loud path."""
    import pandas as pd

    diag_cols = [
        c for c in (
            "n_counts", "original_cell_id", "centroid_x", "centroid_y",
            "xenium_cell_id_nn", "xenium_id_match", "qc_filtered",
        )
        if c in raw_obs.columns
    ]
    df = raw_obs.loc[unaccounted_mask, diag_cols].copy() \
        if diag_cols else pd.DataFrame(index=raw_obs.index[unaccounted_mask])
    df["note"] = (
        "not in filter_status.csv AND xenium-preprocess qc_filtered=True — "
        "investigate split_prep.filter_cells vs. xenium-preprocess qc_filter config"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    df.to_csv(tmp, index=True)
    import os as _os
    _os.replace(tmp, out_path)


def _fold_filter_status(
    raw_adata,
    filter_status_csv: Path,
    qc_filtered_col: str,
    unaccounted_out_path: Path,
) -> tuple[int, int]:
    """Part 1: reindex filter_status.csv onto raw.obs_names, fail-loud
    on unaccounted cells, write three boolean cols onto raw.obs.

    Returns ``(n_in_filter_status, n_missing)``.
    """
    import pandas as pd

    log(f"[writeback_to_raw]   reading filter_status {filter_status_csv}")
    fs = pd.read_csv(filter_status_csv, index_col=0)
    # DTYPE TRAP: filter_status.csv's cell_id may be int64 while
    # raw.obs_names is str. Coerce both sides to str for the join.
    fs.index = fs.index.astype(str)
    raw_obs_names_str = raw_adata.obs_names.astype(str)

    joined = fs.reindex(raw_obs_names_str)
    missing_mask = joined["passed_rctd"].isna()
    n_missing = int(missing_mask.sum())
    n_present = int((~missing_mask).sum())

    if n_present != len(fs):
        # Every filter_status.csv row must land on exactly one raw cell.
        raise SystemExit(
            f"[writeback_to_raw] filter_status has {len(fs)} rows "
            f"but only {n_present} align with raw.obs_names after dtype "
            "coercion. Check for a cell-id-scheme mismatch (proseg int vs. "
            "Xenium string)."
        )

    if qc_filtered_col not in raw_adata.obs.columns:
        raise SystemExit(
            f"[writeback_to_raw] raw h5ad has no .obs[{qc_filtered_col!r}]. "
            "Run xenium-preprocess's `qc_filter` sub-stage first."
        )

    # xenium-preprocess convention: qc_filtered=True means "cell PASSED QC and
    # is in the mtx bundle downstream sub-stages consume". Cells with
    # qc_filtered=False were dropped by xenium-preprocess QC. Cells NOT in
    # filter_status.csv should be exactly those with qc_filtered=False.
    qc_passed = raw_adata.obs[qc_filtered_col].to_numpy(dtype=bool)
    unaccounted = missing_mask.to_numpy() & qc_passed
    n_unaccounted = int(unaccounted.sum())

    if n_unaccounted > 0:
        _emit_unaccounted_summary(
            raw_adata.obs, unaccounted, out_path=unaccounted_out_path,
        )
        raise SystemExit(
            f"[writeback_to_raw] {n_unaccounted} raw cells are neither "
            f"in filter_status.csv nor xenium-preprocess-QC-dropped "
            f"(qc_filtered=False). See summary at {unaccounted_out_path}."
        )

    # For cells NOT in filter_status.csv, the three booleans default to
    # False (they were dropped upstream by xenium-preprocess QC).
    for col in ("passed_rctd", "passed_split_purify", "filtered_by_purification"):
        vals = joined[col].to_numpy()
        vals = pd.Series(vals).where(~missing_mask.to_numpy(), other=False)
        raw_adata.obs[col] = vals.astype(bool).to_numpy()

    return n_present, n_missing


def _fold_unpurified_obs(
    raw_adata,
    unpurified_h5ad: Path,
    exclude_obs_cols: list[str],
) -> list[str]:
    """Part 2: reindex ``unpurified.obs`` columns onto ``raw.obs_names``
    with NA for missing cells; overwrite matching columns on
    ``raw.obs``. Returns the list of columns actually copied."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    log(f"[writeback_to_raw]   reading unpurified {unpurified_h5ad}")
    unp = ad.read_h5ad(unpurified_h5ad)
    unp_obs = unp.obs.copy()
    unp_obs.index = unp_obs.index.astype(str)

    exclude = set(exclude_obs_cols)
    candidate_cols = [c for c in unp_obs.columns if c not in exclude]

    raw_obs_names_str = raw_adata.obs_names.astype(str)
    copied: list[str] = []

    for col in candidate_cols:
        series = unp_obs[col]
        reindexed = series.reindex(raw_obs_names_str)

        # NA fill convention (roadmap Part C.2):
        #   * bool → nullable "boolean" (distinguishes NA from False)
        #   * numeric → float NaN
        #   * categorical / object → "" (empty string, per mtx_to_h5ad
        #     convention)
        is_bool = pd.api.types.is_bool_dtype(series.dtype)
        if not is_bool and isinstance(series.dtype, pd.CategoricalDtype):
            is_bool = pd.api.types.is_bool_dtype(series.dtype.categories.dtype)

        if is_bool:
            # Two subtleties here:
            # 1. `reindexed.astype("boolean").to_numpy()` (no dtype arg)
            #    returns an OBJECT numpy array of True / False / pd.NA when
            #    any cell is NA. Assigning that to `raw.obs[col]` makes the
            #    column object-dtype; anndata's writer then dispatches to
            #    `write_vlen_string_array`, which fails with
            #    "Can't implicitly convert non-string objects to strings"
            #    (internal issue review).
            # 2. Semantically, a cell absent from unpurified.obs was dropped
            #    upstream — treating that as `False` matches the Part-1
            #    filter-status convention (cells not in filter_status.csv
            #    default to False on the three bool cols).
            new_col = reindexed.astype("boolean").fillna(False).astype(bool)
        elif pd.api.types.is_numeric_dtype(series.dtype):
            arr = reindexed.to_numpy(dtype=object)
            out = np.empty(len(arr), dtype=float)
            for i, v in enumerate(arr):
                try:
                    out[i] = float(v) if v is not None else np.nan
                except (TypeError, ValueError):
                    out[i] = np.nan
            new_col = pd.Series(out, index=raw_adata.obs_names)
        else:
            new_col = reindexed.astype(object).fillna("").astype(str)
            new_col.index = raw_adata.obs_names

        raw_adata.obs[col] = (
            new_col.to_numpy() if hasattr(new_col, "to_numpy") else np.asarray(new_col)
        )
        copied.append(col)

    # Reword: "excluded" was ambiguous — it was read as "removed from
    # source" (internal issue review). These
    # columns are NEITHER removed from unpurified.obs NOR removed from
    # raw.obs; they simply are not re-copied from unpurified onto raw
    # because raw.obs already carries the full (all-cells) version from
    # xenium-preprocess's proseg / qc_filter / enrich_xenium_id sub-stages.
    # Skipping the copy preserves xenium-preprocess's authoritative values
    # (folding unpurified would overwrite them with a QC-passed subset).
    skipped_present = sorted(c for c in exclude if c in unp_obs.columns)
    log(f"[writeback_to_raw]   folded {len(copied)} unpurified.obs "
        f"columns onto raw.obs; skipped {len(skipped_present)} columns "
        f"already present on raw.obs (from xenium-preprocess) — not "
        f"re-copied to preserve xenium-preprocess values: "
        f"{skipped_present}")
    return copied


def _fold_passed_purification(
    raw_adata,
    purified_h5ad: Path,
) -> int:
    """Part 3: set ``raw.obs['passed_purification'] =
    raw.obs_names.isin(purified.obs_names)``.

    ``proseg_purified.h5ad`` after sub-stage 6 ``postprocess`` contains only
    the cells that survived ``sc.pp.filter_cells(min_counts=...)``, so
    this bool column captures the effect of BOTH ``SPLIT::purify`` AND
    the postprocess min-counts filter. Returns the count of True cells
    (which must equal ``purified.n_obs`` when the filter dropped no
    non-raw cells).
    """
    import anndata as ad

    log(f"[writeback_to_raw]   reading purified {purified_h5ad}")
    purified = ad.read_h5ad(purified_h5ad)
    raw_names_str = raw_adata.obs_names.astype(str)
    purified_names_str = purified.obs_names.astype(str)
    passed = raw_names_str.isin(purified_names_str)
    raw_adata.obs["passed_purification"] = passed.astype(bool)
    return int(passed.sum())


def run_writeback_to_raw(
    sample_id: str,
    run_id: str,
    output_root: Path,
    raw_h5ad: Path | None,
    qc_filtered_col: str,
    exclude_obs_cols: list[str],
    h5ad_compression: str | None,
    force_rerun: bool,
) -> Path:
    """Write filter-status + unpurified.obs back onto xenium-preprocess's
    raw h5ad.

    Idempotent: overwrites obs columns in place. Sentinel path returned.
    """
    import anndata as ad

    if raw_h5ad is None:
        raw_h5ad = spatial_adata_path(
            output_root, sample_id, run_id, "proseg_raw",
        )
    raw_h5ad = Path(raw_h5ad)

    filter_status_csv = intermediate_path(
        output_root, sample_id, run_id, "filter_status_csv",
    )
    unpurified_h5ad = intermediate_path(
        output_root, sample_id, run_id, "unpurified_h5ad",
    )
    purified_h5ad = spatial_adata_path(
        output_root, sample_id, run_id, "proseg_purified",
    )
    sentinel = intermediate_path(
        output_root, sample_id, run_id, "writeback_sentinel",
    )
    unaccounted_out_path = intermediate_path(
        output_root, sample_id, run_id, "writeback_unaccounted",
    )

    if sentinel_exists(sentinel, force_rerun):
        log(f"[writeback_to_raw] sentinel exists: {sentinel} — "
            "skipping (pass --force-rerun to re-run).")
        return sentinel

    if not raw_h5ad.exists():
        raise SystemExit(
            f"[writeback_to_raw] xenium-preprocess raw h5ad not found: "
            f"{raw_h5ad}. Run xenium-preprocess first (or pass "
            f"--xenium-preprocess-raw-h5ad)."
        )
    if not filter_status_csv.exists():
        raise SystemExit(
            f"[writeback_to_raw] filter_status.csv not found: "
            f"{filter_status_csv}. This file lives under intermediate/ "
            f"and (pre-fix) was dropped after celltype_writeback "
            f"completed. Regenerate by re-running the upstream chain: "
            f"`--stages split_purify export_mtx mtx_to_h5ad "
            f"filter_status writeback_to_raw` (add "
            f"`--keep-intermediate` to persist for future re-runs)."
        )
    if not unpurified_h5ad.exists():
        raise SystemExit(
            f"[writeback_to_raw] unpurified h5ad not found: "
            f"{unpurified_h5ad}. This file lives under intermediate/ "
            f"and (pre-fix) was dropped after celltype_writeback "
            f"completed. Regenerate by re-running the upstream chain: "
            f"`--stages split_purify export_mtx mtx_to_h5ad "
            f"writeback_to_raw` (add `--keep-intermediate`)."
        )
    if not purified_h5ad.exists():
        raise SystemExit(
            f"[writeback_to_raw] purified h5ad not found: "
            f"{purified_h5ad}. Run the postprocess stage first."
        )

    log(f"[writeback_to_raw] reading {raw_h5ad}")
    raw = ad.read_h5ad(raw_h5ad)

    n_present, n_missing = _fold_filter_status(
        raw_adata=raw,
        filter_status_csv=filter_status_csv,
        qc_filtered_col=qc_filtered_col,
        unaccounted_out_path=unaccounted_out_path,
    )
    log(f"[writeback_to_raw] filter-status: {n_present} matched, "
        f"{n_missing} raw cells absent (marked False on all three cols)")

    copied = _fold_unpurified_obs(
        raw_adata=raw,
        unpurified_h5ad=unpurified_h5ad,
        exclude_obs_cols=list(exclude_obs_cols),
    )

    n_passed_purification = _fold_passed_purification(
        raw_adata=raw,
        purified_h5ad=purified_h5ad,
    )
    log(f"[writeback_to_raw] passed_purification: "
        f"{n_passed_purification}/{raw.n_obs} raw cells present in "
        f"proseg_purified.h5ad (post-min_counts filter)")

    write_kwargs = {"compression": h5ad_compression} if h5ad_compression else {}
    atomic_write_h5ad(raw, raw_h5ad, **write_kwargs)
    log(f"[writeback_to_raw] wrote augmented {raw_h5ad} "
        f"(+3 filter-status cols +{len(copied)} unpurified obs cols "
        f"+1 passed_purification col)")

    sentinel.parent.mkdir(parents=True, exist_ok=True)
    sentinel.write_text(
        f"writeback_to_raw complete: {sample_id}\n"
        f"n_raw={raw.n_obs} n_filter_status_matched={n_present} "
        f"n_absent_from_filter_status={n_missing} "
        f"n_unpurified_obs_cols_copied={len(copied)} "
        f"n_passed_purification={n_passed_purification}\n"
    )
    log(f"[writeback_to_raw] wrote sentinel {sentinel}")
    return sentinel
