"""celltype_writeback: NN-map celltype from rctd-split purified onto xenium-preprocess xenium-ranger.

Reads:
  * ``<S>_proseg_purified.h5ad`` (this pipeline's ``mtx_to_h5ad`` output;
    source of celltype via ``.obs[reference_label_col]`` with
    ``.obs[reference_x_col]`` / ``.obs[reference_y_col]`` as reference
    centroids — ``centroid_x`` / ``centroid_y`` by default).
  * ``<S>_xenium_ranger.h5ad`` (xenium-preprocess output; query cells with
    ``.obsm['spatial']`` as ``(n, 2)`` query centroids).
  * ``<S>_rctd_results.rds`` (optional, v0.2.9) — when supplied, per-cell
    RCTD results + SPLIT's post-process columns are written back onto
    xenium_ranger.obs via the ``_internal.rctd_results`` Rscript bridge.

Writes ``.obs[celltype_col]`` and ``.obs[distance_col]`` back on the
xenium-ranger adata via ``atomic_write_h5ad``. Pure NN — no
direct-cell-id-assign, no ``mark_unassigned`` default (per (request on
(internal issue review).

When ``rctd_results_rds`` is non-None (v0.2.9, for the
xenium_ranger-only case where no ``proseg_unpurified.h5ad`` exists to
carry the RCTD results onto, settylab/TracyY123-nexus#40), ALSO writes
these obs columns on xenium_ranger.h5ad.

The provenance schema (v0.2.9 Part 1, c6755625444 — the
investigation-grounded split between "submitted to RCTD" and
"returned a result") tracks every cell through the 7 pipeline
filter stages, counts verified end-to-end on
``MH8_2_BR_2453`` v0.3:

  Stage                                          | N cells
  -----------------------------------------------|---------
  xenium_ranger_to_anndata                       | 322,690
  split_prep.qc_filter (min_counts=10)           | 299,585    <- in_rctd_input=True
  rctd_prep                                      | 299,585
  rctd_run (spacexr create.RCTD + run.RCTD,      |
     drops cells with <10 common-gene counts     |
     after the 5001->4837 gene restriction)      | 255,661    <- has_rctd_result=True
  split_purify.AddMetaData -> unpurified         | 299,585    (NaN meta for the 43,924 w/o rctd_result)
  SPLIT::purify (DO_purify_singlets=TRUE)        | 187,836    (purification_status populated here)
  postprocess qc_filter (min_counts cfg)         | 186,080    <- in_purified=True

Obs columns written on xenium_ranger.h5ad:

  RCTD per-cell (``rctd_`` prefix, 11 cols — sourced from the
  ``@results$results_df`` of the rctd_results.rds AFTER
  ``SPLIT::run_post_process_RCTD``):
    ``rctd_spot_class``, ``rctd_first_type``, ``rctd_second_type``,
    ``rctd_first_class``, ``rctd_second_class``, ``rctd_min_score``,
    ``rctd_singlet_score``, ``rctd_conv_all``, ``rctd_conv_doublet``,
    ``rctd_w1_larger_w2``, ``rctd_same_class``. Columns absent from
    the .rds are silently omitted.

  RCTD purification outcome (``rctd_`` prefix, sourced from the
  purified h5ad's ``obs['purification_status']`` which was
  populated by ``SPLIT::purify``):
    ``rctd_purification_status`` — ``"singlet"`` / etc. for the
    186,080 purified cells; empty string for every other cell.
    LIMITATION: cells with has_rctd_result=True but in_purified=False
    (= 69,581 cells on MH8_2) have LOST ``purification_status``
    because the ``intermediate/split/unpurified.rds`` that would
    carry it is auto-pruned by the terminal qc_report stage. The
    hotfix can't recover this.

  Filter provenance (``in_`` prefix, 3 bools):
    ``in_rctd_input``    : cell was SUBMITTED to RCTD (i.e. in the
                           test_object.rds / split_prep barcodes).
                           Sourced from split_prep barcodes.tsv.gz
                           when ``rctd_input_barcodes_path`` is
                           supplied. Fall-back when not supplied:
                           equals ``has_rctd_result`` (OVERSTATES
                           by the spacexr-dropped cells — see Part 2
                           of settylab/TracyY123-nexus#40 c6755625444).
    ``has_rctd_result``  : cell has a row in rctd_results_df (=
                           what v0.2.9-prior called in_rctd_input).
    ``in_purified``      : cell survived SPLIT::purify + postprocess
                           qc_filter.

Fill behavior on cells NOT in the rctd_input set / NOT in the
results_df:
  - categorical rctd_* cols -> empty string
  - numeric    rctd_* cols -> NaN
  - in_* bools -> False

Idempotent: re-running overwrites the columns in place (never appends
suffixes, never grows adata dims).
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
from rctd_split._internal.nn import nn_map_celltype


def run_celltype_writeback(
    sample_id: str,
    run_id: str,
    output_root: Path,
    xenium_ranger_h5ad: Path | None,
    reference_label_col: str,
    reference_x_col: str,
    reference_y_col: str,
    k: int,
    algorithm: str,
    metric: str,
    tiebreak: str,
    distance_threshold: float | None,
    unmatched_policy: str,
    celltype_col: str,
    distance_col: str,
    h5ad_compression: str | None,
    force_rerun: bool,
    rctd_results_rds: Path | None = None,
    rscript_bin: str = "Rscript",
    r_lib_paths: list[str] | list[Path] | None = None,
    rctd_obs_prefix: str = "rctd_",
    rctd_input_barcodes_path: Path | None = None,
) -> Path:
    """NN-map celltype from purified.h5ad onto xenium_ranger.h5ad.

    Returns the sentinel path.
    """
    import anndata as ad
    import numpy as np

    if xenium_ranger_h5ad is None:
        xenium_ranger_h5ad = spatial_adata_path(
            output_root, sample_id, run_id, "xenium_ranger",
        )
    xenium_ranger_h5ad = Path(xenium_ranger_h5ad)
    purified_h5ad = spatial_adata_path(
        output_root, sample_id, run_id, "proseg_purified",
    )
    sentinel = intermediate_path(
        output_root, sample_id, run_id, "celltype_sentinel",
    )

    if sentinel_exists(sentinel, force_rerun):
        log(f"[celltype_writeback] sentinel exists: {sentinel} — skipping "
            "(pass --force-rerun to re-run).")
        return sentinel

    if not xenium_ranger_h5ad.exists():
        raise SystemExit(
            f"[celltype_writeback] xenium_ranger h5ad not found: "
            f"{xenium_ranger_h5ad}. This stage depends on xenium-preprocess's "
            "xenium_ranger_to_anndata output — run xenium-preprocess first."
        )
    if not purified_h5ad.exists():
        raise SystemExit(
            f"[celltype_writeback] purified h5ad not found: {purified_h5ad}. "
            "Run this pipeline's mtx_to_h5ad + postprocess stages first."
        )

    log(f"[celltype_writeback] reading reference {purified_h5ad}")
    ref = ad.read_h5ad(purified_h5ad)
    log(f"[celltype_writeback] reading query {xenium_ranger_h5ad}")
    query = ad.read_h5ad(xenium_ranger_h5ad)

    for col in (reference_label_col, reference_x_col, reference_y_col):
        if col not in ref.obs.columns:
            raise SystemExit(
                f"[celltype_writeback] reference .obs is missing {col!r}. "
                f"Present columns: {list(ref.obs.columns)}"
            )
    if "spatial" not in query.obsm:
        raise SystemExit(
            "[celltype_writeback] query xenium_ranger.h5ad is missing "
            "obsm['spatial']. Was xenium-preprocess's xenium_ranger_to_anndata run?"
        )

    reference_xy = np.column_stack(
        (
            ref.obs[reference_x_col].to_numpy(dtype=float),
            ref.obs[reference_y_col].to_numpy(dtype=float),
        )
    )
    reference_labels = ref.obs[reference_label_col].to_numpy(dtype=object)
    query_xy = np.asarray(query.obsm["spatial"], dtype=float)

    # Drop reference rows with NaN centroids (would poison NN.fit).
    ref_finite = np.isfinite(reference_xy).all(axis=1)
    if not ref_finite.all():
        n_drop = int((~ref_finite).sum())
        log(f"[celltype_writeback]   reference: dropping {n_drop} rows "
            "with non-finite centroids")
        reference_xy = reference_xy[ref_finite]
        reference_labels = reference_labels[ref_finite]

    # Query cells with NaN spatial: still assign — mark by NaN distance.
    query_finite = np.isfinite(query_xy).all(axis=1)
    if not query_finite.all():
        n_bad = int((~query_finite).sum())
        log(f"[celltype_writeback]   query: {n_bad} cells have non-finite "
            "spatial; they receive NaN distance + '' celltype")

    log(f"[celltype_writeback] fitting NN: k={k}, metric={metric}, "
        f"tiebreak={tiebreak}, unmatched_policy={unmatched_policy}")
    n_query = query.n_obs
    labels = np.array([""] * n_query, dtype=object)
    distances = np.full(n_query, np.nan, dtype=float)
    if reference_xy.shape[0] == 0:
        log("[celltype_writeback] WARN: no reference cells to map from; "
            "all query labels remain empty.")
    elif query_finite.any():
        result = nn_map_celltype(
            reference_xy=reference_xy,
            reference_labels=reference_labels,
            query_xy=query_xy[query_finite],
            k=k,
            algorithm=algorithm,
            metric=metric,
            tiebreak=tiebreak,
            distance_threshold=distance_threshold,
            unmatched_policy=unmatched_policy,
        )
        labels[query_finite] = result["labels"]
        distances[query_finite] = result["distances"]

    # Log a distance summary — the user asked for p50/p90/p99 so a per-
    # sample sanity gate is visible in logs.
    if np.isfinite(distances).any():
        finite = distances[np.isfinite(distances)]
        p50, p90, p99 = np.percentile(finite, [50, 90, 99])
        log(f"[celltype_writeback] distance percentiles μm: "
            f"p50={p50:.2f} p90={p90:.2f} p99={p99:.2f}")

    query.obs[celltype_col] = labels
    query.obs[distance_col] = distances

    # v0.2.9: when rctd_results_rds is supplied, fold per-cell RCTD
    # results + filter provenance onto query.obs. This is the xr-only
    # equivalent of what `writeback_to_raw` does on proseg_raw.h5ad —
    # xenium_ranger-only runs have no proseg_unpurified.h5ad to carry
    # these onto, so they land here on xenium_ranger.h5ad itself.
    # Back-compat: if `rctd_results_rds=None`, behavior is unchanged
    # from v0.2.8 (writes only celltype + distance).
    if rctd_results_rds is not None:
        _fold_rctd_results_onto_query(
            query=query,
            reference_h5ad=purified_h5ad,
            rctd_results_rds=Path(rctd_results_rds),
            rscript_bin=rscript_bin,
            r_lib_paths=r_lib_paths,
            rctd_obs_prefix=rctd_obs_prefix,
            rctd_input_barcodes_path=rctd_input_barcodes_path,
        )

    log(f"[celltype_writeback] writing augmented {xenium_ranger_h5ad}")
    write_kwargs = {"compression": h5ad_compression} if h5ad_compression else {}
    atomic_write_h5ad(query, xenium_ranger_h5ad, **write_kwargs)

    sentinel.parent.mkdir(parents=True, exist_ok=True)
    sentinel.write_text(
        f"celltype_writeback complete: {sample_id}\n"
        f"n_query={n_query} n_reference={reference_xy.shape[0]}\n"
    )
    log(f"[celltype_writeback] wrote sentinel {sentinel}")
    return sentinel


# v0.2.9 ---------------------------------------------------------------
# RCTD + filter-provenance writeback to xenium_ranger.h5ad's obs.


def _fold_rctd_results_onto_query(
    *,
    query,
    reference_h5ad: Path,
    rctd_results_rds: Path,
    rscript_bin: str,
    r_lib_paths,
    rctd_obs_prefix: str,
    rctd_input_barcodes_path: Path | None = None,
) -> None:
    """Fold per-cell RCTD results + filter-provenance provenance
    columns onto query.obs IN PLACE. v0.2.9 Part 1 schema (see the
    module docstring for the full semantics table + the c6755625444
    investigation).

    Writes 11 ``rctd_*`` payload columns + ``rctd_purification_status``
    (from the purified h5ad) + 3 provenance booleans:
        * ``in_rctd_input``   — cell submitted to RCTD
        * ``has_rctd_result`` — cell has a row in rctd_results_df
        * ``in_purified``     — cell in the purified h5ad

    Fill behavior (same as v0.2.9-prior):
        - categorical rctd_ cols -> "" (empty string)
        - numeric    rctd_ cols -> NaN
        - in_* bools            -> False
    """
    import anndata as ad
    import numpy as np
    import pandas as pd

    from rctd_split._internal.rctd_results import extract_rctd_results_df

    log(f"[celltype_writeback]   extracting RCTD results from {rctd_results_rds}")
    rctd_df = extract_rctd_results_df(
        rctd_results_rds=rctd_results_rds,
        rscript_bin=rscript_bin,
        r_lib_paths=list(r_lib_paths) if r_lib_paths else None,
    )

    # Reindex onto query.obs_names (as str — defensive against mixed
    # int64 / str indices, same DTYPE TRAP writeback_to_raw._fold_filter_status
    # calls out).
    query_ids = query.obs_names.astype(str)
    rctd_df.index = rctd_df.index.astype(str)
    joined = rctd_df.reindex(query_ids)

    # Decide dtype per rctd_ column:
    #   - categorical / object -> "" fill
    #   - numeric -> NaN fill (default pandas semantics)
    for col in rctd_df.columns:
        series = joined[col]
        if pd.api.types.is_numeric_dtype(series.dtype):
            query.obs[f"{rctd_obs_prefix}{col}"] = series.to_numpy(dtype=float)
        else:
            # categorical / object / string -> "" fill
            filled = series.astype(object).fillna("").astype(str)
            query.obs[f"{rctd_obs_prefix}{col}"] = filled.to_numpy()

    n_in_rctd_df = int(rctd_df.index.isin(query_ids).sum())
    log(f"[celltype_writeback]   wrote {len(rctd_df.columns)} "
        f"{rctd_obs_prefix}* columns onto query.obs "
        f"(n_cells_with_rctd_result={n_in_rctd_df}/{query.n_obs})")

    # ----- provenance booleans -------------------------------------
    # v0.2.9 Part 1 (c6755625444): the semantic split between
    # "submitted to RCTD" (test_object universe, pre-spacexr-drops)
    # and "returned a result" (rctd_results_df, post-spacexr-drops).
    # See the module docstring's end-to-end counts table.

    # has_rctd_result: cell has a row in rctd_results_df.
    # On MH8_2 v0.3 this = 255,661.
    # `query_ids.isin(...)` returns a plain numpy bool array — assigns
    # directly to adata.obs, no `.to_numpy()` needed.
    query.obs["has_rctd_result"] = query_ids.isin(rctd_df.index)

    # in_rctd_input: cell was SUBMITTED to RCTD = cell_id in
    # split_prep barcodes.tsv.gz (= test_object.rds barcodes).
    # On MH8_2 v0.3 this = 299,585. Fall-back when the caller didn't
    # supply the path: equal has_rctd_result (OVERSTATES by the
    # spacexr-dropped cells on xr-only runs; is correct on
    # all-cells-reach-RCTD runs).
    if rctd_input_barcodes_path is not None:
        rctd_input_barcodes_path = Path(rctd_input_barcodes_path)
        log(f"[celltype_writeback]   loading rctd input barcodes "
            f"{rctd_input_barcodes_path}")
        rctd_input_ids = _read_barcodes_tsv(rctd_input_barcodes_path)
        query.obs["in_rctd_input"] = query_ids.isin(rctd_input_ids)
    else:
        log("[celltype_writeback]   WARN: no rctd_input_barcodes_path "
            "supplied; in_rctd_input falls back to has_rctd_result "
            "(see TracyY123-nexus#40 c6755625444 — this OVERSTATES "
            "by the spacexr-dropped cells on xenium_ranger runs).")
        query.obs["in_rctd_input"] = query.obs["has_rctd_result"]

    # in_purified: cell survived SPLIT::purify + postprocess qc_filter
    # (i.e. is in the purified reference h5ad). On MH8_2 v0.3 this
    # = 186,080.
    log(f"[celltype_writeback]   loading purified h5ad {reference_h5ad}")
    pur_adata = ad.read_h5ad(reference_h5ad, backed="r")
    pur_ids_set = set(pur_adata.obs_names.astype(str))
    query.obs["in_purified"] = np.asarray([cid in pur_ids_set for cid in query_ids])

    # rctd_purification_status: sourced from purified.h5ad's obs
    # (populated by SPLIT::purify, preserved through mtx_to_h5ad +
    # postprocess). Only the 186,080 purified cells carry a value;
    # every other cell gets an empty string. LIMITATION, documented
    # in the module docstring + the v0.2.9 release notes: cells with
    # has_rctd_result=True but in_purified=False have LOST
    # purification_status (intermediate/split/unpurified.rds is
    # pruned by the terminal qc_report stage), so the hotfix can't
    # recover it.
    pstatus_col_name = f"{rctd_obs_prefix}purification_status"
    if "purification_status" in pur_adata.obs.columns:
        pstatus_series = pd.Series(
            pur_adata.obs["purification_status"].astype(str).to_numpy(),
            index=pur_adata.obs_names.astype(str),
        )
        joined_pstatus = (
            pstatus_series.reindex(query_ids)
            .astype(object).fillna("").astype(str)
        )
        query.obs[pstatus_col_name] = joined_pstatus.to_numpy()
        log(f"[celltype_writeback]   wrote {pstatus_col_name} from purified "
            f"h5ad's obs['purification_status'] "
            f"(n_non_empty={int((joined_pstatus != '').sum())}/{query.n_obs})")
    else:
        log(f"[celltype_writeback]   purified h5ad has no "
            f"obs['purification_status']; skipping {pstatus_col_name}")

    log(f"[celltype_writeback]   in_rctd_input=True   : "
        f"{int(query.obs['in_rctd_input'].sum())}/{query.n_obs}")
    log(f"[celltype_writeback]   has_rctd_result=True : "
        f"{int(query.obs['has_rctd_result'].sum())}/{query.n_obs}")
    log(f"[celltype_writeback]   in_purified=True     : "
        f"{int(query.obs['in_purified'].sum())}/{query.n_obs}")


def _read_barcodes_tsv(path: Path) -> set[str]:
    """Read a 10x-style ``barcodes.tsv[.gz]`` → set of cell_id strings.

    One cell_id per line; blank lines ignored. Gzip is picked by
    suffix.
    """
    import gzip

    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return {line.strip() for line in f if line.strip()}
