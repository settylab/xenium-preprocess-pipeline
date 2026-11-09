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
these obs columns on xenium_ranger.h5ad:

  RCTD per-cell (``rctd_`` prefix, carries SPLIT post-process cols too):
    ``rctd_spot_class``, ``rctd_first_type``, ``rctd_second_type``,
    ``rctd_first_class``, ``rctd_second_class``, ``rctd_min_score``,
    ``rctd_singlet_score``, ``rctd_conv_all``, ``rctd_conv_doublet``,
    ``rctd_purification_status``, ``rctd_w1_larger_w2``,
    ``rctd_same_class`` — the subset of these actually present in the
    .rds (older runs may lack some SPLIT-emitted columns).

  Filter provenance (``in_`` prefix):
    ``in_rctd_input`` (bool): cell was in RCTD's input set.
    ``in_purified``   (bool): cell survived all the way through to the
                              post-postprocess purified h5ad.

Cells NOT in the RCTD results carry empty strings for categorical RCTD
columns, NaN for numerics, False for in_rctd_input / in_purified.

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
) -> None:
    """Fold per-cell RCTD results + in_rctd_input / in_purified bool
    provenance columns onto query.obs IN PLACE.

    Cells NOT in the RCTD results:
        - categorical rctd_ cols -> "" (empty string)
        - numeric    rctd_ cols -> NaN
        - in_rctd_input -> False
    Cells NOT in the purified reference h5ad:
        - in_purified -> False
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

    n_in_rctd = int(rctd_df.index.isin(query_ids).sum())
    log(f"[celltype_writeback]   wrote {len(rctd_df.columns)} "
        f"{rctd_obs_prefix}* columns onto query.obs "
        f"(n_cells_with_rctd={n_in_rctd}/{query.n_obs})")

    # in_rctd_input: cell is in the RCTD input set (i.e. has a row in
    # rctd_results_df).
    # `query_ids` is a pandas Index (from obs_names.astype(str)); its
    # `.isin(...)` returns a plain numpy bool array, which assigns
    # directly to adata.obs — no `.to_numpy()` needed.
    query.obs["in_rctd_input"] = query_ids.isin(rctd_df.index)

    # in_purified: cell survived SPLIT::purify + postprocess qc_filter
    # (i.e. is in the purified reference h5ad).
    log(f"[celltype_writeback]   loading purified index {reference_h5ad}")
    pur_ids = set(ad.read_h5ad(reference_h5ad, backed="r").obs_names.astype(str))
    query.obs["in_purified"] = np.asarray([cid in pur_ids for cid in query_ids])

    log(f"[celltype_writeback]   in_rctd_input=True : "
        f"{int(query.obs['in_rctd_input'].sum())}/{query.n_obs}")
    log(f"[celltype_writeback]   in_purified=True   : "
        f"{int(query.obs['in_purified'].sum())}/{query.n_obs}")
