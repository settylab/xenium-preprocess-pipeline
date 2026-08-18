"""celltype_writeback: NN-map celltype from rctd-split purified onto xenium-preprocess xenium-ranger.

Reads:
  * ``<S>_proseg_purified.h5ad`` (this pipeline's ``mtx_to_h5ad`` output;
    source of celltype via ``.obs[reference_label_col]`` with
    ``.obs[reference_x_col]`` / ``.obs[reference_y_col]`` as reference
    centroids — ``centroid_x`` / ``centroid_y`` by default).
  * ``<S>_xenium_ranger.h5ad`` (xenium-preprocess output; query cells with
    ``.obsm['spatial']`` as ``(n, 2)`` query centroids).

Writes ``.obs[celltype_col]`` and ``.obs[distance_col]`` back on the
xenium-ranger adata via ``atomic_write_h5ad``. Pure NN — no
direct-cell-id-assign, no ``mark_unassigned`` default (per (request on
(internal issue review).

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
