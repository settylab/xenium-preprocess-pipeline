"""Shared nearest-neighbor celltype-mapping core.

Extracted from ``he-xenium-registration/nn_celltype_mapping.py::_run_nn``
so both the H&E-side ``nn_celltype_mapping`` stage and this package's
``celltype_writeback`` stage share one kernel with one signature:

    nn_map_celltype(*, reference_xy, reference_labels, query_xy,
                    k, algorithm, metric, tiebreak,
                    distance_threshold, unmatched_policy)
        -> {"labels": ..., "distances": ..., "reference_idx": ...}

``unmatched_policy="nearest_label"`` (default) passes the nearest
neighbor's label through even when ``distance > distance_threshold``;
``"mark_unassigned"`` puts label ``"Unassigned"`` on those rows. Per
(request on (internal issue review), ``nearest_label`` is
the default in step-4 writeback — proseg's ``original_cell_id`` bug
would spuriously flag real cells otherwise.
"""
from __future__ import annotations

from collections import Counter

import numpy as np


def _tiebreak_labels(
    labels_kn: np.ndarray,  # (n_query, k), object dtype
    dists_kn: np.ndarray,   # (n_query, k), float — unused for majority
    *,
    tiebreak: str,
) -> np.ndarray:
    """Reduce K neighbours per query to one label per query.

    ``min_dist``: take the first (already sorted by distance ascending
                  in sklearn output).
    ``majority``: mode across the K, break majority ties by nearest.
    """
    if tiebreak == "min_dist":
        return labels_kn[:, 0]
    if tiebreak == "majority":
        out = np.empty(len(labels_kn), dtype=object)
        for i in range(len(labels_kn)):
            counts = Counter(labels_kn[i].tolist())
            top_count = max(counts.values())
            for j in range(labels_kn.shape[1]):
                lab = labels_kn[i, j]
                if counts[lab] == top_count:
                    out[i] = lab
                    break
        return out
    raise ValueError(
        f"unknown tiebreak={tiebreak!r} (expected 'min_dist' or 'majority')"
    )


def nn_map_celltype(
    *,
    reference_xy: np.ndarray,
    reference_labels: np.ndarray,
    query_xy: np.ndarray,
    k: int = 1,
    algorithm: str = "auto",
    metric: str = "euclidean",
    tiebreak: str = "min_dist",
    distance_threshold: float | None = None,
    unmatched_policy: str = "nearest_label",
) -> dict:
    """Fit on reference centroids, query with query centroids.

    Returns ``{"labels": (n_query,) object, "distances": (n_query,)
    float, "reference_idx": (n_query,) int}``. ``reference_idx`` is the
    index into ``reference_labels`` of the tie-break-winning neighbour.

    Parameters
    ----------
    reference_xy, reference_labels
        ``(n_ref, 2)`` float and ``(n_ref,)`` object arrays holding
        the reference centroids and their labels.
    query_xy
        ``(n_query, 2)`` float array of query centroids.
    k, algorithm, metric, tiebreak
        Threaded through to ``sklearn.neighbors.NearestNeighbors`` /
        the tiebreaker.
    distance_threshold, unmatched_policy
        When ``distance_threshold`` is set, queries whose winning
        distance exceeds the threshold are handled per
        ``unmatched_policy``: ``"nearest_label"`` (default, keep the
        neighbor label) or ``"mark_unassigned"`` (label becomes
        ``"Unassigned"``).
    """
    from sklearn.neighbors import NearestNeighbors

    reference_xy = np.asarray(reference_xy, dtype=float)
    query_xy = np.asarray(query_xy, dtype=float)
    reference_labels = np.asarray(reference_labels, dtype=object)

    if reference_xy.ndim != 2 or reference_xy.shape[1] != 2:
        raise ValueError(
            f"reference_xy must be (n_ref, 2); got shape {reference_xy.shape}"
        )
    if query_xy.ndim != 2 or query_xy.shape[1] != 2:
        raise ValueError(
            f"query_xy must be (n_query, 2); got shape {query_xy.shape}"
        )
    if len(reference_labels) != len(reference_xy):
        raise ValueError(
            f"reference_labels length {len(reference_labels)} != "
            f"reference_xy rows {len(reference_xy)}"
        )
    if k < 1:
        raise ValueError(f"k must be >= 1; got {k}")
    if unmatched_policy not in ("nearest_label", "mark_unassigned"):
        raise ValueError(
            f"unknown unmatched_policy={unmatched_policy!r}; "
            "expected 'nearest_label' or 'mark_unassigned'"
        )

    nbrs = NearestNeighbors(n_neighbors=k, algorithm=algorithm, metric=metric)
    nbrs.fit(reference_xy)
    dists_kn, idx_kn = nbrs.kneighbors(query_xy)

    labels_kn = reference_labels[idx_kn]
    labels = _tiebreak_labels(labels_kn, dists_kn, tiebreak=tiebreak)

    if tiebreak == "min_dist" or k == 1:
        winning_idx = idx_kn[:, 0]
        winning_dist = dists_kn[:, 0]
    else:
        winning_idx = np.empty(len(labels), dtype=idx_kn.dtype)
        winning_dist = np.empty(len(labels), dtype=dists_kn.dtype)
        for i in range(len(labels)):
            for j in range(k):
                if labels_kn[i, j] == labels[i]:
                    winning_idx[i] = idx_kn[i, j]
                    winning_dist[i] = dists_kn[i, j]
                    break

    if distance_threshold is not None and unmatched_policy == "mark_unassigned":
        far = winning_dist > float(distance_threshold)
        labels = labels.copy()
        labels[far] = "Unassigned"

    return {
        "labels": labels,
        "distances": winning_dist,
        "reference_idx": winning_idx,
    }
