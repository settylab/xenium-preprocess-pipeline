"""Stage: enrich the xenium-ranger AnnData with an NN-mapped proseg cell_id.

Motivation (Tracy request 2026-08-11, TracyY123-nexus#26 comment
5260213932): after the big pipeline reorg, we need each xenium-ranger
cell to carry the ID of its spatially-nearest proseg cell — so
downstream analyses that live on the xenium-ranger side can join back
to proseg segmentation without a separate lookup.

This is the REVERSE of the pre-2026-08-11 direction (which added an
independent NN-computed Xenium cell_id to the proseg h5ad). Tracy
flipped it on 2026-08-11: mapping now lives on the xenium-ranger
adata, keyed by proseg_cell_id.

Strictly additive: the stage NEVER mutates any existing ``.obs``
column on the xenium-ranger adata.

Algorithm — same nearest-neighbor plumbing as
``hexenium.stages.nn_celltype_mapping`` (KDTree via sklearn's
``NearestNeighbors``), but here it's ID-mapping only (no celltype
propagation):

  1. Read the xenium-ranger h5ad written by ``xenium_ranger_to_anndata``.
  2. Pull xenium centroids from ``.obsm['spatial']`` (col 0 = x,
     col 1 = y — the stable surface that stage guarantees).
  3. Read the proseg h5ad written by ``proseg_to_anndata`` — the
     source of proseg IDs (``obs_names``) and centroids
     (``.obsm['spatial']``).
  4. Fit ``NearestNeighbors(k=1)`` on proseg centroids, query with
     xenium centroids. For each xenium cell we get: the nearest
     proseg cell_id + the distance.
  5. Add three new ``.obs`` columns to the xenium-ranger adata:

       - ``proseg_cell_id_nn`` (str) — NN-matched proseg cell id.
       - ``proseg_id_nn_distance`` (float) — NN distance in whatever
                                              frame the centroids are in.
       - ``proseg_id_nn_note`` (str)  — ``"nn_match"`` /
                                        ``"unassigned:distance_over_threshold"``.

  6. Log the NN-distance summary (min/median/p95/max) so a `grep`
     on the log immediately shows the mapping shape.
  7. Rewrite the xenium-ranger h5ad in place (default) or to a sidecar
     ``_xenium_enriched.h5ad`` (config knob).

Sentinel-file resume: separate sentinel next to the h5ad so re-runs
of the xenium-ranger stage alone don't force re-runs of this enrichment.
``--force-rerun`` nukes it.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from xenium_preprocess._internal.compat import sentinel_exists
from xenium_preprocess._internal.layout import atomic_write_h5ad
from xenium_preprocess._internal.logging import log


_NOTE_MATCH = "nn_match"
_NOTE_UNASSIGNED_THRESH = "unassigned:distance_over_threshold"


def _extract_centroids(adata, label: str) -> np.ndarray:
    """Pull (n_cells, 2) centroids out of an h5ad.

    Prefers ``.obsm['spatial']`` — the stable interface both
    ``proseg_to_anndata`` and ``xenium_ranger_to_anndata`` guarantee.
    Falls back to ``.obs['centroid_x' / 'centroid_y']`` if for some
    reason obsm is absent — with a warning.
    """
    if "spatial" in adata.obsm:
        arr = np.asarray(adata.obsm["spatial"])
        if arr.ndim != 2 or arr.shape[1] < 2:
            raise SystemExit(
                f"[enrich_xenium_id] {label}: .obsm['spatial'] has "
                f"unexpected shape {arr.shape} — expected (n_cells, >=2). "
                f"Refusing to guess the centroid layout."
            )
        return arr[:, :2].astype(float)

    log(f"[enrich_xenium_id] WARN: {label} .obsm['spatial'] absent; "
        f"falling back to .obs['centroid_x' / 'centroid_y'].")
    for c in ("centroid_x", "centroid_y"):
        if c not in adata.obs.columns:
            raise SystemExit(
                f"[enrich_xenium_id] {label} h5ad is missing obs column "
                f"{c!r} AND obsm['spatial']. Cannot recover centroids."
            )
    return np.stack(
        [
            adata.obs["centroid_x"].to_numpy(dtype=float),
            adata.obs["centroid_y"].to_numpy(dtype=float),
        ],
        axis=1,
    )


def run_enrich_xenium_id(
    sample_id: str,
    run_id: str,
    xenium_ranger_h5ad: Path,
    proseg_h5ad: Path,
    output_root: Path,
    nn_k: int,
    distance_threshold: float | None,
    nn_id_col: str,
    distance_col: str,
    note_col: str,
    write_back_h5ad: bool,
    force_rerun: bool,
) -> Path:
    """Enrich the xenium-ranger h5ad with NN-mapped proseg cell_ids.

    Returns the path to the enriched h5ad (in-place when
    ``write_back_h5ad=True``, else a sibling ``_xenium_enriched.h5ad``).
    """
    import anndata
    from sklearn.neighbors import NearestNeighbors

    if not xenium_ranger_h5ad.exists():
        raise SystemExit(
            f"[enrich_xenium_id] xenium-ranger h5ad not found: "
            f"{xenium_ranger_h5ad}. Run the xenium_ranger_to_anndata "
            f"stage first."
        )
    if not proseg_h5ad.exists():
        raise SystemExit(
            f"[enrich_xenium_id] proseg h5ad not found: {proseg_h5ad}. "
            f"Run the proseg_to_anndata stage first."
        )

    if write_back_h5ad:
        out_h5ad = xenium_ranger_h5ad
    else:
        out_h5ad = xenium_ranger_h5ad.with_name(
            xenium_ranger_h5ad.stem.replace(
                "_xenium_ranger", "_xenium_enriched"
            ) + xenium_ranger_h5ad.suffix
        )
        out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    sentinel = out_h5ad.parent / ".enrich_xenium_id_done.sentinel"

    if sentinel_exists(sentinel, force_rerun):
        log(f"[enrich_xenium_id] sentinel exists: {sentinel} — skipping "
            f"(pass --force-rerun to re-run).")
        return out_h5ad

    log(f"[enrich_xenium_id] xenium-ranger h5ad: {xenium_ranger_h5ad}")
    log(f"[enrich_xenium_id] proseg h5ad: {proseg_h5ad}")
    log(f"[enrich_xenium_id] nn_k={nn_k} distance_threshold={distance_threshold}")

    adata_xen = anndata.read_h5ad(xenium_ranger_h5ad)
    log(f"[enrich_xenium_id] xenium adata loaded: {adata_xen.n_obs} cells × "
        f"{adata_xen.n_vars} genes; obs cols: "
        f"{list(adata_xen.obs.columns)[:12]}"
        f"{'...' if len(adata_xen.obs.columns) > 12 else ''}")

    for new_col in (nn_id_col, distance_col, note_col):
        if new_col in adata_xen.obs.columns and not force_rerun:
            raise SystemExit(
                f"[enrich_xenium_id] xenium adata.obs already carries "
                f"{new_col!r} but the sentinel is missing — inconsistent "
                f"state. Pass --force-rerun to overwrite."
            )

    xenium_xy = _extract_centroids(adata_xen, label="xenium-ranger")
    log(f"[enrich_xenium_id] xenium centroids: n={len(xenium_xy)} "
        f"x=[{xenium_xy[:,0].min():.1f},{xenium_xy[:,0].max():.1f}] "
        f"y=[{xenium_xy[:,1].min():.1f},{xenium_xy[:,1].max():.1f}]")

    # Load proseg h5ad only for its centroids + obs_names (the proseg
    # cell IDs). We don't need the matrix or anything else.
    adata_proseg = anndata.read_h5ad(proseg_h5ad)
    proseg_xy = _extract_centroids(adata_proseg, label="proseg")
    proseg_ids = np.asarray(adata_proseg.obs_names, dtype=object)
    log(f"[enrich_xenium_id] proseg cells: n={len(proseg_xy)} "
        f"x=[{proseg_xy[:,0].min():.1f},{proseg_xy[:,0].max():.1f}] "
        f"y=[{proseg_xy[:,1].min():.1f},{proseg_xy[:,1].max():.1f}]")

    # Drop any proseg rows with NaN centroids — sklearn's KDTree
    # rejects NaN inputs and we'd rather skip cleanly than crash.
    finite_mask = np.isfinite(proseg_xy).all(axis=1)
    n_dropped = int((~finite_mask).sum())
    if n_dropped:
        log(f"[enrich_xenium_id] proseg: dropped {n_dropped} cells with "
            f"NaN centroid before NN fit (kept {int(finite_mask.sum())})")
        proseg_xy = proseg_xy[finite_mask]
        proseg_ids = proseg_ids[finite_mask]

    if len(proseg_xy) == 0:
        raise SystemExit(
            "[enrich_xenium_id] proseg h5ad has 0 usable centroids after "
            "NaN cleanup — nothing to NN-match against."
        )

    if nn_k != 1:
        # Majority-vote across K makes sense for celltypes (many-to-one)
        # but not for IDs (bijective). Refuse to guess.
        raise SystemExit(
            f"[enrich_xenium_id] nn_k={nn_k} is not supported — ID mapping "
            f"is 1:1 by construction. Set nn_k=1."
        )

    nbrs = NearestNeighbors(n_neighbors=1, algorithm="auto", metric="euclidean")
    nbrs.fit(proseg_xy)
    dists, idx = nbrs.kneighbors(xenium_xy)
    dists = dists[:, 0]
    idx = idx[:, 0]

    nn_id_arr = proseg_ids[idx]
    dist_arr = dists.astype(float)
    note_arr = np.full(len(xenium_xy), _NOTE_MATCH, dtype=object)

    if distance_threshold is not None:
        far_mask = dist_arr > float(distance_threshold)
        n_far = int(far_mask.sum())
        if n_far:
            nn_id_arr = nn_id_arr.astype(object).copy()
            nn_id_arr[far_mask] = None
            note_arr[far_mask] = _NOTE_UNASSIGNED_THRESH
            log(f"[enrich_xenium_id] distance_threshold={distance_threshold}: "
                f"{n_far} xenium cells beyond threshold marked "
                f"{_NOTE_UNASSIGNED_THRESH!r} "
                f"({100.0*n_far/len(xenium_xy):.1f}%)")

    # NN distance stats (unfiltered — always logged).
    log(f"[enrich_xenium_id] NN distances: min={dist_arr.min():.3f} "
        f"median={np.median(dist_arr):.3f} p95={np.percentile(dist_arr, 95):.3f} "
        f"max={dist_arr.max():.3f}")

    # Write columns into .obs (order preserved: adata_xen.obs is
    # index-aligned to the same xenium cells the centroids came from).
    adata_xen.obs[nn_id_col] = nn_id_arr
    adata_xen.obs[distance_col] = dist_arr
    adata_xen.obs[note_col] = note_arr

    atomic_write_h5ad(adata_xen, out_h5ad)
    log(f"[enrich_xenium_id] wrote {out_h5ad}")
    sentinel.write_text("ok\n")
    log(f"[enrich_xenium_id] sentinel -> {sentinel.name}")
    return out_h5ad
