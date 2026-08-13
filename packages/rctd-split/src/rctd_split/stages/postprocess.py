"""Stage 6: QC + clipped-log-normalize + PCA + UMAP + Leiden.

Post-SPLIT QC + preprocessing on the purified h5ad. Runs the pipeline's
"Daniel-approach" recipe from
an internal proseg-integer SPLIT pipeline notebook
verbatim: ``sc.pp.filter_cells(min_counts=...)`` then, for each expression
matrix (purified ``.X`` and the ``layers['unpurified_counts']`` attached
by Stage 4), row-normalize → ``log(clip(X, minprop, 1))`` → PCA (all
genes) → umap-package KNN → UMAP → Leiden.

Augments the persisted ``<S>_proseg_purified.h5ad`` (under
``spatial_adata/``) in place. Sentinel is a marker file under
``intermediate/adata/``.
"""
from __future__ import annotations

import os
from pathlib import Path

from rctd_split._internal.compat import sentinel_exists
from rctd_split._internal.layout import (
    atomic_write_h5ad,
    intermediate_path,
    spatial_adata_path,
)
from rctd_split._internal.logging import log


DEFAULT_LEIDEN_RESOLUTIONS = (0.5, 0.7)


def _clipped_log_normalize(X, minprop: float, positive_x: bool):
    """Row-normalize X to sum=1, then log(clip(., minprop, 1))."""
    import numpy as np
    import scipy.sparse as sp

    if sp.issparse(X):
        arr = X.toarray().astype(np.float32)
    else:
        arr = np.asarray(X, dtype=np.float32)

    row_sums = arr.sum(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        arr = np.divide(arr, row_sums, out=np.zeros_like(arr),
                        where=row_sums > 0)
    arr = np.log(np.clip(arr, minprop, 1.0)).astype(np.float32)
    if positive_x:
        arr = arr - np.float32(np.log(minprop))
    return arr


def _pca_knn_umap(
    adata,
    layer_key: str,
    pca_key: str,
    umap_obsm_key: str,
    knn_indices_key: str,
    knn_dists_key: str,
    n_neighbors: int,
    random_state: int,
):
    import scanpy as sc
    import umap

    if pca_key == "X_pca":
        sc.pp.pca(adata, layer=layer_key, mask_var=None)
    else:
        sc.pp.pca(adata, layer=layer_key, mask_var=None, key_added=pca_key)

    knn = umap.umap_.nearest_neighbors(
        adata.obsm[pca_key],
        n_neighbors=n_neighbors,
        metric="euclidean",
        metric_kwds=None,
        angular=False,
        random_state=random_state,
    )
    Xumap = umap.UMAP(
        n_neighbors=n_neighbors, precomputed_knn=knn,
    ).fit_transform(adata.obsm[pca_key])

    adata.obsm[umap_obsm_key] = Xumap
    adata.obsm[knn_indices_key] = knn[0]
    adata.obsm[knn_dists_key] = knn[1]


def _leiden_from_knn(
    adata,
    knn_indices_key: str,
    resolution: float,
    obs_key: str,
    random_state: int,
    use_weights: bool = False,
    knn_dists_key: str | None = None,
):
    import igraph
    import leidenalg

    if knn_indices_key not in adata.obsm:
        raise SystemExit(
            f"[postprocess] {knn_indices_key!r} not found in adata.obsm; "
            f"available keys: {list(adata.obsm.keys())}"
        )
    knn_indices = adata.obsm[knn_indices_key]
    n_cells = adata.n_obs

    if use_weights:
        if knn_dists_key is None or knn_dists_key not in adata.obsm:
            raise SystemExit(
                "[postprocess] use_weights=True but knn_dists key missing"
            )
        knn_dists = adata.obsm[knn_dists_key]
    else:
        knn_dists = None

    edge_weights = {}
    for i in range(knn_indices.shape[0]):
        for k, j in enumerate(knn_indices[i, :]):
            j = int(j)
            if j < 0 or j >= n_cells or i == j:
                continue
            a, b = sorted((i, j))
            weight = 1.0 / (1.0 + float(knn_dists[i, k])) if use_weights else 1.0
            if (a, b) not in edge_weights:
                edge_weights[(a, b)] = weight
            else:
                edge_weights[(a, b)] = max(edge_weights[(a, b)], weight)

    edges = list(edge_weights.keys())
    if not edges:
        raise SystemExit("[postprocess] no edges built from KNN graph")

    G = igraph.Graph(n=n_cells, edges=edges, directed=False)
    weights = list(edge_weights.values()) if use_weights else None
    if use_weights:
        G.es["weight"] = weights

    partition = leidenalg.find_partition(
        G,
        leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=resolution,
        weights=weights,
        seed=random_state,
    )
    adata.obs[obs_key] = [str(x) for x in partition.membership]
    adata.obs[obs_key] = adata.obs[obs_key].astype("category")


def _write_summary_csv(
    csv_path: Path,
    sample_id: str,
    cfg_params: dict,
    n_cells_in: int,
    n_cells_kept: int,
) -> None:
    import pandas as pd

    row = {
        "sample_id": sample_id,
        "n_cells_in": n_cells_in,
        "n_cells_kept": n_cells_kept,
        "min_counts": cfg_params["min_counts"],
        "minprop": cfg_params["minprop"],
        "positive_x": cfg_params["positive_x"],
        "n_neighbors": cfg_params["n_neighbors"],
        "random_state": cfg_params["random_state"],
        "leiden_resolutions": ",".join(
            f"{r:g}" for r in cfg_params["leiden_resolutions"]
        ),
        "process_unpurified_layer": cfg_params["process_unpurified_layer"],
    }
    df = pd.DataFrame([row])
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = csv_path.with_suffix(csv_path.suffix + ".tmp")
    df.to_csv(tmp, index=False)
    os.replace(tmp, csv_path)


def run_postprocess(
    sample_id: str,
    run_id: str,
    output_root: Path,
    min_counts: int,
    minprop: float,
    positive_x: bool,
    n_neighbors: int,
    random_state: int,
    leiden_resolutions,
    process_unpurified_layer: bool,
    h5ad_compression: str | None,
    force_rerun: bool,
) -> Path:
    """Augment ``spatial_adata/<S>_proseg_purified.h5ad`` in place with
    QC-filtered counts + clipped-log-normalized preprocessing + PCA +
    UMAP + Leiden. Returns the sentinel path.
    """
    import anndata as ad
    import scanpy as sc

    purified_h5ad = spatial_adata_path(
        output_root, sample_id, run_id, "proseg_purified",
    )
    sentinel = intermediate_path(
        output_root, sample_id, run_id, "postprocess_sentinel",
    )
    summary_csv = intermediate_path(
        output_root, sample_id, run_id, "postprocess_summary",
    )

    leiden_resolutions = tuple(float(r) for r in leiden_resolutions)

    if sentinel_exists(sentinel, force_rerun):
        log(f"[postprocess] sentinel exists: {sentinel} — skipping "
            "(pass --force-rerun to re-run).")
        return sentinel

    if not purified_h5ad.exists():
        raise SystemExit(
            f"[postprocess] purified h5ad not found: {purified_h5ad}. "
            "Run the mtx_to_h5ad stage first."
        )

    log(f"[postprocess] reading {purified_h5ad}")
    adata = ad.read_h5ad(purified_h5ad)
    n_cells_in = adata.n_obs

    if process_unpurified_layer and "unpurified_counts" not in adata.layers:
        raise SystemExit(
            f"[postprocess] purified.h5ad missing layers['unpurified_counts'] "
            f"but process_unpurified_layer=True. Re-run mtx_to_h5ad with the "
            f"unpurified bundle available, or set process_unpurified_layer=false."
        )

    log(f"[postprocess] filter_cells(min_counts={min_counts}) "
        f"— n_cells_in={n_cells_in}")
    sc.pp.filter_cells(adata, min_counts=min_counts)
    n_cells_kept = adata.n_obs
    log(f"[postprocess]   kept {n_cells_kept}/{n_cells_in} cells")
    if n_cells_kept == 0:
        raise SystemExit(
            f"[postprocess] filter_cells dropped every cell "
            f"(min_counts={min_counts} too aggressive for this sample)"
        )

    adata.layers["counts"] = adata.X.copy()

    log(f"[postprocess] clipped-log-normalize (minprop={minprop}, "
        f"positive_x={positive_x}) on purified .X")
    adata.layers["X_clipped"] = _clipped_log_normalize(
        adata.X, minprop=minprop, positive_x=positive_x,
    )

    if process_unpurified_layer:
        log("[postprocess] clipped-log-normalize on layers['unpurified_counts']")
        adata.layers["X_clipped_unpurified"] = _clipped_log_normalize(
            adata.layers["unpurified_counts"],
            minprop=minprop, positive_x=positive_x,
        )

    log(f"[postprocess] PCA + KNN + UMAP on X_clipped "
        f"(n_neighbors={n_neighbors}, random_state={random_state})")
    _pca_knn_umap(
        adata,
        layer_key="X_clipped",
        pca_key="X_pca",
        umap_obsm_key="X_umap_clipped_norm",
        knn_indices_key="knn_indices",
        knn_dists_key="knn_dists",
        n_neighbors=n_neighbors,
        random_state=random_state,
    )

    if process_unpurified_layer:
        log("[postprocess] PCA + KNN + UMAP on X_clipped_unpurified")
        _pca_knn_umap(
            adata,
            layer_key="X_clipped_unpurified",
            pca_key="X_pca_unpurified",
            umap_obsm_key="X_umap_clipped_norm_unpurified",
            knn_indices_key="knn_indices_unpurified",
            knn_dists_key="knn_dists_unpurified",
            n_neighbors=n_neighbors,
            random_state=random_state,
        )

    for res in leiden_resolutions:
        obs_key = f"leiden_{res:g}"
        log(f"[postprocess] leiden(resolution={res}) → obs[{obs_key!r}]")
        _leiden_from_knn(
            adata,
            knn_indices_key="knn_indices",
            resolution=res,
            obs_key=obs_key,
            random_state=random_state,
        )

    log(f"[postprocess] writing augmented {purified_h5ad}")
    write_kwargs = {"compression": h5ad_compression} if h5ad_compression else {}
    atomic_write_h5ad(adata, purified_h5ad, **write_kwargs)

    _write_summary_csv(
        summary_csv,
        sample_id=sample_id,
        cfg_params={
            "min_counts": min_counts,
            "minprop": minprop,
            "positive_x": positive_x,
            "n_neighbors": n_neighbors,
            "random_state": random_state,
            "leiden_resolutions": leiden_resolutions,
            "process_unpurified_layer": process_unpurified_layer,
        },
        n_cells_in=n_cells_in,
        n_cells_kept=n_cells_kept,
    )
    log(f"[postprocess] wrote summary {summary_csv}")

    sentinel.parent.mkdir(parents=True, exist_ok=True)
    sentinel.write_text(
        f"postprocess complete: {sample_id}\n"
        f"n_cells_in={n_cells_in} n_cells_kept={n_cells_kept}\n"
    )
    log(f"[postprocess] wrote sentinel {sentinel}")
    return sentinel
