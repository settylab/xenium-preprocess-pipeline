"""Stage 2: raw proseg AnnData → preprocessed unpurified AnnData + UMAP plots.

Adapted from
    the internal preprocessing reference notebook
using the helpers imported from the internal source
    the internal Xenium data-processing routines

The reference notebook's imports live at:
    from XeniumDataProcessing import (
        remove_control_probes, daniel_approach_hvg_pca,
        leiden_clustering_via_knn_graph, read_marker_genes_group_names,
        find_potential_cell_types, generate_dict_for_celltype_annotation_global,
        map_celltypes, detect_low_count, read_marker_genes,
        calculate_module_score, calculate_cluster_based_mean_value,
    )

We port these into standalone helpers below so the pipeline does NOT
depend on that internal package being on sys.path in the runtime env.

Dual-matrix support (user request 2026-07-10, (internal issue review):
when `dual_matrix_mode=True` AND both `expected_counts` + `maxpost_counts`
layers are present in the raw h5ad, this stage runs TWO independent
preprocessing passes — one on each layer — and writes layer-suffixed
outputs. The two passes share only the QC filter (control-probe drop +
per-cell min-counts); PCA, KNN, UMAP, Leiden, and cell-type annotation
run separately on each layer, so the two label sets may disagree —
inspecting that divergence is the whole point.

The pipeline stops BEFORE the manual `label_clusters` fixes in the
reference notebook (those are curation, and areuser-driven). The
output is deliberately named `adata_unpurified.h5ad`.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from xenium_preprocess._internal.compat import sentinel_exists
from xenium_preprocess._internal.layout import atomic_write_h5ad
from xenium_preprocess._internal.logging import log


# ---------------------------------------------------------------------
# Ported helpers — verbatim (behaviour) copies of the functions in
# XeniumDataProcessing/QualityControl.py, Preprocessing.py and
# CellAnnotation.py. Docstrings pared down; the algorithms are unchanged.
# The `suffix` parameter (empty string for the single-pass case) lets
# the two dual-matrix passes coexist without stepping on each other's
# obsm / layers / obs keys.
# ---------------------------------------------------------------------

def remove_control_probes(adata):
    """Drop genes whose name starts with 'Neg' or 'Unassigned'.

    Adapted from XeniumDataProcessing.QualityControl.remove_control_probes.
    """
    adata.var["gene_symbol"] = adata.var.index.tolist()
    control_probes = [
        gene
        for gene in adata.var["gene_symbol"]
        if gene.startswith("Neg") or gene.startswith("Unassigned")
    ]
    keep_genes = list(set(adata.var["gene_symbol"]) - set(control_probes))
    return adata[:, adata.var["gene_symbol"].isin(keep_genes)].copy()


def _matrix_from_source(adata, source_layer: str | None):
    """Return a dense (n_cells × n_genes) float32 view of the source.

    - `source_layer=None` → use `adata.X` (the pre-dual-mode default).
    - `source_layer=<name>` → use `adata.layers[<name>]`.
    """
    if source_layer is None:
        source = adata.X
    else:
        if source_layer not in adata.layers:
            raise KeyError(
                f"source layer {source_layer!r} not present in adata.layers. "
                f"Available layers: {list(adata.layers.keys())}"
            )
        source = adata.layers[source_layer]
    if hasattr(source, "toarray"):
        return source.toarray().astype(np.float32)
    return np.asarray(source, dtype=np.float32)


def daniel_approach_hvg_pca(
    adata,
    umap_layer_name: str,
    positive_X: bool = True,
    n_neighbors: int = 15,
    min_prop: float = 1e-3,
    random_state: int = 0,
    source_layer: str | None = None,
    suffix: str = "",
):
    """Clipped-log normalization → PCA → KNN → UMAP.

    Adapted from XeniumDataProcessing.Preprocessing.daniel_approach_hvg_pca.
    X is clipped-normalized and log-transformed, saved as `X_clipped`
    layer. When positive_X is true, the log matrix is shifted by
    ``-log(min_prop)`` and any residual negatives are clipped to zero
    so ``X_clipped`` is non-negative (user request 2026-07-09; matches
    the make_X_positive(X, MINPROP) recipe). PCA runs on the clipped
    layer. KNN is built by umap-learn's nearest_neighbors and UMAP runs
    on X_pca.

    `source_layer` picks the input matrix — `None` for `.X` (pre-dual-mode
    default), a name for `adata.layers[<name>]`. `suffix` disambiguates
    the output keys (`X_pca{suffix}`, `X_clipped{suffix}`, etc.) so two
    passes on two layers coexist.
    """
    import scanpy as sc
    import umap

    x_clipped_key = f"X_clipped{suffix}"
    pca_key = f"X_pca{suffix}"
    knn_indices_key = f"knn_indices{suffix}"
    knn_dists_key = f"knn_dists{suffix}"
    umap_key = f"X_{umap_layer_name}{suffix}"

    src = _matrix_from_source(adata, source_layer)

    # clipped normalization
    X = src.copy()
    X /= X.sum(axis=1, keepdims=True)
    X = np.log(np.clip(X, min_prop, 1.0))

    if positive_X:
        # make_X_positive(X, MINPROP): shift by -log(MINPROP) then clip
        # residual negatives (float-roundoff safety net) to zero.
        log_minprop = np.log(min_prop)
        X = X - log_minprop
        X[X < 0] = 0

    adata.layers[x_clipped_key] = X

    # sc.pp.pca writes to adata.obsm["X_pca"] + adata.varm["PCs"] +
    # adata.uns["pca"] with no key_added support. Compute-then-rename
    # so the second dual-mode pass doesn't obliterate the first's PC
    # loadings and variance-ratio metadata.
    sc.pp.pca(adata, layer=x_clipped_key, mask_var=None)
    if pca_key != "X_pca":
        adata.obsm[pca_key] = adata.obsm.pop("X_pca")
        if "PCs" in adata.varm:
            adata.varm[f"PCs{suffix}"] = adata.varm.pop("PCs")
        if "pca" in adata.uns:
            adata.uns[f"pca{suffix}"] = adata.uns.pop("pca")

    knn = umap.umap_.nearest_neighbors(
        adata.obsm[pca_key],
        n_neighbors=n_neighbors,
        metric="euclidean",
        metric_kwds=None,
        angular=False,
        random_state=random_state,
    )
    Xumap = umap.UMAP(
        n_neighbors=n_neighbors, precomputed_knn=knn
    ).fit_transform(adata.obsm[pca_key])

    adata.obsm[umap_key] = Xumap
    adata.obsm[knn_indices_key] = knn[0]
    adata.obsm[knn_dists_key] = knn[1]


def leiden_clustering_via_knn_graph(
    adata,
    resolution: float = 0.5,
    suffix: str = "",
    random_state: int = 0,
    use_weights: bool = False,
    output_key: str | None = None,
):
    """Leiden clustering on the KNN graph saved by daniel_approach_hvg_pca.

    Adapted from XeniumDataProcessing.Preprocessing.leiden_clustering_via_knn_graph.
    """
    import igraph
    import leidenalg

    knn_indices_key = "knn_indices" + suffix
    knn_dists_key = "knn_dists" + suffix

    if knn_indices_key not in adata.obsm:
        raise ValueError(
            f"{knn_indices_key} not found in adata.obsm. "
            f"Available obsm keys are: {list(adata.obsm.keys())}"
        )
    knn_indices = adata.obsm[knn_indices_key]

    if use_weights:
        if knn_dists_key not in adata.obsm:
            raise ValueError(
                f"use_weights=True, but {knn_dists_key} was not found in adata.obsm."
            )
        knn_dists = adata.obsm[knn_dists_key]
    else:
        knn_dists = None

    edge_weights: dict[tuple[int, int], float] = {}
    n_cells = adata.n_obs
    for i in range(knn_indices.shape[0]):
        for k, j in enumerate(knn_indices[i, :]):
            j = int(j)
            if j < 0 or j >= n_cells or i == j:
                continue
            a, b = sorted((i, j))
            if use_weights:
                dist = float(knn_dists[i, k])
                weight = 1.0 / (1.0 + dist)
            else:
                weight = 1.0
            if (a, b) not in edge_weights:
                edge_weights[(a, b)] = weight
            else:
                edge_weights[(a, b)] = max(edge_weights[(a, b)], weight)

    edges = list(edge_weights.keys())
    if len(edges) == 0:
        raise ValueError("No edges were created from the KNN graph.")

    G = igraph.Graph(n=n_cells, edges=edges, directed=False)
    weights = list(edge_weights.values()) if use_weights else None
    if use_weights:
        G.es["weight"] = weights

    log(f"[preprocess] Leiden KNN graph ({suffix or 'single'}): "
        f"n_cells={n_cells} n_edges={len(edges)}")

    partition = leidenalg.find_partition(
        G,
        leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=resolution,
        weights=weights,
        seed=random_state,
    )

    if output_key is None:
        output_key = f"leiden_{resolution}{suffix}"
    adata.obs[output_key] = [str(x) for x in partition.membership]
    adata.obs[output_key] = adata.obs[output_key].astype("category")

    log(f"[preprocess] Leiden clusters saved to adata.obs[{output_key!r}] "
        f"n_clusters={adata.obs[output_key].nunique()}")
    return adata


def read_marker_genes_group_names(json_file_path):
    """Adapted from XeniumDataProcessing.CellAnnotation.read_marker_genes_group_names."""
    with open(json_file_path) as f:
        data = json.load(f)
    return list(data.keys()), data


def read_marker_genes(file_path, celltype_name):
    """Adapted from XeniumDataProcessing.CellAnnotation.read_marker_genes."""
    with open(file_path) as f:
        data = json.load(f)
    for k in data.keys():
        if celltype_name in k:
            return data[k]
    raise KeyError(
        f"Celltype '{celltype_name}' not found in JSON keys. "
        f"Available keys are: {', '.join(data.keys())}"
    )


def calculate_cluster_based_mean_value(adata, obs_column: str, leiden_key: str):
    """Adapted from XeniumDataProcessing.CellAnnotation.calculate_cluster_based_mean_value."""
    cluster_mean = (
        adata.obs.groupby(leiden_key)[obs_column]
        .mean()
        .reset_index(name=f"mean_{obs_column}")
    )
    return cluster_mean


def calculate_module_score(
    adata, gene_list, gene_list_name: str, leiden_key: str, x_layer: str = "X_clipped"
):
    """Adapted from XeniumDataProcessing.CellAnnotation.calculate_module_score."""
    import scanpy as sc

    # HDF5 group keys can't contain '/'; the marker JSON may declare keys
    # like 'B/Plasma_marker' or 'T/NK_marker'. Sanitize before it lands
    # as an adata.obs column, otherwise anndata's write_h5ad blows up at
    # end-of-stage with "Forward slashes are not allowed in keys."
    safe_name = gene_list_name.replace("/", "_")
    score_name = f"module_score_{safe_name}"

    sc.tl.score_genes(
        adata,
        gene_list,
        score_name=score_name,
        ctrl_size=50,
        use_raw=False,
        layer=x_layer,
    )
    return calculate_cluster_based_mean_value(adata, score_name, leiden_key)


def find_potential_cell_types(
    adata, celltype_group_names, data, leiden_key: str, x_layer: str = "X_clipped"
):
    """Adapted from XeniumDataProcessing.CellAnnotation.find_potential_cell_types."""
    clusters = sorted(adata.obs[leiden_key].unique())
    final_module_scores = pd.DataFrame(clusters, columns=[leiden_key])

    for group_name in celltype_group_names:
        marker_genes = data[group_name]
        module_score_df = calculate_module_score(
            adata, marker_genes, gene_list_name=group_name,
            leiden_key=leiden_key, x_layer=x_layer,
        )
        final_module_scores = pd.merge(
            final_module_scores, module_score_df, on=leiden_key, how="right"
        )

    final_module_scores.index = final_module_scores[leiden_key]
    del final_module_scores[leiden_key]
    return final_module_scores


def detect_low_count(adata, leiden_key: str):
    """Adapted from XeniumDataProcessing.CellAnnotation.detect_low_count.

    Returns `(low_count_list, average_count)`. The reference function
    additionally calls sc.pl.embedding for interactive inspection; that
    plotting side-effect is omitted here (the UMAP is written separately).
    """
    average_count = calculate_cluster_based_mean_value(
        adata, "log1p_total_counts", leiden_key
    )
    max_mean = np.max(average_count.iloc[:, 1])
    min_mean = np.min(average_count.iloc[:, 1])

    if (abs(max_mean) - abs(min_mean) >= 2) and (abs(min_mean) <= 4):
        low_cluster = average_count.loc[
            average_count["mean_log1p_total_counts"] == min_mean,
            average_count.columns[0],
        ].iloc[0]
        log(f"[preprocess] low-count cluster candidate: {low_cluster}")
        return {"potential low-count cluster": low_cluster}, average_count
    return [], average_count


def generate_dict_for_celltype_annotation_global(
    final_module_scores, low_count_list, low_count_name: str, threshold: float = 0.1
):
    """Adapted verbatim from
    XeniumDataProcessing.CellAnnotation.generate_dict_for_celltype_annotation_global.

    NOTE: preserves the reference's `low_count_list[i]` indexing shape
    (works when low_count_list is a dict indexed by str, which is what
    detect_low_count returns for a hit; empty otherwise).
    """
    max_cols = final_module_scores.iloc[:, 0:].idxmax(axis=1)
    max_values = final_module_scores.iloc[:, 0:].max(axis=1).tolist()
    max_cols_list = max_cols.tolist()
    celltypes = [col.split("_")[-2] for col in max_cols_list]

    for i, score in enumerate(max_values):
        if score <= threshold:
            celltypes[i] = "unknown_maybe_tumor"

    max_row = final_module_scores.iloc[:, 0:].idxmax(axis=0).tolist()
    index_to_int = [final_module_scores.index.get_loc(label) for label in max_row]
    counted = Counter(index_to_int)
    repeats = [idx for idx in index_to_int if counted[idx] > 1]
    module_cts = [c.split("_")[-2] for c in final_module_scores.columns[0:]]

    for j, cluster_idx in enumerate(index_to_int):
        if cluster_idx in repeats:
            js = [jj for jj, rid in enumerate(index_to_int) if rid == cluster_idx]
            combined = "_".join(module_cts[jj] for jj in js)
            celltypes[cluster_idx] = combined

    cluster_ids = final_module_scores.index.tolist()
    dict_celltype = {cid: celltypes[i] for i, cid in enumerate(cluster_ids)}
    if len(low_count_list) > 0:
        for i in low_count_list:
            dict_celltype[low_count_list[i]] = low_count_name

    return dict_celltype


def map_celltypes(adata, dict_celltype, celltype_col_name: str, leiden_key: str):
    """Adapted from XeniumDataProcessing.CellAnnotation.map_celltypes."""
    adata.obs[celltype_col_name] = adata.obs[leiden_key].map(dict_celltype)


# ---------------------------------------------------------------------
# Per-pass driver: run one full PCA→UMAP→Leiden→annotate cycle against
# a chosen source layer, tagging all outputs with `suffix`. Used once
# in single-matrix mode (suffix="") and twice in dual-matrix mode.
# ---------------------------------------------------------------------

def _run_one_pass(
    adata,
    *,
    source_layer: str | None,
    suffix: str,
    plots_dir: Path,
    min_prop: float,
    positive_x: bool,
    n_neighbors: int,
    random_state: int,
    leiden_resolution: float,
    use_leiden_weights: bool,
    global_non_tumor_json: Path | None,
    global_tumor_json: Path | None,
    tumor_type: str | None,
    threshold_global: float,
    low_count_label: str,
    global_level1_celltype_col: str,
    umap_dpi: int,
    umap_figsize,
) -> None:
    """Run one PCA→KNN→UMAP→Leiden→(optional annotate) pass.

    All outputs land under keys suffixed with `suffix` so two passes
    on two layers don't overwrite each other.
    """
    import matplotlib.pyplot as plt
    import scanpy as sc

    pass_label = suffix.lstrip("_") if suffix else "single"
    log(f"[preprocess] pass={pass_label!r} source_layer={source_layer!r} "
        f"suffix={suffix!r}")

    # --- Normalization + PCA + KNN + UMAP
    daniel_approach_hvg_pca(
        adata,
        umap_layer_name="umap",
        positive_X=positive_x,
        n_neighbors=n_neighbors,
        min_prop=min_prop,
        random_state=random_state,
        source_layer=source_layer,
        suffix=suffix,
    )

    # --- Leiden clustering
    leiden_clustering_via_knn_graph(
        adata,
        resolution=leiden_resolution,
        random_state=random_state,
        use_weights=use_leiden_weights,
        suffix=suffix,
    )

    # --- Leiden UMAP plot
    leiden_key = f"leiden_{leiden_resolution}{suffix}"
    umap_key = f"X_umap{suffix}"
    fig, ax = plt.subplots(figsize=tuple(umap_figsize), dpi=umap_dpi)
    sc.pl.embedding(adata, umap_key.removeprefix("X_"), color=leiden_key,
                    ax=ax, show=False)
    plot_name = "umap_leiden" + (f"_{pass_label}" if suffix else "") + ".png"
    leiden_umap_path = plots_dir / plot_name
    fig.savefig(leiden_umap_path, dpi=umap_dpi, bbox_inches="tight")
    plt.close(fig)
    log(f"[preprocess] wrote {leiden_umap_path}")

    # --- Optional celltype inference from marker JSONs
    if global_non_tumor_json is None:
        log(f"[preprocess] pass={pass_label!r}: global_non_tumor_json is null — "
            f"celltype inference skipped.")
        return

    log(f"[preprocess] pass={pass_label!r}: running celltype inference against "
        f"{global_non_tumor_json}")
    celltype_group_names, data_global = read_marker_genes_group_names(
        global_non_tumor_json
    )
    x_clipped_key = f"X_clipped{suffix}"

    low_count_list, _average_count = detect_low_count(adata, leiden_key)
    final_module_scores_global = find_potential_cell_types(
        adata, celltype_group_names, data_global, leiden_key, x_layer=x_clipped_key
    )
    dict_celltype_global = generate_dict_for_celltype_annotation_global(
        final_module_scores_global, low_count_list, low_count_label, threshold_global
    )
    celltype_col = f"{global_level1_celltype_col}{suffix}"
    map_celltypes(adata, dict_celltype_global, celltype_col, leiden_key)

    # Celltype UMAP plot
    fig, ax = plt.subplots(figsize=tuple(umap_figsize), dpi=umap_dpi)
    sc.pl.embedding(adata, umap_key.removeprefix("X_"), color=celltype_col,
                    ax=ax, show=False)
    if suffix:
        celltype_plot_name = f"umap_celltype_{pass_label}.png"
    else:
        # Preserve the historical single-pass plot name so anyone
        # grepping for it in older notebooks still finds it.
        celltype_plot_name = "umap_global_level1_celltype.png"
    celltype_umap_path = plots_dir / celltype_plot_name
    fig.savefig(celltype_umap_path, dpi=umap_dpi, bbox_inches="tight")
    plt.close(fig)
    log(f"[preprocess] wrote {celltype_umap_path}")

    # Optional tumor-marker overlay (reference notebook plots the
    # top-N cancer-marker genes on the UMAP for QC).
    if global_tumor_json is not None and tumor_type is not None:
        try:
            cancer_marker = read_marker_genes(global_tumor_json, tumor_type)
            fig = sc.pl.embedding(
                adata,
                umap_key.removeprefix("X_"),
                layer=x_clipped_key,
                color=cancer_marker + [celltype_col],
                title=[f"{g} expression" for g in cancer_marker]
                + [celltype_col],
                ncols=4,
                vmax="p99",
                cmap="Blues",
                show=False,
                return_fig=True,
            )
            tumor_plot_name = f"umap_tumor_markers_{tumor_type}"
            if suffix:
                tumor_plot_name += f"_{pass_label}"
            tumor_umap_path = plots_dir / f"{tumor_plot_name}.png"
            fig.savefig(tumor_umap_path, dpi=umap_dpi, bbox_inches="tight")
            plt.close(fig)
            log(f"[preprocess] wrote {tumor_umap_path}")
        except KeyError as e:
            log(f"[preprocess] WARN: tumor-marker overlay skipped ({e}).")


# ---------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------

def run_preprocess(
    sample_id: str,
    raw_h5ad: Path,
    output_root: Path,
    qc_percent_top,
    min_counts_cell: int,
    min_prop: float,
    positive_x: bool,
    n_neighbors: int,
    random_state: int,
    leiden_resolution: float,
    use_leiden_weights: bool,
    global_non_tumor_json: Path | None,
    global_tumor_json: Path | None,
    tumor_type: str | None,
    threshold_global: float,
    low_count_label: str,
    global_level1_celltype_col: str,
    umap_dpi: int,
    umap_figsize,
    force_rerun: bool,
    legacy_symlinks: bool = True,   # accepted for backward compat; unused (no more legacy symlinks)
    dual_matrix_mode: bool = True,
    run_id: str | None = None,       # accepted for uniform signature; unused (preprocess is deprecated)
) -> Path:
    """Run the preprocessing stage. Returns the h5ad path.

    `dual_matrix_mode=True` runs two independent PCA→UMAP→Leiden→annotate
    passes — one on `layers["expected_counts"]`, one on
    `layers["maxpost_counts"]` — with per-pass suffixes so results
    don't overlap. Falls back to a single-pass on `.X` when either
    layer is missing (with a WARN log).
    """
    import matplotlib
    matplotlib.use("Agg")
    import scanpy as sc

    # `preprocess` is DEPRECATED — not in DEFAULT_STAGES. When opted-in
    # via `--stages ... preprocess ...` it writes to a `legacy_preprocess/`
    # subfolder next to the raw h5ad so nothing collides with the locked
    # `spatial_adata/` / `rctd/` shape.
    legacy_dir = raw_h5ad.parent.parent / "legacy_preprocess"
    legacy_dir.mkdir(parents=True, exist_ok=True)
    out_h5ad = legacy_dir / f"{sample_id}_preprocessed.h5ad"
    plots_dir = legacy_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    if sentinel_exists(out_h5ad, force_rerun):
        log(f"[preprocess] sentinel exists: {out_h5ad} — skipping "
            f"(pass --force-rerun to re-run).")
        return out_h5ad

    log(f"[preprocess] reading {raw_h5ad}")
    adata = sc.read_h5ad(raw_h5ad)
    log(f"[preprocess] initial adata: {adata.n_obs} cells × {adata.n_vars} genes; "
        f"layers={list(adata.layers.keys())}")

    # --- QC (mirrors the notebook's calculate_qc_metrics + filter_cells)
    adata = remove_control_probes(adata)
    sc.pp.calculate_qc_metrics(adata, percent_top=tuple(qc_percent_top), inplace=True)
    sc.pp.filter_cells(adata, min_counts=min_counts_cell)
    log(f"[preprocess] after QC filter: {adata.n_obs} cells × {adata.n_vars} genes")

    # Decide single vs dual pass. Dual requires BOTH layers present; a
    # missing layer falls back to single-pass on `.X` with a WARN.
    has_expected = "expected_counts" in adata.layers
    has_maxpost = "maxpost_counts" in adata.layers
    if dual_matrix_mode and has_expected and has_maxpost:
        log("[preprocess] dual_matrix_mode: running TWO passes — expected + maxpost.")
        _run_one_pass(
            adata,
            source_layer="expected_counts",
            suffix="_expected",
            plots_dir=plots_dir,
            min_prop=min_prop,
            positive_x=positive_x,
            n_neighbors=n_neighbors,
            random_state=random_state,
            leiden_resolution=leiden_resolution,
            use_leiden_weights=use_leiden_weights,
            global_non_tumor_json=global_non_tumor_json,
            global_tumor_json=global_tumor_json,
            tumor_type=tumor_type,
            threshold_global=threshold_global,
            low_count_label=low_count_label,
            global_level1_celltype_col=global_level1_celltype_col,
            umap_dpi=umap_dpi,
            umap_figsize=umap_figsize,
        )
        _run_one_pass(
            adata,
            source_layer="maxpost_counts",
            suffix="_maxpost",
            plots_dir=plots_dir,
            min_prop=min_prop,
            positive_x=positive_x,
            n_neighbors=n_neighbors,
            random_state=random_state,
            leiden_resolution=leiden_resolution,
            use_leiden_weights=use_leiden_weights,
            global_non_tumor_json=global_non_tumor_json,
            global_tumor_json=global_tumor_json,
            tumor_type=tumor_type,
            threshold_global=threshold_global,
            low_count_label=low_count_label,
            global_level1_celltype_col=global_level1_celltype_col,
            umap_dpi=umap_dpi,
            umap_figsize=umap_figsize,
        )
    else:
        if dual_matrix_mode and not (has_expected and has_maxpost):
            log(f"[preprocess] WARN: dual_matrix_mode=true but layers missing "
                f"(expected_counts={has_expected}, maxpost_counts={has_maxpost}). "
                f"Falling back to single-pass on `.X`. Rerun proseg_to_anndata "
                f"with maxpost_matrix_glob set to populate both layers.")
        else:
            log("[preprocess] single-matrix mode: one pass on `.X`.")
        _run_one_pass(
            adata,
            source_layer=None,
            suffix="",
            plots_dir=plots_dir,
            min_prop=min_prop,
            positive_x=positive_x,
            n_neighbors=n_neighbors,
            random_state=random_state,
            leiden_resolution=leiden_resolution,
            use_leiden_weights=use_leiden_weights,
            global_non_tumor_json=global_non_tumor_json,
            global_tumor_json=global_tumor_json,
            tumor_type=tumor_type,
            threshold_global=threshold_global,
            low_count_label=low_count_label,
            global_level1_celltype_col=global_level1_celltype_col,
            umap_dpi=umap_dpi,
            umap_figsize=umap_figsize,
        )

    atomic_write_h5ad(adata, out_h5ad)
    log(f"[preprocess] wrote {out_h5ad}")
    return out_h5ad
