# Methods

`xenium-preprocess` is a stage-orchestrated wrapper around the internal SPLIT/Proseg workflow's proseg → SPLIT → RCTD preprocessing workflow. Behaviour is a verbatim port of the reference scripts and notebooks; this document points each stage back at its source.

## Stage 1: `proseg_to_anndata`

**Source:** the internal proseg-to-anndata reference script.

Reads a proseg count matrix (rows = cells, columns = genes) and cell metadata table, auto-detecting `.parquet` vs `.csv` from the extension. Centroids (`centroid_x`, `centroid_y`) are stacked into `adata.obsm['spatial']`; metadata rows become `adata.obs` indexed by `cell` (renamed to `cell_d` to sidestep anndata's index-name-collision rule); counts are stored as a CSR sparse matrix (`scipy.sparse.csr_matrix`). No filtering, no normalization — this stage is a straight port of the original 40-line script.

### Dual-matrix support (2026-07-10)

Proseg emits two count matrices:

- **`expected-counts.csv.gz`** — continuous, per-cell expected transcript counts.
- **`maxpost_counts.csv.gz`** — integer, maximum-a-posteriori transcript assignments.

When `proseg_to_anndata.maxpost_matrix_glob` is non-null (default `"maxpost_counts*"`), both are loaded and stored side-by-side:

| Slot | Contents |
|---|---|
| `adata.X` | picked by `x_source` — default `maxpost_counts` (integer MAP counts, 2026-07-23 default flip). Pass `x_source: expected_counts` to opt back into the pre-2026-07-23 fractional-expected layout. |
| `adata.layers["expected_counts"]` | expected-counts (CSR, float — the continuous per-cell expected transcript matrix, kept as an auxiliary layer for downstream analyses that want the fractional matrix). |
| `adata.layers["maxpost_counts"]` | maxpost integer counts (CSR, float dtype — split_prep casts to int32 on export). |

The two matrices are validated for exact-match cells × genes and identical gene column order. A shape/order mismatch fails loud — a silent gene shuffle would corrupt every downstream layer-scoped result. Setting `maxpost_matrix_glob: null` skips the maxpost load and reverts to the original single-matrix behaviour; with the default `x_source: maxpost_counts` in that setup, `.X` silently falls back to `expected_counts` with a WARN log (preserves back-compat for pre-dual-mode configs).

**Which layer `.X` mirrors matters for QC** (`sc.pp.filter_cells`, `sc.pp.calculate_qc_metrics`) and for `preprocess`'s single-pass mode — both read `.X`. Under the 2026-07-23 default, QC and single-pass normalization see integer counts (the scanpy-conventional shape). Dual-pass mode reads `.layers[...]` directly and is unaffected.

## Stage 2: `preprocess`

**Sources:**
- the internal preprocessing reference notebook — the cell-by-cell reference workflow.
- the internal Xenium data-processing routines — the internal source whose routines the notebook imports.

### Design choice: vendored, not imported

The reference notebook does `from XeniumDataProcessing import (…)` — that package lives outside the pipeline and its version is not pinned. Rather than add a fragile `sys.path.insert` in the runtime env (option (a) in the task spec) or vendor a copy of the package under `_vendored/` (option (b)), the eight routines actually used by the preprocessing workflow have been ported into `xenium_preprocess.stages.preprocess` as standalone functions. Each carries an `# Adapted from …` docstring pointing back at the source file. This keeps the pipeline self-contained and lets it evolve without breaking the upstream package.

The eight ported routines are:

| Ported function | Source |
|---|---|
| `remove_control_probes` | `XeniumDataProcessing/QualityControl.py` |
| `daniel_approach_hvg_pca` | `XeniumDataProcessing/Preprocessing.py` |
| `leiden_clustering_via_knn_graph` | `XeniumDataProcessing/Preprocessing.py` |
| `read_marker_genes_group_names` | `XeniumDataProcessing/CellAnnotation.py` |
| `read_marker_genes` | `XeniumDataProcessing/CellAnnotation.py` |
| `calculate_cluster_based_mean_value` | `XeniumDataProcessing/CellAnnotation.py` |
| `calculate_module_score` | `XeniumDataProcessing/CellAnnotation.py` |
| `find_potential_cell_types` | `XeniumDataProcessing/CellAnnotation.py` |
| `detect_low_count` | `XeniumDataProcessing/CellAnnotation.py` |
| `generate_dict_for_celltype_annotation_global` | `XeniumDataProcessing/CellAnnotation.py` |
| `map_celltypes` | `XeniumDataProcessing/CellAnnotation.py` |

### Steps

1. **`remove_control_probes(adata)`** — drop genes whose name starts with `Neg` or `Unassigned`.
2. **`sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150))`** followed by **`sc.pp.filter_cells(adata, min_counts=10)`** — this is the ONE step shared across both dual-matrix passes.
3. **`daniel_approach_hvg_pca`** — clipped-log normalization
   ```
   X /= X.sum(axis=1, keepdims=True)
   X = np.log(np.clip(X, min_prop, 1.0))       # min_prop = 1e-3
   ```
   stored as the `X_clipped{suffix}` layer; PCA on that layer; `umap.umap_.nearest_neighbors(X_pca{suffix}, n_neighbors=15)` for the KNN; `umap.UMAP(precomputed_knn=…).fit_transform(X_pca{suffix})` for the UMAP.
4. **`leiden_clustering_via_knn_graph`** — build an unweighted igraph from `obsm['knn_indices{suffix}']`; run `leidenalg.find_partition(G, RBConfigurationVertexPartition, resolution_parameter=0.4, seed=0)`.
5. **UMAP plot** — `sc.pl.embedding` coloured by `leiden_0.4{suffix}`, saved as `plots/umap_leiden{_pass}.png`.
6. **Optional celltype inference** (fires when `--global-non-tumor-json` is set): as before, but the inference runs against each pass's own `leiden_0.4{suffix}` partition and `X_clipped{suffix}` layer, and writes into `adata.obs['global_level1_celltype{suffix}']`.

### Dual-matrix mode (2026-07-10)

When `preprocess.dual_matrix_mode: true` (default) AND both `expected_counts` + `maxpost_counts` layers are present in the raw h5ad, steps 3–6 run TWICE — once per layer — with per-pass suffixes:

| Pass | `source_layer` | `suffix` | Key outputs |
|---|---|---|---|
| Expected | `expected_counts` | `_expected` | `X_pca_expected`, `X_umap_expected`, `X_clipped_expected`, `leiden_0.4_expected`, `global_level1_celltype_expected`, `plots/umap_leiden_expected.png`, `plots/umap_celltype_expected.png` |
| Maxpost | `maxpost_counts` | `_maxpost` | `X_pca_maxpost`, `X_umap_maxpost`, `X_clipped_maxpost`, `leiden_0.4_maxpost`, `global_level1_celltype_maxpost`, `plots/umap_leiden_maxpost.png`, `plots/umap_celltype_maxpost.png` |

The two passes are fully independent — annotation runs against each pass's own Leiden partition, so cell-type labels may disagree between the two. That divergence is exactly what the user asked to be able to inspect (internal issue review).

Setting `preprocess.dual_matrix_mode: false` reverts to a single pass on `.X`, writing unsuffixed keys (`X_pca`, `X_umap`, `leiden_0.4`, `global_level1_celltype`, `plots/umap_leiden.png`, `plots/umap_global_level1_celltype.png`) — this preserves the pre-2026-07-10 behaviour. Combined with `proseg_to_anndata.maxpost_matrix_glob: null`, the pipeline reverts to the exact single-matrix behaviour it originally shipped with.

The pipeline stops *before* the manual `label_clusters` fixes in the reference notebook (`label_clusters(adata, "unknown_maybe_tumor", "Tumor", …)`, `label_clusters(adata, "9", "NonTumor", …)`, etc.). Those are curation decisions the user makes per sample; the pipeline output is deliberately named `adata_unpurified.h5ad` to reflect this.

## Stage 3: `split_prep`

**Source:** the internal SPLIT preparation notebook.

Materialises the export layer (`adata.layers['maxpost_counts']` by default — the integer maxpost counts written by `proseg_to_anndata`; falls back to `adata.X` when the requested layer is missing, preserving pre-2026-07-10 back-compat), casts to `int32`, transposes to `genes × cells`, and writes the 10x-style triple plus two sidecar files:

- `<sample>{suffix}_counts.mtx.gz` — MatrixMarket, gzipped, `scipy.io.mmwrite`.
- `<sample>{suffix}_features.tsv.gz` — one gene name per line (single-column form; matches Seurat `ReadMtx(feature.column=1)`).
- `<sample>{suffix}_barcodes.tsv.gz` — one cell id per line.
- `<sample>{suffix}_metadata.csv` — `adata.obs`, cell-indexed.
- `<sample>{suffix}_spatial_coords.csv.gz` — `(x, y)` per cell.

Optionally re-runs the QC filter (`sc.pp.filter_cells(min_counts=10)`) before export — matching the reference notebook, which runs QC on the raw h5ad even though the same filter has already been applied in stage 2. Skip by setting `split_prep.min_counts_cell: null` in a user YAML.

## Stage 4: `rctd_prep`

**Source:** the internal spatial-RCTD preparation reference (first half; the second half — building a `spacexr::Reference` from a scRNA mtx triple — is out of scope for step 1).

Python (`xenium_preprocess.stages.rctd_prep`) shells out to `Rscript src/xenium_preprocess/r/rctd_prep.R`. The R script:

1. `Seurat::ReadMtx(mtx=…, features=…, cells=…, feature.column=1, cell.column=1)` on the SPLIT triple.
2. `SpatialExperiment(assays=list(counts=mtx), colData=DataFrame(meta[colnames(mtx), , drop=FALSE]))`.
3. `spatialCoords(spe) <- as.matrix(coords)` (from `<sample>{suffix}_spatial_coords.csv.gz`).
4. `CreateSeuratObject(counts=counts(spe), assay="Proseg", meta.data=as.data.frame(colData(spe)))`.
5. Attach the spatial coordinates as a `DimReduc` keyed `ST_` (columns renamed to `ST_1`, `ST_2` via `CreateDimReducObject`).
6. `saveRDS(seu, "test_object.rds")`.

The resulting RDS is the RCTD "test object" — downstream `spacexr::create.RCTD(spatial_seurat, reference)` calls take it as the query side of the deconvolution. RCTD's `create.RCTD(..., require_int=TRUE)` accepts integer counts only, which is why `rctd_prep.source_layer` defaults to `maxpost_counts` (the integer proseg output) rather than `expected_counts` (continuous). The layer selection lives in one config knob (`rctd_prep.source_layer`); `pipeline.py` threads it through into `split_prep.layer` so there is a single source of truth for what feeds RCTD.

## Reproducibility

Every run writes a `resolved_config.yaml` at the run root, capturing the exact merged config (default YAML + user YAML + CLI overrides) that the run used. The pipeline's random seeds — `random_state` in stage 2 for KNN + UMAP + Leiden — live in the config, so a re-run with the same seed + same input reproduces the same clustering.
