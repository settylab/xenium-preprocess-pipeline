# Methods

Per-stage attribution to the source scripts + notebooks + summary. Every RCTD / SPLIT parameter this pipeline defaults to has a citation to the exact line in the reference Rmd (see the internal SPLIT/Proseg workflow).

## Stage 1: `rctd_run` — RCTD deconvolution

**Ported from:** the internal SPLIT/Proseg workflow lines 305-363.

```r
# Rmd line 306-307: sanitise reference labels — spacexr rejects '/'
sc$Final_level1_celltype_annotation <- gsub("/", "_",
    as.character(sc$Final_level1_celltype_annotation))

# Rmd line 311: intersect gene panels
common_genes <- intersect(rownames(downsample_proseg), rownames(sc))

# Rmd lines 324-327: build SpatialRNA test object
test.obj <- SpatialRNA(
    coords      = downsample_proseg@reductions$spatial@cell.embeddings %>% as.data.frame(),
    counts      = GetAssayData(downsample_proseg, assay="Proseg", layer="counts")[common_genes,],
    require_int = TRUE
)

# Rmd lines 344-353: create.RCTD
RCTD <- create.RCTD(
    test.obj,
    ref.obj,
    UMI_min           = 10,      # Rmd line 347
    counts_MIN        = 10,      # Rmd line 348
    UMI_min_sigma     = 100,     # Rmd line 349
    max_cores         = 4,       # Rmd line 350
    CELL_MIN_INSTANCE = 20,      # Pipeline default; Rmd line 351 used 25.
                                 # Lowered 2026-07-28 (internal issue review).
    class_df          = NULL     # Rmd line 331: class_df <- NULL
)

# Rmd line 362: run.RCTD
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
```

The pipeline preserves every one of the above parameters as defaults in `config/default.yaml` — see `rctd_run.*`. Overriding a parameter should be motivated by sample characteristics (e.g. larger `max_cores` when a bigger Slurm allocation permits) or by a downstream analysis need (e.g. `doublet_mode="full"` for a mode-comparison study).

### `max_cores` and PSOCK worker library paths

When `max_cores > 1`, `spacexr` runs the per-spot RCTD computation in `process_beads_batch` / `decompose_batch` via `foreach %dopar%` on a `parallel::makeCluster(numCores)` PSOCK cluster it creates itself. PSOCK workers are fresh `Rscript` subprocesses whose `.libPaths()` is derived at startup from the inherited `R_LIBS_USER` environment variable — the parent's in-memory `.libPaths(c(user_lib, .libPaths()))` prepend is **not** visible to them. If `spacexr` (or any package the worker deserialises a closure against) lives on a path that is only in the parent's in-memory `.libPaths()`, workers die with `there is no package called 'spacexr'` during closure deserialisation.

`rctd_run.R` therefore, before invoking `run.RCTD`, mirrors its full `.libPaths()` back into `R_LIBS_USER` so PSOCK workers inherit the same library search path. Pre-registering a cluster with `doParallel::registerDoParallel(cl)` after `parallel::clusterCall(cl, .libPaths, .libPaths())` does **not** work: `spacexr` unconditionally calls `parallel::makeCluster` itself and its `registerDoParallel(cl)` overwrites any pre-registered backend.

## Stage 2: `split_purify` — SPLIT post-process + purify

**Ported from:** the internal SPLIT/Proseg workflow lines 494-523.

```r
# Rmd lines 495-497: run_post_process_RCTD → unpurified variant
rctd_result <- SPLIT::run_post_process_RCTD(rctd_result)
xe <- AddMetaData(xe, rctd_result@results$results_df)
# saveRDS(xe, "unpurified.rds")

# Rmd lines 509-513: SPLIT::purify → purified variant
res_split <- SPLIT::purify(
    counts             = GetAssayData(xe, assay = "Xenium", layer = "counts"),
    rctd               = rctd_result,
    DO_purify_singlets = TRUE          # Rmd line 512
)

# Rmd lines 517-522: rebuild Seurat from purified counts
xe_purified <- CreateSeuratObject(
    counts    = res_split$purified_counts,
    meta.data = res_split$cell_meta,
    assay     = "Xenium"                # rctd-split uses assay_name="Proseg" (default)
)
# saveRDS(xe_purified, "purified.rds")
```

The `assay` swap (Xenium → Proseg) mirrors the xenium-preprocess pipeline's assay convention; the xenium-preprocess `rctd_prep` stage writes the spatial test object with `assay="Proseg"`. All three R scripts read this from `config.<stage>.assay_name`.

**Metadata carried into both variants** (from RCTD's `results_df` + SPLIT's `cell_meta`):
- `first_type` — RCTD's top cell-type call.
- `second_type` — RCTD's second cell-type call.
- `spot_class` — RCTD's spot classification (`singlet`, `doublet_certain`, `doublet_uncertain`, `reject`).
- `purification_status` — SPLIT's per-cell purification status.
- `w1_larger_w2` — SPLIT-emitted boolean.
- `same_class` — SPLIT-emitted boolean.
- `nCount_Proseg` — Seurat-computed total counts for the Proseg assay.
- `x`, `y` — spatial coordinates.

The purified variant carries the same shape but only for cells SPLIT retained after purification.

## Stage 3: `export_mtx` — Seurat RDS → 10X-style bundle

**Convention:** matches xenium-preprocess's `split_prep` and ref-build's `export_mtx` filename shape:
- `<sample>_<variant>_counts.mtx.gz` — MatrixMarket, genes × cells, gzipped.
- `<sample>_<variant>_features.tsv.gz` — one gene per line (single-column form; Seurat `ReadMtx(feature.column=1)`).
- `<sample>_<variant>_barcodes.tsv.gz` — one cell id per line.
- `<sample>_<variant>_metadata.csv` — cell-indexed CSV, all metadata columns from the RDS.
- `<sample>_<variant>_spatial_coords.csv.gz` — `(x, y)` per cell.

The R side reads each of `unpurified.rds` and `purified.rds` and writes both bundles in one Rscript invocation.

## Stage 4: `mtx_to_h5ad` — 10X bundle → AnnData

Pure Python. Reads the mtx (transposes to cells × genes), attaches features → `.var_names`, barcodes → `.obs_names`, `metadata.csv` → `.obs`, `spatial_coords.csv.gz` → `.obsm["spatial"]`, and writes an `.h5ad` per variant with gzip compression by default.

The consumer notebook (Stage E of the summary) uses `sc.pp.filter_cells(min_counts=50)`, `daniel_approach_hvg_pca` for HVG + PCA + UMAP, and `xdp.leiden_clustering_via_knn_graph` at resolutions `0.5` and `0.7`. It colours by `first_type`, `purification_status`, `spot_class`, and `log_count = log10(nCount_Proseg)`. rctd-split preserves every one of those columns in `.obs`.

## Related summary reading

- `the internal reference summary` §Stage D (lines 282-298) — the summary's description of this stage.
- §Stage E (lines 299-330) — how the purified `.h5ad` is consumed downstream.
- §Caveats §4 — the float-counts pattern seen in some purified outputs (informational; RCTD is called with `require_int=TRUE` on the reference side, but the spatial-side `SPLIT::purify` can emit fractional values on some samples).
