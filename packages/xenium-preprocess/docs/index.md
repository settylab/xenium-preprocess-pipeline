# xenium-preprocess

`xenium-preprocess` takes a **proseg output directory** (the segmented-cell output of the proseg cell-segmentation tool on a Xenium in-situ transcriptomics run) and produces four artifacts, in order:

- a **raw AnnData** with `X` = proseg counts, `.obs` = proseg cell metadata, `.obsm['spatial']` = centroids,
- a **preprocessed "unpurified" AnnData** (`adata_unpurified.h5ad`) with clipped-log normalization, PCA on `X_clipped`, KNN + UMAP, Leiden clustering, and optional inferred celltype labels — plus a **UMAP plot** for each colouring,
- a **10x-style mtx/features/barcodes triple** derived from the raw proseg counts (for SPLIT downstream analyses),
- an **RCTD "test object" RDS** built by the R script (`rctd_prep.R`) from the mtx triple — the query side of `spacexr::create.RCTD`.

Every stage is resumable via a sentinel-file existence check; every run snapshots its resolved config to disk for traceability.

This is **step 1** of the Xenium spatial-data preprocessing pipeline. Later steps (celltype-curation, deconvolution, differential expression, spatial pattern discovery) will land in sibling packages.

## Pipeline shape

```
   proseg output directory
                │
                ▼
   ┌─────────────────────────┐
   │ 1. proseg_to_anndata    │  count matrix + cell metadata → raw.h5ad
   └─────────────────────────┘
                │
                ▼
   ┌─────────────────────────┐
   │ 2. preprocess           │  QC + normalize + Leiden + celltype infer + UMAP plots
   └─────────────────────────┘
                │
                ▼
   ┌─────────────────────────┐
   │ 3. split_prep           │  mtx + features + barcodes + metadata + spatial coords
   └─────────────────────────┘
                │
                ▼
   ┌─────────────────────────┐
   │ 4. rctd_prep            │  Rscript rctd_prep.R → test_object.rds
   └─────────────────────────┘
```

## Reading order

1. [`install.md`](install.md) — set up the `xeniumPreprocess` conda env, load `fhR/…`, `pip install -e .`.
2. [`usage.md`](usage.md) — CLI reference, config surface, worked examples.
3. [`methods.md`](methods.md) — one-paragraph algorithm summary per stage, in methods-section shape, honest citation to the reference scripts and notebooks.
