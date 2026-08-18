# ref-build

`ref-build` takes a **primary sample's preprocessed scRNA AnnData** (from step 2, `flex-preprocess`) and a **pool of same-primary-tumor-type donor AnnData files**, applies the per-celltype adaptive migration rules, and produces two artifacts:

- a **10X-style mtx bundle** — `{counts.mtx.gz, features.tsv.gz, barcodes.tsv.gz, metadata.csv}` — the file convention xenium-preprocess's `split_prep` writes for the spatial side, so downstream R tooling sees a consistent layout across pipelines,
- a **spacexr `Reference` `.rds`** — the RCTD-ready reference object that `spacexr::create.RCTD(spatial, reference)` takes as the reference side of the deconvolution.

Along the way it emits `census.csv`, an operator-visible audit trail that records the per-celltype migration decision (`primary_only` / `balanced` / `borrowed` / `missing_no_donor`) and where each celltype's cells came from.

This is **ref-build** of the Xenium spatial-data preprocessing pipeline. xenium-preprocess (`xenium-preprocess`) prepares the spatial side; step 2 (`flex-preprocess`) prepares each sample's `_preprocessed_scRNA.h5ad`; ref-build (this package) merges the primary + donor pool into the reference; downstream RCTD/SPLIT drives the deconvolution.

## Pipeline shape

```
   primary + donor preprocessed scRNA h5ads
                        │
                        ▼
    ┌────────────────────────────────────┐
    │ 1. load_primary_and_donors         │  outer-join concat + sample_ID
    └────────────────────────────────────┘
                        │
                        ▼
    ┌────────────────────────────────────┐
    │ 2. census                          │  the rules → census.csv audit trail
    └────────────────────────────────────┘
                        │
                        ▼
    ┌────────────────────────────────────┐
    │ 3. assemble                        │  apply per-celltype policy → reference h5ad
    └────────────────────────────────────┘
                        │
                        ▼
    ┌────────────────────────────────────┐
    │ 4. export_mtx                      │  10X-style mtx + features + barcodes + metadata
    └────────────────────────────────────┘
                        │
                        ▼
    ┌────────────────────────────────────┐
    │ 5. rctd_reference_build            │  Rscript → spacexr::Reference → .rds
    └────────────────────────────────────┘
```

## Reading order

1. [`install.md`](install.md) — set up the `refBuild` conda env, load `fhR/…`, install `spacexr` into `~/.claude/r_libs/4.4.1`, `pip install -e .`.
2. [`usage.md`](usage.md) — CLI reference, config surface, worked examples for each tumor-type template.
3. [`methods.md`](methods.md) — one-paragraph algorithm summary per stage, in methods-section shape, honest citation to the reference notebook family and the summary document.
