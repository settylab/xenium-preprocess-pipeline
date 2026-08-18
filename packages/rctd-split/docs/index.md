# rctd-split

**rctd-split of the Xenium spatial-data preprocessing pipeline.**

Given a spatial test object RDS (from xenium-preprocess's `rctd_prep` stage) and a scRNA reference RDS (from ref-build's `rctd_reference_build` stage), rctd-split:

1. Runs **RCTD** (`create.RCTD` + `run.RCTD(doublet_mode="doublet")`) — saves `rctd_results.rds`.
2. Applies **SPLIT::run_post_process_RCTD** to produce the "unpurified" Seurat spatial object with per-cell RCTD calls attached (`first_type`, `second_type`, `spot_class`, `purification_status`, `w1_larger_w2`, `same_class`, `nCount_Proseg`, `x`, `y`) — saves `unpurified.rds`.
3. Applies **SPLIT::purify** with `DO_purify_singlets=TRUE` to produce the "purified" Seurat object (confident subset with purified counts) — saves `purified.rds`.
4. Writes matching 10X-style mtx bundles for each variant (`counts.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`, `metadata.csv`, `spatial_coords.csv.gz`).
5. Converts each bundle to an `.h5ad` AnnData ready for downstream Python analysis (per Stage E of the summary).

## Table of contents

- [Installation](install.md) — environments (Python + R side).
- [Usage](usage.md) — CLI, config knobs, sample invocations.
- [Methods](methods.md) — per-stage attribution to the user's reference Rmd + summary.

## Quickstart

```bash
micromamba activate rctdSplit
ml fhR/4.4.1-foss-2023b

rctd-split run \
    --sample-id       SAMPLE1 \
    --test-object     /data/SAMPLE1/xenium_preprocess/rctd_prep/test_object.rds \
    --reference-rds   /data/SAMPLE1/ref_build/rctd_reference/SAMPLE1_scRNA_ref.rds \
    --output-root     /data/rctd_split_runs
```

## Outputs

Under `<output_root>/<sample_id>/`:

```
rctd/rctd_results.rds
split/unpurified.rds
split/purified.rds
split/unpurified/<sample>_unpurified_{counts.mtx.gz, features.tsv.gz, barcodes.tsv.gz, metadata.csv, spatial_coords.csv.gz}
split/purified/<sample>_purified_{...same shape...}
h5ad/<sample>_unpurified.h5ad
h5ad/<sample>_purified.h5ad
config.yaml
```

Downstream analysis (Stage E of `ref-build-summary-v3.md`) consumes the purified `.h5ad`.
