# rctd-split

**rctd-split of the Xenium spatial-data preprocessing pipeline: spatial test object (from xenium-preprocess) + scRNA reference (from ref-build) → RCTD run → SPLIT post-process (unpurified) → SPLIT purify (purified) → 10X-style mtx bundles → .h5ad AnnData files.**

At a glance:

- **`rctd_run`** — shells out to `Rscript r/rctd_run.R`. Loads the spatial test object RDS + the `spacexr::Reference` RDS, calls `create.RCTD(...)` with the reference-workflow parameters and `run.RCTD(..., doublet_mode="doublet")`. Saves `rctd_results.rds`.
- **`split_purify`** — shells out to `Rscript r/split_purify.R`. Loads `rctd_results.rds` and the spatial test object; calls `SPLIT::run_post_process_RCTD(...)` (unpurified variant) and `SPLIT::purify(..., DO_purify_singlets=TRUE)` (purified variant). Saves `unpurified.rds` and `purified.rds` — each is a Seurat spatial object carrying RCTD's per-cell metadata (`first_type`, `second_type`, `spot_class`, `purification_status`, `w1_larger_w2`, `same_class`, `nCount_Proseg`, spatial `x`/`y`).
- **`export_mtx`** — shells out to `Rscript r/export_mtx.R`. Reads both `unpurified.rds` and `purified.rds` and writes matched 10X-style bundles per variant: `{counts.mtx.gz, features.tsv.gz, barcodes.tsv.gz, metadata.csv, spatial_coords.csv.gz}`. Same filename convention as xenium-preprocess's `split_prep` and ref-build's `export_mtx`.
- **`mtx_to_h5ad`** — pure Python. Reads each bundle, builds an AnnData with `.var_names` = features, `.obs_names` = barcodes; merges `metadata.csv` into `.obs`; attaches `spatial_coords.csv.gz` into `.obsm['spatial']`. Writes `<sample>_unpurified.h5ad` and `<sample>_purified.h5ad`.
- Every stage is idempotent — sentinel-file resume — and every run snapshots its resolved configuration to disk.

Table of contents: [What it does](#what-it-does) · [Installation](#installation) · [Quickstart](#quickstart) · [Outputs](#outputs) · [Configuration](#configuration) · [HPC (Slurm)](#hpc-usage-slurm) · [Citation](#citation) · [License](#license) · [Acknowledgements](#acknowledgements)

## What it does

Four stages, run in order (`rctd_run → split_purify → export_mtx → mtx_to_h5ad`). Restrict a run to a subset with `--stages`; each stage is guarded by a sentinel file, so re-running with the same `--sample-id` and `--output-root` picks up where the last run left off. Nuke sentinels with `--force-rerun`.

### `rctd_run` — spatial + reference → RCTD deconvolution result

Adapted verbatim from the internal SPLIT/Proseg workflow lines 344-363. Python shells out to `Rscript src/rctd_split/r/rctd_run.R`. The R script:

1. Loads the spatial test object RDS (Seurat with `assay=Proseg`, spatial DimReduc keyed `ST_`, `x`/`y` metadata — the output of xenium-preprocess's `rctd_prep`).
2. Loads the scRNA reference RDS (`spacexr::Reference` — the output of ref-build's `rctd_reference_build`).
3. Intersects the gene panel between spatial and reference.
4. Builds a `spacexr::SpatialRNA` object from the spatial coords + counts.
5. Calls `create.RCTD(test.obj, ref.obj, UMI_min=10, counts_MIN=10, UMI_min_sigma=100, max_cores=<config>, CELL_MIN_INSTANCE=20, class_df=NULL)`. (Pipeline default; Rmd used 25 — lowered on 2026-07-28 per `(internal issue review)`.)
6. Calls `run.RCTD(RCTD, doublet_mode="doublet")`.
7. `saveRDS(RCTD, "rctd_results.rds")`.

**Writes:** `<output_root>/<sample_id>/rctd/rctd_results.rds`.

### `split_purify` — RCTD result → unpurified + purified Seurat objects

Adapted verbatim from the internal SPLIT/Proseg workflow lines 494-523. The R script:

1. Loads `rctd_results.rds` and the spatial test object RDS.
2. Calls `SPLIT::run_post_process_RCTD(rctd_result)`.
3. Attaches `rctd_result@results$results_df` into the spatial Seurat's metadata (`AddMetaData(...)`) — this is the "unpurified" variant (all cells + RCTD calls, before purification).
4. `saveRDS(xe_unpurified, "unpurified.rds")`.
5. Calls `SPLIT::purify(counts=GetAssayData(xe, layer="counts"), rctd=rctd_result, DO_purify_singlets=TRUE)`.
6. Builds a new Seurat from `purified$purified_counts` + `purified$cell_meta`.
7. `saveRDS(xe_purified, "purified.rds")`.

**Writes:** `<output_root>/<sample_id>/split/unpurified.rds`, `<output_root>/<sample_id>/split/purified.rds`.

### `export_mtx` — Seurat RDS → 10X-style bundle

Reads both `unpurified.rds` and `purified.rds`, extracts counts + metadata + spatial coords, and writes matched bundles. Same naming convention as xenium-preprocess's `split_prep` and ref-build's `export_mtx`:

- `<sample>_{unpurified,purified}_counts.mtx.gz` — genes × cells, gzipped MatrixMarket.
- `<sample>_{unpurified,purified}_features.tsv.gz` — one gene name per line (single-column).
- `<sample>_{unpurified,purified}_barcodes.tsv.gz` — one cell id per line.
- `<sample>_{unpurified,purified}_metadata.csv` — cell-indexed CSV.
- `<sample>_{unpurified,purified}_spatial_coords.csv.gz` — `(x, y)` per cell.

**Writes:** `<output_root>/<sample_id>/split/{unpurified,purified}/` with the five files above.

### `mtx_to_h5ad` — 10X bundle → AnnData `.h5ad`

Pure Python. Reads each of the two bundles, builds an AnnData: rows = cells (from `barcodes.tsv.gz`), cols = genes (from `features.tsv.gz`), counts from the transposed mtx. Merges `metadata.csv` into `.obs`; attaches `spatial_coords.csv.gz` as an `(n_obs, 2)` array into `.obsm["spatial"]`. Writes each as an `.h5ad` (compressed).

**Writes:** `<output_root>/<sample_id>/h5ad/<sample>_unpurified.h5ad`, `<output_root>/<sample_id>/h5ad/<sample>_purified.h5ad`.

## Installation

`rctd-split` targets Python 3.10+ (validated on 3.11). The R stages depend on `Seurat`, `Matrix`, `readr`, `spacexr` (RCTD), and `SPLIT`.

`fhR/4.4.1-foss-2023b` carries Seurat + Matrix + readr; `spacexr` and `SPLIT` need a one-time user-local install at `~/.claude/r_libs/4.4.1` — see `docs/install.md`.

```bash
# 1. Install micromamba (skip if you already have it)
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# 2. Create the env
micromamba env create -f environment.yml
micromamba activate rctdSplit

# 3. Editable install
pip install -e /path/to/rctd-split

# 4. Verify
rctd-split --help
rctd-split run --help
```

## Quickstart

```bash
micromamba activate rctdSplit
ml fhR/4.4.1-foss-2023b   # for the R stages; skip if you already have Rscript with Seurat + spacexr + SPLIT

rctd-split run \
    --sample-id       SAMPLE1 \
    --test-object     /data/SAMPLE1/xenium_preprocess/rctd_prep/test_object.rds \
    --reference-rds   /data/SAMPLE1/ref_build/rctd_reference/SAMPLE1_scRNA_ref.rds \
    --output-root     /data/rctd_split_runs
```

Skip a stage set with `--stages`:

```bash
# Regenerate only the h5ad files (assumes mtx bundles exist).
rctd-split run --sample-id SAMPLE1 \
    --test-object /… --reference-rds /… --output-root /… \
    --stages mtx_to_h5ad
```

## Outputs

Under `<output_root>/<sample_id>/`:

```
rctd/rctd_results.rds                              — stage 1
split/unpurified.rds                               — stage 2
split/purified.rds                                 — stage 2
split/unpurified/<sample>_unpurified_counts.mtx.gz — stage 3
split/unpurified/<sample>_unpurified_features.tsv.gz
split/unpurified/<sample>_unpurified_barcodes.tsv.gz
split/unpurified/<sample>_unpurified_metadata.csv
split/unpurified/<sample>_unpurified_spatial_coords.csv.gz
split/purified/<sample>_purified_counts.mtx.gz     — stage 3
split/purified/<sample>_purified_features.tsv.gz
split/purified/<sample>_purified_barcodes.tsv.gz
split/purified/<sample>_purified_metadata.csv
split/purified/<sample>_purified_spatial_coords.csv.gz
h5ad/<sample>_unpurified.h5ad                      — stage 4
h5ad/<sample>_purified.h5ad                        — stage 4
config.yaml                               — snapshot of the resolved config
```

## Configuration

Every knob lives in `config/default.yaml`. Override precedence: CLI flag > user YAML (`--config …`) > default YAML.

The four required-at-runtime keys (`sample_id`, `test_object`, `reference_rds`, `output_root`) are `null` in the default YAML — the pipeline refuses to run until they're supplied. Everything else has a sensible default matched to the reference Rmd + summary.

Full CLI reference: `rctd-split run --help`. Semantics of each knob: `docs/usage.md`.

## HPC usage (Slurm)

```bash
sbatch scripts/submit.slurm.sh SAMPLE1 \
    /data/SAMPLE1/xenium_preprocess/rctd_prep/test_object.rds \
    /data/SAMPLE1/ref_build/rctd_reference/SAMPLE1_scRNA_ref.rds \
    --output-root /data/rctd_split_runs
```

The submit wrapper carries the Slurm-tee/pipefail deadlock fix (see `scripts/submit.slurm.sh` for the comment block). Under sbatch it uses a plain `exec >> "$LOG_FILE" 2>&1`; interactive `bash submit.slurm.sh` still gets the tee for live console echo. Skipping this split re-introduces the silent-freeze deadlock the H&E pipeline hit on 2026-07.

## Citation

If you use `rctd-split` in published work, cite it via the metadata in [`CITATION.cff`](CITATION.cff). Also cite the upstream tools this pipeline glues together:

- **RCTD / spacexr** — Cable *et al.* Nat. Biotechnol. 40, 517–526 (2022).
- **SPLIT** — Sokolov, A. *et al.* bdsc-tds/SPLIT (2024) — a post-processing layer for RCTD that reassigns doublet cell contributions.
- **Seurat** — Hao *et al.* Cell 184, 3573–3587.e29 (2021).
- **anndata / scanpy** — Wolf, F.A. *et al.* Genome Biology 19, 15 (2018); Virshup et al. bioRxiv (2021).

## License

MIT — see [`LICENSE`](LICENSE). See the "Third-party attribution" block in `LICENSE` for the internal SPLIT/Proseg workflow this package is a port of.

## Acknowledgements

Adapted from an internal SPLIT/Proseg workflow; every RCTD/SPLIT parameter is traced to a line in the internal SPLIT/Proseg workflow in the inline comments. See `docs/methods.md` for a per-stage attribution to the source scripts and notebooks.
