# Usage

## Basic invocation

```bash
rctd-split run \
    --sample-id       SAMPLE1 \
    --test-object     /path/to/step1/rctd_prep/test_object.rds \
    --reference-rds   /data/SAMPLE1/ref_build/rctd_reference/SAMPLE1_scRNA_ref.rds \
    --output-root     /path/to/output/root
```

The three required inputs are:
- `--sample-id` — used as the per-sample sub-directory under `--output-root`.
- `--test-object` — the spatial test object RDS, output of step 1's `rctd_prep` stage. A Seurat object with `assay="Proseg"`, a spatial `DimReduc` keyed `ST_`, and `x`/`y` metadata columns.
- `--reference-rds` — the scRNA reference RDS, output of step 3's `rctd_reference_build` stage. A `spacexr::Reference` object.

## Stages

Restrict a run to a subset of stages with `--stages`:

```bash
# Regenerate only the h5ad files (assumes mtx bundles exist).
rctd-split run --sample-id SAMPLE1 --test-object /… --reference-rds /… --output-root /… \
    --stages mtx_to_h5ad
```

Valid choices: `rctd_run`, `split_purify`, `export_mtx`, `mtx_to_h5ad`.

Sentinel-file resume: each stage's anchor output (e.g. `rctd_results.rds`, `purified.rds`, `<sample>_purified_counts.mtx.gz`, `<sample>_purified.h5ad`) is checked before the stage runs; if present, the stage is skipped. Nuke sentinels with `--force-rerun`.

## Configuration

Every knob lives in `config/default.yaml`. Override precedence: CLI flag > user YAML (`--config …`) > default YAML.

### The knobs the user is most likely to change

- `--max-cores` — RCTD parallelism (default `4`, matches Rmd line 350). Increase for larger samples if your Slurm allocation permits.
- `--umi-min`, `--counts-min`, `--umi-min-sigma`, `--cell-min-instance` — the `create.RCTD` gating parameters. Defaults are `10 / 10 / 100 / 20`. The first three match the Rmd (lines 347-350); `--cell-min-instance` was lowered from the Rmd's `25` (line 351) to `20` on 2026-07-28 (`(internal issue review)`) so detailed-annotation references with smaller per-type counts retain rare celltypes at RCTD time.
- `--doublet-mode` — `"doublet"` (Rmd line 362) is the only mode validated by this pipeline; `"full"` and `"multi"` would need a config-file override.
- `--do-purify-singlets` — `SPLIT::purify(DO_purify_singlets=…)`. Defaults `true` (Rmd line 512).
- `--assay-name` — the Seurat assay carrying spatial counts on the test object. Default `Proseg` (matches step 1's `rctd_prep`).

### Full config surface

Everything in `config/default.yaml`:

```yaml
sample_id: null
test_object: null
reference_rds: null
output_root: null
force_rerun: false

rctd_run:
  rscript_bin: Rscript
  assay_name: Proseg
  spatial_reduction: spatial
  UMI_min: 10
  counts_MIN: 10
  UMI_min_sigma: 100
  max_cores: 4
  CELL_MIN_INSTANCE: 20    # Pipeline default; Rmd used 25.
  doublet_mode: doublet
  label_slash_replacement: "_"
  output_filename: rctd_results.rds

split_purify:
  rscript_bin: Rscript
  assay_name: Proseg
  DO_purify_singlets: true
  unpurified_filename: unpurified.rds
  purified_filename: purified.rds

export_mtx:
  rscript_bin: Rscript
  gzip_outputs: true
  assay_name: Proseg

mtx_to_h5ad:
  h5ad_compression: gzip
```

Override an individual knob with a user YAML:

```yaml
# my_overrides.yaml
rctd_run:
  max_cores: 16
  UMI_min: 20
```

```bash
rctd-split run --config my_overrides.yaml --sample-id SAMPLE1 --test-object … --reference-rds … --output-root …
```

## HPC (Slurm)

The submit wrapper is at `scripts/submit.slurm.sh`:

```bash
sbatch scripts/submit.slurm.sh SAMPLE1 \
    /data/SAMPLE1/xenium_preprocess/rctd_prep/test_object.rds \
    /data/SAMPLE1/ref_build/rctd_reference/SAMPLE1_scRNA_ref.rds \
    --output-root /data/rctd_split_runs
```

Defaults: `--cpus-per-task=8`, `--mem=128G`, `--time=1-00:00:00`, `--partition=YOUR_PARTITION`. Override in place if the sample needs more.

The wrapper:
- Redirects logs to `<output_root>/<sample_id>/logs/<job_tag>.log` after validating inputs (splits the sbatch-vs-interactive tee/pipefail deadlock the same way as step 1).
- Activates `rctdSplit` micromamba env (override with `ENV_NAME=your_env`).
- `ml fhR/4.4.1-foss-2023b` before invoking the pipeline (override with `R_MODULE=…`).

## What the outputs look like

Under `<output_root>/<sample_id>/`:

```
rctd/rctd_results.rds                                — stage 1
split/unpurified.rds                                 — stage 2 (unpurified Seurat)
split/purified.rds                                   — stage 2 (purified Seurat)
split/unpurified/<sample>_unpurified_counts.mtx.gz   — stage 3
split/unpurified/<sample>_unpurified_features.tsv.gz
split/unpurified/<sample>_unpurified_barcodes.tsv.gz
split/unpurified/<sample>_unpurified_metadata.csv
split/unpurified/<sample>_unpurified_spatial_coords.csv.gz
split/purified/<sample>_purified_counts.mtx.gz       — stage 3
split/purified/<sample>_purified_features.tsv.gz
split/purified/<sample>_purified_barcodes.tsv.gz
split/purified/<sample>_purified_metadata.csv
split/purified/<sample>_purified_spatial_coords.csv.gz
h5ad/<sample>_unpurified.h5ad                        — stage 4
h5ad/<sample>_purified.h5ad                          — stage 4
resolved_config.yaml                                 — provenance snapshot
logs/<job_tag>.log                                   — from Slurm submit wrapper
```

## Downstream use of the h5ad

Stage E of `ref-build-summary-v3.md` (`the noGeneFilter reference notebooks`) reads the purified `.h5ad`, runs `sc.pp.filter_cells(min_counts=50)` re-QC, then `daniel_approach_hvg_pca` (clipped log-norm + HVG + PCA) and `xdp.leiden_clustering_via_knn_graph`. The `.obs` columns rctd-split preserves — `first_type`, `purification_status`, `w1_larger_w2`, `same_class`, `nCount_Proseg`, `x`, `y` — are exactly what those notebooks colour by.
