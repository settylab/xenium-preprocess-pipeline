# rctd-split

**Function.** Apply RCTD deconvolution to spatial cells, refine singlets via
SPLIT, then write cell-type labels back onto the proseg + xenium-ranger AnnData
files (plus render a per-run QC summary HTML).

Part of the [`xenium-preprocess-pipeline`](../../README.md) chain. Runs after
`ref-build` under a shared `--run-id`.

## Sub-stages (default order)

| Stage | Purpose |
|---|---|
| `rctd_run`             | R shell-out. `spacexr::create.RCTD(...) + run.RCTD(doublet_mode="doublet")`. Writes `<sample>_rctd_results.rds`. |
| `split_purify`         | R shell-out. `SPLIT::run_post_process_RCTD(...)` (unpurified) then `SPLIT::purify(DO_purify_singlets=TRUE)` (purified). Writes two Seurat RDS variants. |
| `export_mtx`           | R shell-out. Each RDS variant → 10X-style bundle (`counts.mtx.gz` + `features.tsv.gz` + `barcodes.tsv.gz` + `metadata.csv` + `spatial_coords.csv.gz`). |
| `mtx_to_h5ad`          | Pure Python. Each bundle → AnnData; writes `<sample>_proseg_purified.h5ad` and an intermediate unpurified adata. |
| `filter_status`        | Populate per-cell filter flags (`passed_rctd`, `passed_split_purify`, `filtered_by_purification`) on the intermediate unpurified adata. |
| `postprocess`          | QC + clipped-log normalization + PCA + KNN + UMAP + Leiden on the purified adata (and optionally the unpurified layer). |
| `writeback_to_raw`     | Merge filter-status + purified obs back onto `<sample>_proseg_raw.h5ad`, so every downstream reader gets celltype + QC without a separate join. |
| `celltype_writeback`   | NN-map celltype from the purified adata onto `<sample>_xenium_ranger.h5ad` by nearest centroid. Adds `.obs['celltype']` + `.obs['celltype_source_distance']`. |
| `qc_report`            | Render `<sample>_summary_report.html` — per-cell counts, spatial plots, celltype pie, provenance (config + package versions + invocation). |

Every stage is sentinel-gated (idempotent resume); pass `--force-rerun` to redo.

## Quickstart

```bash
micromamba activate xenium
ml fhR/4.4.1-foss-2023b     # for rctd_run / split_purify / export_mtx

# Consumes xenium-preprocess's rctd/<sample>_test_object.rds
# and ref-build's rctd/<sample>_reference.rds from the same --run-id.
rctd-split run \
    --sample-id       SAMPLE1 \
    --run-id          demo_v1 \
    --output-root     /data/workflow_runs \
    --max-cores       12
```

Run a subset with `--stages`:

```bash
# Re-render only the QC summary HTML from an existing run
rctd-split run --sample-id SAMPLE1 --run-id demo_v1 --output-root /data/workflow_runs \
    --stages qc_report
```

Mix inputs from different samples via explicit `--test-object` / `--reference-rds`
(both-or-neither):

```bash
rctd-split run --sample-id SAMPLE1 --run-id demo_v1 --output-root /data/workflow_runs \
    --test-object   /data/OTHER1/rctd/OTHER1_test_object.rds \
    --reference-rds /data/OTHER2/rctd/OTHER2_reference.rds
```

## Outputs

Under `<output_root>/<sample_id>/<sample_id>_<run_id>/`:

```
rctd/<sample>_rctd_results.rds                 rctd_run
rctd/<sample>_test_object.rds                  input (from xenium-preprocess)
rctd/<sample>_reference.rds                    input (from ref-build)
intermediate/split/{unpurified,purified}.rds   split_purify (Seurat objects)
intermediate/adata/<sample>_unpurified.h5ad    mtx_to_h5ad — intermediate unpurified adata
spatial_adata/<sample>_proseg_purified.h5ad    mtx_to_h5ad — SPLIT-purified adata
spatial_adata/<sample>_proseg_raw.h5ad         writeback_to_raw augments in place
spatial_adata/<sample>_xenium_ranger.h5ad      celltype_writeback augments in place
summary/<sample>_summary_report.html           qc_report — the per-run summary
config.yaml                                    resolved config (this stage writes the rctd_split: key)
logs/rctd-split.log                            stage log
```

## Configuration

Every knob lives in `config/default.yaml`. Override precedence:
**CLI flag > user YAML (`--config …`) > default YAML**.

Required-at-runtime: `sample_id`, `run_id`, `output_root`. `test_object` +
`reference_rds` auto-discover under the `--run-id` layout when omitted (from
`xenium-preprocess` and `ref-build` respectively); supply both explicitly to mix
inputs from different runs.

Key optional knobs (defaults in parentheses):

| Key | Default | Meaning |
|---|---|---|
| `rctd_run.max_cores`          | `12`             | RCTD parallelism (sbatch cap is 16). |
| `rctd_run.UMI_min`            | `10`             | `create.RCTD(UMI_min=...)`. |
| `rctd_run.CELL_MIN_INSTANCE`  | `20`             | `create.RCTD(CELL_MIN_INSTANCE=...)`. |
| `rctd_run.doublet_mode`       | `doublet`        | `run.RCTD(doublet_mode=...)`. |
| `split_purify.DO_purify_singlets` | `true`       | `SPLIT::purify(DO_purify_singlets=...)`. |
| `postprocess.qc.min_counts`   | `50`             | QC gate on the purified adata. |
| `keep_intermediate`           | `false`          | Retain the intermediate mtx bundles + unpurified adata after `mtx_to_h5ad`. |
| `force_rerun`                 | `false`          | Nuke every sentinel and redo from scratch. |

Full CLI reference: `rctd-split run --help`.

## Installation

```bash
micromamba create -n xenium -f environment.yml   # once
micromamba activate xenium
uv pip install -e .
```

For the R stages, load an R module with Seurat + spacexr + SPLIT + Matrix — on
Fred Hutch:

```bash
ml fhR/4.4.1-foss-2023b
```

## Citation

`CITATION.cff`. Upstream tools to cite alongside:

- **RCTD / spacexr** — Cable *et al.* Nat. Biotechnol. 40, 517–526 (2022).
- **SPLIT** — Sokolov, A. *et al.* bdsc-tds/SPLIT (2024) — RCTD post-processing.
- **Seurat** — Hao *et al.* Cell 184, 3573–3587.e29 (2021).
- **anndata / scanpy** — Wolf, F.A. *et al.* Genome Biology 19, 15 (2018); Virshup *et al.* bioRxiv (2021).

## License

MIT — see [`LICENSE`](LICENSE).
