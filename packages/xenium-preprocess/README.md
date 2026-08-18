# xenium-preprocess

**Function.** Convert proseg + Xenium Ranger outputs into AnnData objects with QC
filter flags, plus an RCTD-ready test-object RDS for downstream deconvolution.

Part of the [`xenium-preprocess-pipeline`](../../README.md) chain. Runs before
`ref-build` and `rctd-split` under a shared `--run-id`.

## Sub-stages (default order)

| Stage | Purpose |
|---|---|
| `proseg_to_anndata`      | proseg output → `<sample>_proseg_raw.h5ad` (expected + maxpost count layers, centroids in `.obsm['spatial']`). |
| `qc_filter`              | Annotate `<sample>_proseg_raw.h5ad` with `.obs['qc_filtered']`. Additive — no cells removed. |
| `xenium_ranger_to_anndata` | Xenium Ranger bundle → `<sample>_xenium_ranger.h5ad` (counts, spatial, standard xenium QC obs). |
| `enrich_xenium_id`       | Add `proseg_cell_id_nn` / `proseg_id_nn_distance` / `proseg_id_nn_note` to the xenium-ranger h5ad by nearest-centroid match against the proseg h5ad. |
| `rctd_prep`              | R shell-out that builds the RCTD spatial test-object RDS from the maxpost counts (required-integer input for RCTD). |

Legacy stages `preprocess`, `split_prep` (10x-style bundle) are kept in
`VALID_STAGES` for reversibility but not in `DEFAULT_STAGES`; opt in with
`--stages ... preprocess split_prep ...`.

Every stage is sentinel-gated (idempotent resume); pass `--force-rerun` to redo.

## Quickstart

```bash
micromamba activate xenium
ml fhR/4.4.1-foss-2023b     # for the R stage

xenium-preprocess run \
    --sample-id       SAMPLE1 \
    --run-id          demo_v1 \
    --output-root     /data/workflow_runs \
    --proseg-dir      /data/SAMPLE1/proseg \
    --xenium-ranger-dir /data/SAMPLE1/xenium_ranger
```

Run a subset with `--stages`:

```bash
xenium-preprocess run --sample-id SAMPLE1 --run-id demo_v1 \
    --output-root /data/workflow_runs \
    --stages xenium_ranger_to_anndata enrich_xenium_id
```

## Outputs

Under `<output_root>/<sample_id>/<sample_id>_<run_id>/`:

```
spatial_adata/<sample>_proseg_raw.h5ad         proseg counts + qc_filter flag + enrich_xenium_id join
spatial_adata/<sample>_xenium_ranger.h5ad      xenium-ranger counts + proseg cell-id nearest-neighbour join
rctd/<sample>_test_object.rds                  RCTD-ready spatial test object (Seurat + spacexr contract)
config.yaml                                    resolved config (this stage writes the xenium_preprocess: key)
logs/xenium-preprocess.log                     stage log
```

## Configuration

Every knob lives in `config/default.yaml`. Override precedence:
**CLI flag > user YAML (`--config …`) > default YAML**.

Required-at-runtime: `sample_id`, `run_id`, `output_root`, `proseg_dir`.
Full CLI reference: `xenium-preprocess run --help`.

## Installation

```bash
micromamba create -n xenium -f environment.yml   # once
micromamba activate xenium
uv pip install -e .
```

For the R stage (`rctd_prep`), also load an R module carrying Seurat,
SpatialExperiment, and spacexr — on Fred Hutch:

```bash
ml fhR/4.4.1-foss-2023b
```

## Citation

`CITATION.cff`. Upstream tools to cite alongside:

- **proseg** — Jones et al., "Cell segmentation of Xenium in-situ transcriptomics data with proseg."
- **scanpy** — Wolf, F.A. *et al.* Genome Biology 19, 15 (2018).
- **RCTD / spacexr** — Cable *et al.* Nat. Biotechnol. 40, 517–526 (2022).

## License

MIT — see [`LICENSE`](LICENSE).
