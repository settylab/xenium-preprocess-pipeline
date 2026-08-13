# xenium-preprocess-pipeline

End-to-end spatial-genomics workflow for 10x Xenium slides: proseg
cell-segmentation preprocessing → celltype-marker reference building → RCTD
deconvolution + SPLIT processing.

## Motivation

This pipeline standardises the workflow that start from xenium cells to cleaned-up proseg segmented cells via a scRNA reference-based purification method called SPLIT. `xenium-preprocess-pipeline` bundles the chain — proseg → xenium-ranger → RCTD
reference → RCTD + SPLIT processing — into three publishable packages plus one
shell driver, so a new sample runs end-to-end from a single command, with
every intermediate product traceable to the exact config that built it.

## Quick start

```bash
git clone https://github.com/settylab/xenium-preprocess-pipeline
cd xenium-preprocess-pipeline

micromamba create -n xenium -f environments/xenium.yml
micromamba activate xenium

uv pip install -e packages/xenium-preprocess
uv pip install -e packages/ref-build
uv pip install -e packages/rctd-split

./scripts/submit_workflow.sh \
    --sample-id             SAMPLE1 \
    --output-root           /data/workflow_runs \
    --run-id                demo_v1 \
    --flex-h5ad             /data/SAMPLE1/scRNA/SAMPLE1_flex.h5ad \
    --celltype-marker-json  /data/markers/markers.json \
    --proseg-dir            /data/SAMPLE1/proseg \
    --xenium-ranger-dir     /data/SAMPLE1/xenium_ranger
```

Outputs land under `/data/workflow_runs/SAMPLE1/SAMPLE1_demo_v1/`.
See [Outputs](#outputs) for the tree.

## Pipeline overview

Three publishable packages + one shell driver, chained by Slurm `afterok:`
dependencies:

```
    Xenium raw            step 1                        spatial_adata/
    ┌──────────┐        ┌──────────────────┐          ┌──────────────────┐
    │ proseg   │───────▶│ xenium-preprocess│─────────▶│ *_proseg_raw.h5ad│
    │ xranger  │        │ (5 stages)       │          │ *_xenium_ranger  │
    └──────────┘        └──────────────────┘          │   .h5ad          │
                                 │                    └──────────────────┘
                                 │                             │
                                 ▼                             │
                          rctd/*_test_object.rds               │
                                                               │
    Flex scRNA          step 3                                 │
    ┌──────────┐        ┌──────────────────┐                   │
    │ flex.h5ad│───────▶│ ref-build        │                   │
    │ donors   │        │ (celltype-marker │                   │
    │ markers  │        │  driven ref)     │                   │
    └──────────┘        └──────────────────┘                   │
                                 │                             │
                                 ▼                             │
                          rctd/*_reference.rds                 │
                                                               │
                        step 4                                 │
                        ┌──────────────────┐                   │
                        │ rctd-split       │◀──────────────────┘
                        │ (RCTD + SPLIT    │
                        │  typing)         │
                        └──────────────────┘
                                 │
                                 ▼
                          rctd/*_rctd_split.rds
                          spatial_adata/ (typed)
```

Each package is independently installable, testable, and reusable outside
this pipeline. The shell driver (`scripts/submit_workflow.sh`) is a thin
Slurm layer on top.

## Installation

### System requirements

- **OS**: Linux (developed and tested on a Slurm cluster).
- **Scheduler**: Slurm (`sbatch`, `--dependency=afterok:` support).
- **RAM**: ~64 GB per step for typical samples; step 4 benefits from 16 CPUs.
- **Disk**: ~50 GB per sample per run for intermediate + final outputs.
- **Micromamba** (or `conda`), **`uv`**, and **R 4.4+** with `spacexr`,
  `SPLIT`, and `Seurat`.

### Python environment

```bash
# Once per user / per machine
micromamba create -n xenium -f environments/xenium.yml
micromamba activate xenium

# Editable installs of the three packages
uv pip install -e packages/xenium-preprocess
uv pip install -e packages/ref-build
uv pip install -e packages/rctd-split

# Verify
xenium-preprocess --help
ref-build         --help
rctd-split        --help
```

### R dependencies (for `ref-build` and `rctd-split`)

`ref-build` and `rctd-split` shell out to `Rscript`. Both stages need
Seurat, Matrix, spacexr (RCTD), and SPLIT.

On the cluster's `fhR/4.4.1-foss-2023b` module provides
Seurat + Matrix + SpatialExperiment; `spacexr` and `SPLIT` install into
a user library — see [`docs/installation.md`](docs/installation.md) for
the recipe.

## Usage

### Run the full pipeline (steps 1 → 3 → 4)

```bash
./scripts/submit_workflow.sh \
    --sample-id             SAMPLE1 \
    --output-root           /data/workflow_runs \
    --run-id                demo_v1 \
    --flex-h5ad             /data/SAMPLE1/scRNA/SAMPLE1_flex.h5ad \
    --celltype-marker-json  /data/markers/markers.json \
    --proseg-dir            /data/SAMPLE1/proseg \
    --xenium-ranger-dir     /data/SAMPLE1/xenium_ranger \
    --donor-h5ad            /data/SAMPLE2/scRNA/SAMPLE2_flex.h5ad \
    --donor-h5ad            /data/SAMPLE3/scRNA/SAMPLE3_flex.h5ad
```

Submits three chained Slurm jobs via `--dependency=afterok:` and prints
their job ids + dependency chain. All jobs share the `--run-id` so
outputs colocate under one folder.

### Run one step directly

Each package installs its own CLI:

```bash
xenium-preprocess run \
    --sample-id   SAMPLE1 --run-id demo_v1 \
    --output-root /data/workflow_runs \
    --proseg-dir  /data/SAMPLE1/proseg --xenium-ranger-dir /data/SAMPLE1/xenium_ranger

ref-build run \
    --sample-id            SAMPLE1 --run-id demo_v1 \
    --output-root          /data/workflow_runs \
    --flex-h5ad            /data/SAMPLE1/scRNA/SAMPLE1_flex.h5ad \
    --celltype-marker-json /data/markers/markers.json \
    --donor-h5ad           /data/SAMPLE2/scRNA/SAMPLE2_flex.h5ad

rctd-split run \
    --sample-id   SAMPLE1 --run-id demo_v1 \
    --output-root /data/workflow_runs \
    --max-cores   12
```

### Run a subset of sub-stages

Each `run` subcommand takes a `--stages` flag; pass a subset to run only
part of a step:

```bash
xenium-preprocess run \
    --sample-id SAMPLE1 --run-id demo_v1 \
    --output-root /data/workflow_runs \
    --stages preprocess mtx

ref-build run \
    --sample-id SAMPLE1 --run-id demo_v1 \
    --output-root /data/workflow_runs \
    --flex-h5ad /data/SAMPLE1/scRNA/SAMPLE1_flex.h5ad --celltype-marker-json /data/markers/markers.json \
    --stages assemble_reference validate_reference
```

Stages are sentinel-gated — reruns skip already-done work; pass
`--force-rerun` to redo.

### Resume a run from a specific step

```bash
# Replace the YOUR_PARTITION with your own cluster's partition. At here, campus-new is the example.

grep -rn YOUR_PARTITION scripts/                       # inventory hits
sed -i 's/YOUR_PARTITION/campus-new/g' scripts/*.sbatch scripts/*.sh
grep -rn YOUR_PARTITION scripts/                       # should print nothing now

# Step 1 already ran; resume the chain at step 3
./scripts/submit_workflow.sh \
    --sample-id             SAMPLE1 \
    --output-root           /data/workflow_runs \
    --run-id                demo_v1 \
    --flex-h5ad             /data/SAMPLE1/scRNA/SAMPLE1_flex.h5ad \
    --celltype-marker-json  /data/markers/markers.json \
    --start-step            3
# --start-step > 1 implies --reuse-run-dir
```

Resume policy:

- `--reuse-run-dir` (implied by `--start-step 3` or `4`) — proceed
  against an existing run folder; each step overwrites only the files it
  writes, everything else is preserved. This is the everyday "resume the
  chain" flag.
- `--force` — `rm -rf` the run folder, then run from step 1. Destructive;
  for a from-scratch re-run under an already-used `--run-id`. Rejected
  together with `--start-step 3` / `4`.

## Configuration

Each package ships a `config/default.yaml` with sensible defaults.
Overrides in precedence order (highest first):

1. CLI flag (e.g. `--x-source expected_counts`).
2. `--config <user.yaml>` — user YAML overriding default keys.
3. `config/default.yaml` — package-shipped defaults.

The full effective config for a run is recorded at
`<run-dir>/resolved_config.yaml`.

Per-step docs live under each package's `docs/` directory.

## Outputs

All three steps write into a single run folder:

```
<output-root>/<sample>/<sample>_<run-id>/
├── spatial_adata/                        # step 1 (step 4 augments in place)
│   ├── <sample>_proseg_raw.h5ad          # proseg cell×gene
│   ├── <sample>_xenium_ranger.h5ad       # xenium-ranger cell×gene
│   └── provenance/                       # verbatim copies of proseg-run scripts + configs
├── rctd/
│   ├── <sample>_test_object.rds          # step 1 — RCTD test object (spatial query)
│   ├── <sample>_reference.rds            # step 3 — RCTD reference (celltype pool)
│   └── <sample>_rctd_split.rds           # step 4 — RCTD + SPLIT typing result
├── resolved_config.yaml                  # merged effective config across the three steps
└── logs/
    ├── slurm-<jobid>-xenium-preprocess.out
    ├── slurm-<jobid>-ref-build.out
    ├── slurm-<jobid>-rctd-split.out
    ├── step1.log  step3.log  step4.log   # per-step app logs
    └── workflow-submit.log               # authoritative record of jobids + dep chain
```

## Testing

```bash
# Per-package pytest suites
pytest packages/xenium-preprocess/tests
pytest packages/ref-build/tests
pytest packages/rctd-split/tests

# Shell driver smoke test
pytest scripts/tests
```

## Contributing

Contributions welcome:

- Open an issue on this repo before starting non-trivial work.
- Branch from `main`; open a PR into `main`.
- Keep per-package changes in per-package PRs when possible.
- `pytest` must pass locally before pushing.

## Citation

If you use `xenium-preprocess-pipeline` in a publication, please cite:

> _TBD — placeholder pending the associated manuscript._

Per-package citations live in each `packages/*/CITATION.cff`.

## License

MIT — see [`LICENSE`](./LICENSE).
