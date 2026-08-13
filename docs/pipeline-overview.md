# Pipeline overview

`xenium-preprocess-pipeline` chains three packages via a shell driver:

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

## Step 1 — `xenium-preprocess`

Turns a proseg cell-segmentation run + a Xenium Ranger bundle into
per-sample h5ads and an RCTD test object.

Stages (sentinel-gated):

1. `preprocess` — proseg output → `<sample>_proseg_raw.h5ad`.
2. `xenium_ranger` — Xenium Ranger output → `<sample>_xenium_ranger.h5ad`.
3. `provenance` — copies proseg run scripts + configs into
   `spatial_adata/provenance/`.
4. `mtx` — writes matrix/barcodes/features for downstream SPLIT.
5. `rctd_prep` — builds the RCTD test object (spatial query).

Outputs: `spatial_adata/`, `rctd/<sample>_test_object.rds`.

Full CLI: `xenium-preprocess --help`; docs under
`packages/xenium-preprocess/docs/`.

## Step 3 — `ref-build`

Builds a celltype-marker-driven RCTD reference from a Flex scRNA
h5ad + optional per-donor h5ads.

Key inputs:

- `--flex-h5ad` — 10x Flex scRNA reference.
- `--celltype-marker-json` — marker gene JSON.
- `--donor-h5ad` (repeatable) — per-donor h5ads mixed into the reference.
- `--fallback-donor-h5ad` — fallback for missing donors.

Outputs: `rctd/<sample>_reference.rds`.

Full CLI: `ref-build --help`; docs under `packages/ref-build/docs/`.

## Step 4 — `rctd-split`

Runs RCTD deconvolution + SPLIT typing against the test object (step 1)
+ reference (step 3).

Key knobs:

- `--max-cores` — parallelism cap (benefits from 16 CPUs).

Outputs: `rctd/<sample>_rctd_split.rds`, plus a typed `spatial_adata/`
h5ad written in-place.

Full CLI: `rctd-split --help`; docs under `packages/rctd-split/docs/`.

## Workflow driver — `scripts/submit_workflow.sh`

Chains the three steps via Slurm `--dependency=afterok:`. Steps 1
and 3 are independent (both fan out from the same driver call);
step 4 depends on both.

`scripts/submit_step{1,3,4}.sbatch` are the per-step Slurm templates
the driver submits. Each writes a per-step log to
`<run-dir>/logs/slurm-<jobid>-<step>.out`; the driver itself writes
`workflow-submit.log` recording the dependency chain.

Resume:

- `--start-step 3` or `--start-step 4` — reuse an existing run folder;
  each step overwrites only the files it writes.
- `--force` — `rm -rf` the run folder and restart from step 1
  (destructive; incompatible with `--start-step 3`/`4`).

## Configuration

Each package ships `config/default.yaml`. Overrides in precedence order
(highest first): CLI flag, `--config user.yaml`, package default.
The full effective config for a run is recorded at
`<run-dir>/resolved_config.yaml`.
