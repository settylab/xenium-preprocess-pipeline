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
    Xenium raw            xenium-preprocess             spatial_adata/
    ┌──────────┐        ┌──────────────────┐          ┌──────────────────┐
    │ proseg   │───────▶│ xenium-preprocess│─────────▶│ *_proseg_raw.h5ad│
    │ xranger  │        │ (5 stages)       │          │ *_xenium_ranger  │
    └──────────┘        └──────────────────┘          │   .h5ad          │
                                 │                    └──────────────────┘
                                 │                             │
                                 ▼                             │
                          rctd/*_test_object.rds               │
                                                               │
    Flex scRNA          ref-build                              │
    ┌──────────┐        ┌──────────────────┐                   │
    │ flex.h5ad│───────▶│ ref-build        │                   │
    │ donors   │        │ (celltype-marker │                   │
    │ markers  │        │  driven ref)     │                   │
    └──────────┘        └──────────────────┘                   │
                                 │                             │
                                 ▼                             │
                          rctd/*_reference.rds                 │
                                                               │
                        rctd-split                             │
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
- **RAM**: ~64 GB per step for typical samples; rctd-split benefits from 16 CPUs.
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

`scripts/submit_workflow.sh` is the primary entry point — it chains
xenium-preprocess → ref-build → rctd-split via Slurm
`--dependency=afterok:` and exposes step + sub-stage selection flags.
For running a single package's CLI in isolation (no Slurm; quick
sanity checks on a laptop; debugging one step in isolation), see
[Advanced: single-package usage without the driver](#advanced-single-package-usage-without-the-driver)
below.

### Via the driver (recommended)

#### Run the full pipeline (xenium-preprocess → ref-build → rctd-split)

```bash
# Replace the YOUR_PARTITION with your own cluster's partition. At here, campus-new is the example.

grep -rn YOUR_PARTITION scripts/                       # inventory hits
sed -i 's/YOUR_PARTITION/campus-new/g' scripts/*.sbatch scripts/*.sh
grep -rn YOUR_PARTITION scripts/                       # should print nothing now

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

#### Resume a run from a specific step

```bash

# xenium-preprocess already ran; resume the chain at ref-build
./scripts/submit_workflow.sh \
    --sample-id             SAMPLE1 \
    --output-root           /data/workflow_runs \
    --run-id                demo_v1 \
    --flex-h5ad             /data/SAMPLE1/scRNA/SAMPLE1_flex.h5ad \
    --celltype-marker-json  /data/markers/markers.json \
    --start-step            ref-build
# --start-step past xenium-preprocess implies --reuse-run-dir
```

Resume behaviour depends on which flag you pass. Three levels of
"force", coarsest first:

- `--force` — **DESTRUCTIVE**: `rm -rf` the entire
  `<output-root>/<sample>/<sample>_<run-id>/` folder, then re-run
  the whole pipeline from `xenium-preprocess`. Use ONLY for a
  clean-slate re-run under an already-used `--run-id`. Rejected
  together with `--start-step ref-build` / `--start-step
  rctd-split` (which would wipe the very outputs those flags need
  to resume from), and with `--reuse-run-dir`.
- `--reuse-run-dir` (implied by `--start-step ref-build` or
  `--start-step rctd-split`) — proceed against the existing run
  folder; each step overwrites only the files it writes,
  everything else preserved. This is the everyday "resume the
  chain" flag.
- **Per-step force-rerun** — narrower than `--force`; wipes only
  ONE step's sentinels + outputs and re-runs that step, leaving
  other steps' outputs intact:
  - `--xenium-preprocess-force-rerun` — re-run xenium-preprocess only.
  - `--ref-build-force-rerun` — re-run ref-build only.
  - `--rctd-split-force-rerun` — re-run rctd-split only.

  Combines with `--start-step` and `--<pkg>-stages` (see next
  section) for even narrower re-runs — e.g. re-render only
  rctd-split's `qc_report` sub-stage against an existing run:

  ```bash
  ./scripts/submit_workflow.sh \
      --sample-id             SAMPLE1 \
      --output-root           /data/workflow_runs \
      --run-id                demo_v1 \
      --flex-h5ad             /data/SAMPLE1/scRNA/SAMPLE1_flex.h5ad \
      --celltype-marker-json  /data/markers/markers.json \
      --start-step            rctd-split \
      --rctd-split-stages     qc_report \
      --rctd-split-force-rerun
  ```

#### Run a subset of sub-stages

The driver forwards a per-step sub-stage subset via three semantic
flags — `--xenium-preprocess-stages`, `--ref-build-stages`,
`--rctd-split-stages` — each taking a comma-separated stage list.
Composes with `--start-step` and the per-step force-rerun flags:

```bash
# ref-build was interrupted mid-pipeline; resume from `census`
# (assumes load_primary_and_donors already wrote intermediate/loaded/concat.h5ad
# on the previous run under the same --run-id).
./scripts/submit_workflow.sh \
    --sample-id            SAMPLE1 \
    --output-root          /data/workflow_runs \
    --run-id               demo_v1 \
    --flex-h5ad            /data/SAMPLE1/scRNA/SAMPLE1_flex.h5ad \
    --celltype-marker-json /data/markers/markers.json \
    --start-step           ref-build \
    --ref-build-stages     census,assemble,export_mtx,rctd_reference_build
```

Omit `--<pkg>-stages` to fall back to that package's
`DEFAULT_STAGES` (the full list). Sub-stages are sentinel-gated —
re-runs skip work whose sentinel already exists — so if you want
to force a specific sub-stage to redo, combine `--<pkg>-stages
<name>` with the matching `--<pkg>-force-rerun` (see above).

### Advanced: single-package usage without the driver

Use these paths when you don't have Slurm, or when you want to run ONE
package in isolation (debugging one step, quick sanity check on a
laptop). Each package installs its own CLI via
`pip install -e packages/<name>` and accepts the same run scope
(`--sample-id`, `--run-id`, `--output-root`) as the driver.

#### Run one package's CLI directly

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

#### Run a subset of that package's stages

Each `run` subcommand takes a `--stages` flag; pass a subset to run
only part of a step:

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

## Configuration

Each package ships a `config/default.yaml` with sensible defaults.
Overrides in precedence order (highest first):

1. CLI flag (e.g. `--x-source expected_counts`).
2. `--config <user.yaml>` — user YAML overriding default keys.
3. `config/default.yaml` — package-shipped defaults.

The full effective config for a run is recorded at
`<run-dir>/config.yaml`.

Per-step docs live under each package's `docs/` directory.

### Driver-level `--config` (per-run YAML)

`scripts/submit_workflow.sh` accepts `--config <path>`. The file gathers
the values you'd otherwise pass as CLI flags to the driver — no shared
folder needed; inputs point at any absolute path. CLI flags on the
same invocation win over the YAML, so the file is a defaults source you
can layer overrides on top of.

Recommended location: colocated with the run's outputs at
`<output_root>/<sample>/<sample>_<run_id>/config.yaml`.

Per-step keys use the semantic step names (`xenium_preprocess:`,
`ref_build:`, `rctd_split:`).

A ready-to-copy template listing every supported key with comments
lives at [`configs/example.yaml`](configs/example.yaml). Minimal shape:

```yaml
# runs/MH10/MH10_cap100/config.yaml
sample_id: MH10
run_id: cap100
output_root: /fh/fast/setty_m/user/ryang/workflow_runs

flex_h5ad: /fh/fast/setty_m/user/ryang/data/MH10_flex.h5ad
celltype_marker_json: /fh/fast/setty_m/user/ryang/data/markers.json

# Optional rctd-split explicit inputs — bypasses the run-folder layout
# auto-discovery. Use to mix a test_object from one sample with a
# reference from another.
test_object:   /fh/fast/setty_m/user/ryang/other/MH3_test_object.rds
reference_rds: /fh/fast/setty_m/user/ryang/refs/MH2_scRNA_ref.rds

xenium_preprocess:
  x_source: maxpost_counts
  qc_min_counts_cell: 10
ref_build:
  donor_borrow_cap: 100
rctd_split:
  umi_min: 10
  doublet_mode: doublet
```

```bash
./scripts/submit_workflow.sh --config runs/MH10/MH10_cap100/config.yaml
# CLI still wins — same YAML with a knob overridden:
./scripts/submit_workflow.sh --config runs/MH10/MH10_cap100/config.yaml \
    --ref-build-donor-borrow-cap 200
```

### Cell-type marker JSON

`ref-build` needs a JSON file declaring the celltypes it will build a
reference for, with one list of marker genes per celltype. The
pipeline doesn't ship one — it's data-dependent (which celltypes are
in your Flex scRNA and which genes read as markers in your panel).
Pass the path via `--celltype-marker-json` (driver flag), or under
`celltype_marker_json:` in a YAML config.

**Schema.** One top-level JSON object. Keys are `<CellType>_marker` or
`<CellType>_markers` (both accepted; trailing suffix stripped when
the census matches celltype names in `.obs`). Values are lists of
HGNC gene symbols. Example:

```json
{
  "B_marker":          ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1"],
  "Plasma_marker":     ["MZB1", "XBP1", "DERL3", "FKBP11", "TENT5C"],
  "Myeloid_marker":    ["CD14", "CD68", "MRC1", "CSF1R", "CD163", "SPI1"],
  "Fibroblast_marker": ["POSTN", "THY1", "PDGFRA", "CXCL12", "FAP", "COL5A1"],
  "T/NK_marker":       ["CD8A", "CD3D", "CD3E", "NKG7", "GZMB", "PRF1"],
  "endothelial_marker":["PECAM1", "VWF", "CDH5", "KDR", "ENG"],
  "RBC_markers":       ["HBB", "HBA1", "HBA2", "ALAS2"],
  "epithelial_markers":["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"]
}
```

The keys must cover every celltype present in your Flex data's
`.obs[celltype_col]` (default column name
`Final_level1_celltype_annotation`; override with
`--celltype-col-for-ref-build`). Missing celltypes get a `WARN` in
`census.csv` and are dropped from the reference. Extra celltypes
(present in the JSON but with zero cells across primary + donors +
fallback) get `decision=missing_no_donor` in the census with a
warning.

For a working reference on the shared drive:
`/fh/fast/setty_m/metx_liver_met/supplementary_data/marker_genes/markers_NonTumor_level1.json`
(8 celltypes: B, Plasma, Myeloid, Fibroblast, T/NK, endothelial,
RBC, epithelial). Copy + edit for your own celltype set.

## Outputs

All three steps write into a single run folder:

```
<output-root>/<sample>/<sample>_<run-id>/
├── spatial_adata/                        # xenium-preprocess (rctd-split augments in place)
│   ├── <sample>_proseg_raw.h5ad          # proseg cell×gene
│   ├── <sample>_xenium_ranger.h5ad       # xenium-ranger cell×gene
│   └── <sample>_proseg_purified.h5ad     # rctd-split — SPLIT-purified cell×gene
├── rctd/
│   ├── <sample>_test_object.rds          # xenium-preprocess — RCTD test object (spatial query)
│   ├── <sample>_reference.rds            # ref-build — RCTD reference (celltype pool)
│   └── <sample>_rctd_split.rds           # rctd-split — RCTD + SPLIT typing result
├── config.yaml                  # merged effective config across the three steps
└── logs/
    ├── slurm-<jobid>-xenium-preprocess.log   # slurm-captured stdout+stderr of the xenium-preprocess sbatch job
    ├── slurm-<jobid>-ref-build.log           # same, for ref-build
    ├── slurm-<jobid>-rctd-split.log          # same, for rctd-split
    └── workflow-submit.log                   # authoritative record of jobids + dep chain
```

Each `slurm-<jobid>-<stage>.log` is the FULL per-stage log — the stage's
Python + R shell-outs print progress lines to stdout, slurm captures both
stdout and stderr into that single file. There is no separate
`<stage>.log` sidecar; the sbatch stubs don't write one.

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
