# Installation

Full installation recipe for `xenium-preprocess-pipeline`.

## Prerequisites

- Linux with a Slurm scheduler (`sbatch`, `--dependency=afterok:` support).
- [Micromamba](https://mamba.readthedocs.io/) (or `conda`).
- [`uv`](https://docs.astral.sh/uv/) for the editable Python installs.
- R 4.4+ with Seurat, Matrix, spacexr, SPLIT, SpatialExperiment.

## 1. Clone the repository

```bash
git clone https://github.com/settylab/xenium-preprocess-pipeline
cd xenium-preprocess-pipeline
```

## 2. Create the Python environment

```bash
micromamba create -n xenium -f environments/xenium.yml
micromamba activate xenium
```

## 3. Install the three packages

Each package installs its own CLI entry point:

```bash
uv pip install -e packages/xenium-preprocess
uv pip install -e packages/ref-build
uv pip install -e packages/rctd-split
```

Verify:

```bash
xenium-preprocess --help
ref-build         --help
rctd-split        --help
```

## 4. R side

### On a cluster with Lmod

```bash
ml fhR/4.4.1-foss-2023b
Rscript -e '
  dir.create("~/R/x86_64-pc-linux-gnu-library/4.4", recursive = TRUE, showWarnings = FALSE)
  .libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("dmcable/spacexr")
  remotes::install_github("bdsc-tds/SPLIT")
'
```

### Off-cluster

Uncomment the R block in `environments/xenium.yml` and re-create the
env, or install R 4.4+ separately and install Seurat + Matrix + readr
+ SpatialExperiment via `install.packages()` /
`BiocManager::install()`. spacexr and SPLIT come from GitHub as above.

## 5. Verify the install

```bash
# Python side
pytest packages/xenium-preprocess/tests
pytest packages/ref-build/tests
pytest packages/rctd-split/tests

# Driver dry-run
./scripts/submit_workflow.sh --dry-run \
    --sample-id     TEST \
    --output-root   /tmp/workflow_runs \
    --run-id        v0 \
    --flex-h5ad     /tmp/nowhere.h5ad \
    --celltype-marker-json /tmp/nowhere.json \
    --proseg-dir    /tmp/nowhere \
    --xenium-ranger-dir /tmp/nowhere
```

The `--dry-run` prints the sbatch invocations without submitting.
