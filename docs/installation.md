# Installation

Full installation recipe for `xenium-preprocess-pipeline`.

> **BREAKING CHANGE for existing installs.** The `scripts/submit_*.sbatch` stubs no longer
> activate the Python env by `$HOME/.local/bin` + an env NAME, and no longer hardcode the R
> module. They now REQUIRE a gitignored `scripts/env.local.conf` (resolved absolute paths) and
> refuse to run without one. **Migration: run `scripts/write-env-config.sh` once after
> upgrading** (step 6 below) — that's the whole fix.
>
> **Why:** `$HOME` and a bare env NAME can resolve differently inside a Slurm job than in the
> shell that submitted it — that's what actually broke this pipeline (twice). `sbatch
> --export=ALL` forwards env VAR *values* faithfully, but does nothing to stop `$HOME` itself
> from resolving to a different instance once the job's own shell starts — so anything that
> re-derives a path from `$HOME`, or an env name, at job time is still exposed even with
> `--export=ALL`. A file of already-resolved absolute paths on shared storage is immune to that
> class of failure; that's the only thing that changed, and it's why there's no `$HOME`-derived
> fallback — a silent fallback would just reintroduce the same bug one layer down. See
> `scripts/lib/env_config.sh` for the full rationale and the exact failure mode this replaced.

## Prerequisites

- Linux with a Slurm scheduler (`sbatch`, `--dependency=afterok:` support).
- [Micromamba](https://mamba.readthedocs.io/) (or `conda`). The sbatch
  layer (`scripts/submit_*.sbatch`) activates the Python env by its
  ABSOLUTE PREFIX, resolved once by step 6 below into
  `scripts/env.local.conf` — it does not assume micromamba lives at any
  particular path (e.g. `$HOME/.local/bin`) or that a bare env NAME
  resolves the same way inside a Slurm job as it did in your login shell.
- [`uv`](https://docs.astral.sh/uv/) for the editable Python installs.
- R 4.4+ with Seurat, Matrix, spacexr, SPLIT, SpatialExperiment.

## 1. Clone the repository

```bash
git clone https://github.com/settylab/xenium-preprocess-pipeline
cd xenium-preprocess-pipeline
```

## 2. Create the Python environment

```bash
scripts/create-env.sh -n xenium -f environments/xenium.yml
micromamba activate xenium
```

`scripts/create-env.sh` is a thin wrapper around `micromamba create`.
If you don't set `MAMBA_ROOT_PREFIX`, it's a plain passthrough. If you
**do** set `MAMBA_ROOT_PREFIX` to keep this install fully isolated
(e.g. off `$HOME`, for a scratch/test install, or to keep multiple
installs from sharing state) — for example:

```bash
export MAMBA_ROOT_PREFIX=/abs/path/to/isolated/root
scripts/create-env.sh -n xenium -f environments/xenium.yml
```

— the wrapper also exports `CONDA_PKGS_DIRS="$MAMBA_ROOT_PREFIX/pkgs"`
and asserts afterward that `~/.mamba/pkgs` was not touched. Without
this, micromamba's `pkgs_dirs` silently resolves to
`[$MAMBA_ROOT_PREFIX/pkgs, ~/.mamba/pkgs]` — an undocumented second
entry — so an "isolated" install can still quietly write package-cache
state to `~/.mamba/pkgs`. Nothing errors when that happens; the
wrapper's post-create check exists because this failure is otherwise
invisible. A bare `micromamba create` still works exactly as before if
you don't need isolation.

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

The `pytest` suites need the `[test]` extra — it's gated behind
`[project.optional-dependencies]` in each package's `pyproject.toml` and
is NOT part of `environments/xenium.yml`, so install it explicitly
before running bare `pytest` (omitting this step is the most common
"why does `pytest` fail with `ModuleNotFoundError: No module named
'pytest'`" report):

```bash
uv pip install -e "packages/xenium-preprocess[test]"
uv pip install -e "packages/ref-build[test]"
uv pip install -e "packages/rctd-split[test]"

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

## 6. Record the resolved environment

The final install step: resolve the micromamba binary, the env
**prefix** (not name), the R library dir from step 4, and the R module
into absolute values, and write them to the gitignored
`scripts/env.local.conf`. Every `scripts/submit_*.sbatch` stub sources
this instead of re-deriving `$HOME`-relative paths inside the Slurm
job — `$HOME` can resolve to a different location inside a job than in
your login shell, which is exactly what broke this pipeline twice
before (a missing R library, and micromamba's own root-prefix
resolution) — see `scripts/lib/env_config.sh` for the full rationale.

```bash
./scripts/write-env-config.sh \
    --env-prefix "$(micromamba env list | awk '/xenium/ {print $NF}')" \
    --r-lib-dir  ~/R/x86_64-pc-linux-gnu-library/4.4
```

Then verify it end-to-end (also runs automatically before every
`submit_workflow.sh` dispatch — `--skip-preflight` to opt out):

```bash
./scripts/env-preflight.sh
```
