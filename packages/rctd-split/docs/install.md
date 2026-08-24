# Installation

`rctd-split` has a Python side (stage 4, `mtx_to_h5ad`, plus the three shell-out wrappers) and an R side (stages 1-3, all shelling out to `Rscript`). The two are activated separately: a conda/micromamba env for Python, and either an LMOD module + user-local R lib or a conda-native R install for R.

## Python side

```bash
# 1. Install micromamba (skip if you already have it)
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# 2. Create the env
# (../../scripts/create-env.sh wraps `micromamba env create`; if
# MAMBA_ROOT_PREFIX is set for an isolated install, it also keeps the
# package cache isolated — see docs/installation.md § 2 at the repo root)
../../scripts/create-env.sh -f environment.yml
micromamba activate rctdSplit

# 3. Editable install
pip install -e /path/to/rctd-split

# 4. Verify
rctd-split --help
rctd-split run --help
```

The `environment.yml` in the repo carries the exact conda-side stack that the Python side has been validated against (`anndata=0.11`, `scanpy=1.11`).

## R side (for stages 1-3: rctd_run, split_purify, export_mtx)

The R environment must carry:

- `Seurat` (≥ 5)
- `Matrix`
- `spacexr` (RCTD)
- `SPLIT` (bdsc-tds/SPLIT)

`spacexr` and `SPLIT` both install via `remotes::install_github()`
below, not CRAN — that needs outbound network access to `github.com` /
`api.github.com`. **No GitHub credentials are required** — both source
repos ([`dmcable/spacexr`](https://github.com/dmcable/spacexr),
[`bdsc-tds/SPLIT`](https://github.com/bdsc-tds/SPLIT)) are public — but
GitHub's unauthenticated API is capped at 60 requests/hour per source
IP, easy to exhaust on shared infrastructure. If you already have
`GITHUB_PAT`/`GITHUB_TOKEN` set, or have ever run `gh auth login`,
`remotes` silently uses those credentials (raising the ceiling to
5,000/hour) with no indication the install would otherwise be running
against the smaller anonymous quota. See
[`../../../docs/installation.md`](../../../docs/installation.md) §
GitHub access for the full explanation.

### Recommended: Lmod R module + user-local R lib

`fhR/4.4.1-foss-2023b` module carries `Seurat` + `Matrix`, but NOT `spacexr` or `SPLIT`. Both need a one-time user-local install at `~/.claude/r_libs/4.4.1`:

```bash
ml fhR/4.4.1-foss-2023b

Rscript --vanilla -e '
  lib <- file.path(Sys.getenv("HOME"), ".claude", "r_libs", "4.4.1")
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
  options(timeout = 600, download.file.method = "libcurl")
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", lib = lib)
  remotes::install_github("dmcable/spacexr",
                          build_vignettes = FALSE, upgrade = "never", lib = lib)
  remotes::install_github("bdsc-tds/SPLIT",
                          build_vignettes = FALSE, upgrade = "never", lib = lib)
'

Rscript -e '
  lib <- file.path(Sys.getenv("HOME"), ".claude", "r_libs", "4.4.1")
  .libPaths(c(lib, .libPaths()))
  library(Seurat)
  library(spacexr)
  library(SPLIT)
  cat("ok\n")
'
```

Each R script this pipeline ships (`r/rctd_run.R`, `r/split_purify.R`, `r/export_mtx.R`) prepends `~/.claude/r_libs/4.4.1` to `.libPaths()` on start, so once the install above is done, invocations transparently pick up spacexr + SPLIT.

`scripts/submit.slurm.sh` runs `ml $R_MODULE` (default: `fhR/4.4.1-foss-2023b`) before invoking the pipeline, so under Slurm you get this for free.

### Alternative: conda-native R

If you don't have LMOD access, uncomment the R-side block in `environment.yml`, re-run `../../scripts/create-env.sh`, then finish the install manually because `spacexr` and `SPLIT` are not on conda:

```r
BiocManager::install("SpatialExperiment")   # optional; not required by any R script
remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)
remotes::install_github("bdsc-tds/SPLIT",  build_vignettes = FALSE)
```

Then point the pipeline at the conda-managed Rscript: `--rscript-bin "$CONDA_PREFIX/bin/Rscript"`.

## Verification

```bash
rctd-split --version
rctd-split run --help
# The R side is verified indirectly by a dry-run of stage 1 (it fails
# loud if Rscript is missing or if spacexr / SPLIT aren't installed).
```
