# Installation

`xenium-preprocess` has a Python side (stages 1-3) and an R side (stage 4). The two are activated separately: a conda/micromamba env for Python, and either an LMOD module or a conda-native R install for R.

## Python side

```bash
# 1. Install micromamba (skip if you already have it)
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# 2. Create the env
# (../../scripts/create-env.sh wraps `micromamba env create`; if
# MAMBA_ROOT_PREFIX is set for an isolated install, it also keeps the
# package cache isolated — see docs/installation.md § 2 at the repo root)
../../scripts/create-env.sh -f environment.yml
micromamba activate xeniumPreprocess

# 3. Editable install
pip install -e /path/to/xenium-preprocess

# 4. Verify
xenium-preprocess --help
xenium-preprocess run --help
```

The `environment.yml` in the repo carries the exact conda-side stack that stages 1-3 have been validated against.

## R side (for the `rctd_prep` stage)

Stage 4 shells out to `Rscript`. The R environment must carry:

- `Seurat` (≥ 5)
- `Matrix`
- `readr`
- `SpatialExperiment`
- `spacexr` (RCTD)

`spacexr` installs via `remotes::install_github()` (see the conda-native
path below), not CRAN — that needs outbound network access to
`github.com` / `api.github.com`. **No GitHub credentials are required**
— [`dmcable/spacexr`](https://github.com/dmcable/spacexr) is a public
repo — but GitHub's unauthenticated API is capped at 60 requests/hour
per source IP, easy to exhaust on shared infrastructure. If you already
have `GITHUB_PAT`/`GITHUB_TOKEN` set, or have ever run `gh auth login`,
`remotes` silently uses those credentials (raising the ceiling to
5,000/hour) with no indication the install would otherwise be running
against the smaller anonymous quota. See
[`../../../docs/installation.md`](../../../docs/installation.md) §
GitHub access for the full explanation.

### Recommended: Lmod R module

`fhR/4.4.1-foss-2023b` module carries all of the above:

```bash
ml fhR/4.4.1-foss-2023b
Rscript -e "library(Seurat); library(spacexr); cat('ok\n')"
```

`scripts/submit.slurm.sh` runs `ml $R_MODULE` (default: `fhR/4.4.1-foss-2023b`) before invoking the pipeline, so under Slurm you get this for free.

### Alternative: conda-native R

If you don't have LMOD access, uncomment the R-side block in `environment.yml`, re-run `../../scripts/create-env.sh`, then finish the install manually because `spacexr` is not on conda:

```r
BiocManager::install("SpatialExperiment")
remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)
```

Then point the pipeline at the conda-managed Rscript: `--rscript-bin "$CONDA_PREFIX/bin/Rscript"`.

## Verification

```bash
xenium-preprocess --version
xenium-preprocess run --help
# The R side is verified indirectly by a dry-run of stage 4 (it fails
# loud if Rscript is missing or if the R packages aren't installed).
```
