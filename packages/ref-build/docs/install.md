# Installation

`ref-build` has a Python side (stages 1-4) and an R side (stage 5). The two are activated separately: a conda/micromamba env for Python, and either an LMOD module or a conda-native R install for R.

## Python side

```bash
# 1. Install micromamba (skip if you already have it)
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# 2. Create the env
# (../../scripts/create-env.sh wraps `micromamba env create`; if
# MAMBA_ROOT_PREFIX is set for an isolated install, it also keeps the
# package cache isolated — see docs/installation.md § 2 at the repo root)
../../scripts/create-env.sh -f environment.yml
micromamba activate refBuild

# 3. Editable install
pip install -e /path/to/ref-build

# 4. Verify
ref-build --help
ref-build run --help
```

The `environment.yml` in the repo carries the exact conda-side stack that stages 1-4 have been validated against.

## R side (for the `rctd_reference_build` stage)

Stage 5 shells out to `Rscript`. The R environment must carry:

- `Seurat` (≥ 5)
- `Matrix`
- `spacexr` (RCTD; **not** shipped by the LMOD module — one-time install below)

`spacexr` installs via `remotes::install_github()` below, not CRAN —
that needs outbound network access to `github.com` / `api.github.com`.
**No GitHub credentials are required** —
[`dmcable/spacexr`](https://github.com/dmcable/spacexr) is a public
repo — but GitHub's unauthenticated API is capped at 60 requests/hour
per source IP, easy to exhaust on shared infrastructure. If you already
have `GITHUB_PAT`/`GITHUB_TOKEN` set, or have ever run `gh auth login`,
`remotes` silently uses those credentials (raising the ceiling to
5,000/hour) with no indication the install would otherwise be running
against the smaller anonymous quota. See
[`../../../docs/installation.md`](../../../docs/installation.md) §
GitHub access for the full explanation.

### Recommended: Lmod R module + user-local `spacexr`

`fhR/4.4.1-foss-2023b` module carries `Seurat` + `Matrix`. `spacexr` is not on the module; the ref-build recipe (see `the internal reference summary` lines 279-281) installs it once into a user-local R library at `~/.claude/r_libs/4.4.1` and reuses it across builds:

```bash
ml fhR/4.4.1-foss-2023b
mkdir -p ~/.claude/r_libs/4.4.1
Rscript -e '.libPaths(c("~/.claude/r_libs/4.4.1", .libPaths())); \
            if (!requireNamespace("remotes", quietly=TRUE)) \
                install.packages("remotes", lib="~/.claude/r_libs/4.4.1", \
                                 repos="https://cloud.r-project.org"); \
            remotes::install_github("dmcable/spacexr", \
                                    lib="~/.claude/r_libs/4.4.1", \
                                    upgrade="never", build_vignettes=FALSE)'
Rscript -e '.libPaths(c("~/.claude/r_libs/4.4.1", .libPaths())); library(spacexr); cat("spacexr ok\n")'
```

The `ref-build/src/ref_build/r/rctd_reference_build.R` script prepends `~/.claude/r_libs/4.4.1` to `.libPaths()` on startup (idempotent — if the dir doesn't exist, the prepend is skipped), so subsequent stage-5 runs pick up `spacexr` automatically once installed. See the R script header for details.

`scripts/submit.slurm.sh` runs `ml $R_MODULE` (default: `fhR/4.4.1-foss-2023b`) before invoking the pipeline, so under Slurm you get the module for free — you still need the one-time `spacexr` install above.

### Alternative: conda-native R

If you don't have LMOD access, uncomment the R-side block in `environment.yml`, re-run `../../scripts/create-env.sh`, then finish the install manually because `spacexr` is not on conda:

```r
remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)
```

Then point the pipeline at the conda-managed Rscript: `--rscript-bin "$CONDA_PREFIX/bin/Rscript"`.

## Verification

```bash
ref-build --version
ref-build run --help
# Stage 5 is verified indirectly by a dry-run (Rscript is required to
# exist on PATH; missing Rscript fails LOUD before the R subprocess
# starts).
```
