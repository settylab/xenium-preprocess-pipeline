# environments/

Environment specifications for `xenium-preprocess-pipeline`.

## `xenium.yml` — Python

A single micromamba (or conda) environment spec covering the Python
runtime for all three packages:

```bash
../scripts/create-env.sh -n xenium -f xenium.yml
micromamba activate xenium

uv pip install -e ../packages/xenium-preprocess
uv pip install -e ../packages/ref-build
uv pip install -e ../packages/rctd-split
```

## R

`ref-build` and `rctd-split` shell out to `Rscript`. Both need:

- `Seurat`
- `Matrix`
- `readr`
- `SpatialExperiment` (rctd-split)
- `spacexr` (RCTD)
- `SPLIT`

### On a cluster with Lmod

The `fhR/4.4.1-foss-2023b` module ships Seurat + Matrix +
SpatialExperiment + readr. `spacexr` and `SPLIT` are not in the module
and need a one-time user-library install:

```bash
ml fhR/4.4.1-foss-2023b
../scripts/install-r-packages.sh --r-lib-dir ~/R/x86_64-pc-linux-gnu-library/4.4
```

(`remotes::install_github()` checks the installed SHA across *all*
`.libPaths()` entries, not just the first — merely prepending a fresh
library, as an earlier version of this doc did, doesn't stop a
contaminated default library from silently absorbing the "install".
See [`../docs/installation.md`](../docs/installation.md) § R side for
why `scripts/install-r-packages.sh` exists.)

### Off-cluster

Install R 4.4+ separately, then install the packages listed above
either from the conda-native R block in `xenium.yml` (uncomment) or
from source. See [`../docs/installation.md`](../docs/installation.md).
