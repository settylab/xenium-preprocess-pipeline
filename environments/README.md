# environments/

Environment specifications for `xenium-preprocess-pipeline`.

## `xenium.yml` — Python

A single micromamba (or conda) environment spec covering the Python
runtime for all three packages:

```bash
micromamba create -n xenium -f environments/xenium.yml
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
Rscript -e '
  dir.create("~/R/x86_64-pc-linux-gnu-library/4.4", recursive = TRUE, showWarnings = FALSE)
  .libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("dmcable/spacexr")
  remotes::install_github("bdsc-tds/SPLIT")
'
```

### Off-cluster

Install R 4.4+ separately, then install the packages listed above
either from the conda-native R block in `xenium.yml` (uncomment) or
from source. See [`../docs/installation.md`](../docs/installation.md).
