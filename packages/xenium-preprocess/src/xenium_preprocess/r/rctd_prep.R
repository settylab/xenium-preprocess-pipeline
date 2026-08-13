#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# rctd_prep.R -- build the RCTD "test object" (Seurat spatial) from the
# mtx/features/barcodes triple + metadata + spatial-coords sidecars
# written by xenium_preprocess.stages.split_prep.
#
# Adapted verbatim from the first half of
#   the internal spatial-RCTD preparation reference
#
# The Rmd's second half (building a spacexr Reference from a scRNA
# mtx triple) is OUT OF SCOPE for step 1 -- the scRNA reference is
# user data, not derived from proseg. Later steps will land that.
#
# Invocation (via xenium_preprocess.stages.rctd_prep):
#   Rscript --vanilla rctd_prep.R \
#       --split-dir=<path> \
#       --stem=<sample_id><name_suffix> \
#       --out-rds=<output.rds> \
#       --assay-name=Proseg \
#       --spatial-key=ST_
#
# Files consumed under <split-dir>:
#   <stem>_counts.mtx.gz
#   <stem>_features.tsv.gz
#   <stem>_barcodes.tsv.gz
#   <stem>_metadata.csv
#   <stem>_spatial_coords.csv.gz
# ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(readr)
  library(SpatialExperiment)
})

# ----- Arg parsing (mini kv style; no CRAN dep) ----------------------
args <- commandArgs(trailingOnly = TRUE)
kv <- list()
for (a in args) {
  if (grepl("^--[^=]+=", a)) {
    k <- sub("^--", "", sub("=.*$", "", a))
    v <- sub("^[^=]+=", "", a)
    kv[[k]] <- v
  }
}

need <- c("split-dir", "stem", "out-rds", "assay-name", "spatial-key")
missing <- setdiff(need, names(kv))
if (length(missing)) {
  stop("[rctd_prep.R] missing required args: ",
       paste(missing, collapse = ", "))
}

split_dir    <- kv[["split-dir"]]
stem         <- kv[["stem"]]
out_rds      <- kv[["out-rds"]]
assay_name   <- kv[["assay-name"]]
spatial_key  <- kv[["spatial-key"]]

cat(sprintf("[rctd_prep.R] split_dir   = %s\n", split_dir))
cat(sprintf("[rctd_prep.R] stem        = %s\n", stem))
cat(sprintf("[rctd_prep.R] out_rds     = %s\n", out_rds))
cat(sprintf("[rctd_prep.R] assay_name  = %s\n", assay_name))
cat(sprintf("[rctd_prep.R] spatial_key = %s\n", spatial_key))

fp <- function(suffix) file.path(split_dir, paste0(stem, suffix))
meta_path    <- fp("_metadata.csv")
coords_path  <- fp("_spatial_coords.csv.gz")
mtx_path     <- fp("_counts.mtx.gz")
feats_path   <- fp("_features.tsv.gz")
barcodes_path <- fp("_barcodes.tsv.gz")

for (p in c(meta_path, coords_path, mtx_path, feats_path, barcodes_path)) {
  if (!file.exists(p)) {
    stop("[rctd_prep.R] missing input file: ", p)
  }
}

# ----- Load the SPLIT triple + sidecars (mirrors the Rmd) ------------
meta   <- read.csv(meta_path, row.names = 1, check.names = FALSE)
coords <- as.matrix(readr::read_csv(coords_path, show_col_types = FALSE))

mtx <- ReadMtx(
  mtx           = mtx_path,
  features      = feats_path,
  cells         = barcodes_path,
  feature.column = 1,   # single-column features file (see split_prep.py)
  cell.column    = 1
)

# ----- Build the SpatialExperiment (matches the Rmd) -----------------
spe <- SpatialExperiment(
  assays  = list(counts = mtx),
  colData = DataFrame(meta[colnames(mtx), , drop = FALSE])
)

spatialCoords(spe) <- as.matrix(coords)

# Keep only the (x, y) columns for the spatial coords. The reference
# Rmd assumes the coords CSV has an id column + x + y; readr::read_csv
# reads the id into the first column, so we hard-select "x" and "y".
if (all(c("x", "y") %in% colnames(spatialCoords(spe)))) {
  spatialCoords(spe) <- spatialCoords(spe)[, c("x", "y"), drop = FALSE]
}

# ----- Convert to Seurat + attach spatial DimReduc (matches the Rmd) -
seu <- CreateSeuratObject(
  counts    = counts(spe),
  assay     = assay_name,
  meta.data = as.data.frame(colData(spe))
)

sp_coords <- spatialCoords(spe)
# Match the Rmd's rename to ST_1 / ST_2 for the DimReduc.
colnames(sp_coords) <- c(paste0(spatial_key, 1), paste0(spatial_key, 2))
rownames(sp_coords) <- colnames(seu)

seu[["spatial"]] <- CreateDimReducObject(
  sp_coords, assay = assay_name, key = spatial_key
)

seu$x <- sp_coords[, 1]
seu$y <- sp_coords[, 2]

# ----- Save (drop the SPE to keep the RDS small) --------------------
rm(spe)
gc()

dir.create(dirname(out_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(seu, file = out_rds)

cat(sprintf("[rctd_prep.R] wrote %s (cells=%d, features=%d)\n",
            out_rds, ncol(seu), nrow(seu)))
