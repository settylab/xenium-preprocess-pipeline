#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# export_mtx.R -- Stage 3 of the rctd-split pipeline.
#
# Reads unpurified.rds + purified.rds (from stage 2) and writes matched
# 10X-style bundles per variant:
#   <out-dir>/<sample>_{unpurified,purified}_counts.mtx.gz
#   <out-dir>/<sample>_{unpurified,purified}_features.tsv.gz
#   <out-dir>/<sample>_{unpurified,purified}_barcodes.tsv.gz
#   <out-dir>/<sample>_{unpurified,purified}_metadata.csv
#   <out-dir>/<sample>_{unpurified,purified}_spatial_coords.csv.gz
#
# Same filename convention as step 1's split_prep and step 3's
# export_mtx.
#
# Invocation (via rctd_split.stages.export_mtx):
#   Rscript --vanilla export_mtx.R \
#       --sample-id=<sample> \
#       --unpurified-rds=<path> \
#       --purified-rds=<path> \
#       --unpurified-out-dir=<path> \
#       --purified-out-dir=<path> \
#       --assay-name=Proseg \
#       --gzip-outputs=TRUE
# ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# ----- Arg parsing ---------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
kv <- list()
for (a in args) {
  if (grepl("^--[^=]+=", a)) {
    k <- sub("^--", "", sub("=.*$", "", a))
    v <- sub("^[^=]+=", "", a)
    kv[[k]] <- v
  }
}

need <- c(
  "sample-id",
  "unpurified-rds", "purified-rds",
  "unpurified-out-dir", "purified-out-dir",
  "assay-name", "gzip-outputs"
)
missing <- setdiff(need, names(kv))
if (length(missing)) {
  stop("[export_mtx.R] missing required args: ",
       paste(missing, collapse = ", "))
}

sample_id          <- kv[["sample-id"]]
unpurified_rds     <- kv[["unpurified-rds"]]
purified_rds       <- kv[["purified-rds"]]
unpurified_out_dir <- kv[["unpurified-out-dir"]]
purified_out_dir   <- kv[["purified-out-dir"]]
assay_name         <- kv[["assay-name"]]
gzip_outputs       <- toupper(kv[["gzip-outputs"]]) %in% c("TRUE", "T", "YES", "1")

cat(sprintf("[export_mtx.R] sample_id           = %s\n", sample_id))
cat(sprintf("[export_mtx.R] unpurified_rds      = %s\n", unpurified_rds))
cat(sprintf("[export_mtx.R] purified_rds        = %s\n", purified_rds))
cat(sprintf("[export_mtx.R] unpurified_out_dir  = %s\n", unpurified_out_dir))
cat(sprintf("[export_mtx.R] purified_out_dir    = %s\n", purified_out_dir))
cat(sprintf("[export_mtx.R] assay_name          = %s\n", assay_name))
cat(sprintf("[export_mtx.R] gzip_outputs        = %s\n",
            ifelse(gzip_outputs, "TRUE", "FALSE")))

# ----- Helpers -------------------------------------------------------
# Write a matrix to `path` as MatrixMarket, gzipped when requested. Base
# R's writeMM has no built-in gzip; we route via a temp file + gzip when
# needed (matches step-1 split_prep's on-disk shape).
write_mm_maybe_gz <- function(mat, path, gz) {
  if (gz) {
    tmp <- tempfile(fileext = ".mtx")
    writeMM(mat, file = tmp)
    # Read the plain file back, write gzipped.
    con_in <- file(tmp, "rb")
    con_out <- gzfile(path, "wb")
    on.exit({ close(con_in); close(con_out); unlink(tmp) }, add = TRUE)
    while (length(chunk <- readBin(con_in, "raw", n = 65536))) {
      writeBin(chunk, con_out)
    }
  } else {
    writeMM(mat, file = path)
  }
}

write_lines_maybe_gz <- function(vec, path, gz) {
  if (gz) {
    con <- gzfile(path, "wt")
  } else {
    con <- file(path, "wt")
  }
  on.exit(close(con), add = TRUE)
  # Single-column form; matches Seurat::ReadMtx(feature.column=1).
  writeLines(as.character(vec), con)
}

write_coords_maybe_gz <- function(df, path, gz) {
  if (gz) {
    write.csv(df, gzfile(path), row.names = TRUE)
  } else {
    write.csv(df, path, row.names = TRUE)
  }
}

# ----- Emit one bundle from a Seurat RDS -----------------------------
emit_bundle <- function(rds_path, out_dir, variant) {
  cat(sprintf("[export_mtx.R] --- %s bundle ---\n", variant))
  cat(sprintf("[export_mtx.R] reading %s\n", rds_path))
  seu <- readRDS(rds_path)
  if (!(assay_name %in% Assays(seu))) {
    stop("[export_mtx.R] assay '", assay_name, "' not present in ", rds_path,
         "; available: ", paste(Assays(seu), collapse = ", "))
  }

  counts_mat <- as(GetAssayData(seu, assay = assay_name, layer = "counts"),
                   "CsparseMatrix")
  cat(sprintf("[export_mtx.R]   %s: %d genes x %d cells\n",
              variant, nrow(counts_mat), ncol(counts_mat)))

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  gz_suffix <- if (gzip_outputs) ".gz" else ""
  stem <- file.path(out_dir, sprintf("%s_%s", sample_id, variant))

  mtx_path      <- paste0(stem, "_counts.mtx", gz_suffix)
  feats_path    <- paste0(stem, "_features.tsv", gz_suffix)
  barcodes_path <- paste0(stem, "_barcodes.tsv", gz_suffix)
  meta_path     <- paste0(stem, "_metadata.csv")
  coords_path   <- paste0(stem, "_spatial_coords.csv", gz_suffix)

  # Counts (genes x cells).
  write_mm_maybe_gz(counts_mat, mtx_path, gzip_outputs)
  cat(sprintf("[export_mtx.R]   wrote %s\n", mtx_path))

  # Features + barcodes.
  write_lines_maybe_gz(rownames(counts_mat), feats_path, gzip_outputs)
  cat(sprintf("[export_mtx.R]   wrote %s\n", feats_path))
  write_lines_maybe_gz(colnames(counts_mat), barcodes_path, gzip_outputs)
  cat(sprintf("[export_mtx.R]   wrote %s\n", barcodes_path))

  # Metadata (cell-indexed CSV, matches split_prep / ref-build's shape).
  meta <- seu@meta.data
  meta <- meta[colnames(counts_mat), , drop = FALSE]
  # Plain CSV (not gzipped) — matches the sibling pipelines.
  write.csv(meta, meta_path, row.names = TRUE)
  cat(sprintf("[export_mtx.R]   wrote %s\n", meta_path))

  # Spatial coords — pulled from top-level metadata cols x, y (which
  # split_purify.R keeps on both variants). If absent (rare — e.g. a
  # purified variant that lost x/y in a bespoke re-run), we skip with
  # a warning rather than fail.
  if (all(c("x", "y") %in% colnames(meta))) {
    coords <- data.frame(
      x = meta[, "x"],
      y = meta[, "y"],
      row.names = rownames(meta),
      check.names = FALSE
    )
    write_coords_maybe_gz(coords, coords_path, gzip_outputs)
    cat(sprintf("[export_mtx.R]   wrote %s\n", coords_path))
  } else {
    cat(sprintf("[export_mtx.R]   WARN: %s has no x/y in metadata; ",
                variant))
    cat("skipping spatial_coords export.\n")
  }
}

emit_bundle(unpurified_rds, unpurified_out_dir, "unpurified")
emit_bundle(purified_rds,   purified_out_dir,   "purified")

cat("[export_mtx.R] done.\n")
