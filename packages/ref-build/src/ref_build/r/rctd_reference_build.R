#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# rctd_reference_build.R -- Stage C of the ref-build pipeline.
#
# Reads the 10X-style mtx bundle written by ref_build.stages.export_mtx
# and produces a spacexr::Reference object saved as an RDS. That RDS
# is the RCTD-ready reference — downstream `spacexr::create.RCTD(spatial,
# reference)` calls take it as the reference side.
#
# Recipe ported from ref-build-summary-v3.md Stage C (lines 252-281),
# and from the internal SPLIT/Proseg workflow's rebuild scripts under
#   an internal source path
#   an internal source path
#
# Invocation (via ref_build.stages.rctd_reference_build):
#   Rscript --vanilla rctd_reference_build.R \
#       --mtx-dir=<path> \
#       --sample-id=<sample_id> \
#       --out-rds=<output.rds> \
#       --min-umi=10 \
#       --require-int=TRUE \
#       --celltype-col=Final_level1_celltype_annotation \
#       --label-slash-replacement=_
#
# Files consumed under <mtx-dir>:
#   <sample_id>_counts.mtx.gz
#   <sample_id>_features.tsv.gz
#   <sample_id>_barcodes.tsv.gz
#   <sample_id>_metadata.csv
# ---------------------------------------------------------------------

# ----- Prepend the user-local R lib so spacexr is found ---------------
# The default fhR/4.4.1-foss-2023b module carries Seurat + Matrix + readr
# but NOT spacexr; the ref-build-summary v3 lines 279-281 describe the
# one-time user-local install at `~/.claude/r_libs/4.4.1`. Adding it to
# .libPaths() here is idempotent: if the dir doesn't exist, we skip it.
user_lib <- file.path(Sys.getenv("HOME"), ".claude", "r_libs", "4.4.1")
if (dir.exists(user_lib)) {
  .libPaths(c(user_lib, .libPaths()))
  cat(sprintf("[rctd_reference_build.R] prepended user lib: %s\n", user_lib))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(spacexr)
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

need <- c(
  "mtx-dir", "sample-id", "out-rds",
  "min-umi", "require-int",
  "celltype-col", "label-slash-replacement"
)
missing <- setdiff(need, names(kv))
if (length(missing)) {
  stop("[rctd_reference_build.R] missing required args: ",
       paste(missing, collapse = ", "))
}

mtx_dir      <- kv[["mtx-dir"]]
sample_id    <- kv[["sample-id"]]
out_rds      <- kv[["out-rds"]]
min_umi      <- as.integer(kv[["min-umi"]])
require_int  <- toupper(kv[["require-int"]]) %in% c("TRUE", "T", "YES", "1")
celltype_col <- kv[["celltype-col"]]
slash_repl   <- kv[["label-slash-replacement"]]

cat(sprintf("[rctd_reference_build.R] mtx_dir      = %s\n", mtx_dir))
cat(sprintf("[rctd_reference_build.R] sample_id    = %s\n", sample_id))
cat(sprintf("[rctd_reference_build.R] out_rds      = %s\n", out_rds))
cat(sprintf("[rctd_reference_build.R] min_umi      = %d\n", min_umi))
cat(sprintf("[rctd_reference_build.R] require_int  = %s\n",
            ifelse(require_int, "TRUE", "FALSE")))
cat(sprintf("[rctd_reference_build.R] celltype_col = %s\n", celltype_col))
cat(sprintf("[rctd_reference_build.R] slash_repl   = %s\n", slash_repl))

# ----- Resolve input paths -------------------------------------------
# Prefer .gz over uncompressed; fall back to the plain file if gzip
# outputs were disabled in export_mtx.
resolve_path <- function(stem, suffix_gz, suffix_plain) {
  gz <- file.path(mtx_dir, paste0(stem, suffix_gz))
  pl <- file.path(mtx_dir, paste0(stem, suffix_plain))
  if (file.exists(gz)) return(gz)
  if (file.exists(pl)) return(pl)
  stop("[rctd_reference_build.R] missing input file: neither ", gz, " nor ", pl)
}

mtx_path      <- resolve_path(sample_id, "_counts.mtx.gz", "_counts.mtx")
feats_path    <- resolve_path(sample_id, "_features.tsv.gz", "_features.tsv")
barcodes_path <- resolve_path(sample_id, "_barcodes.tsv.gz", "_barcodes.tsv")
meta_path     <- file.path(mtx_dir, paste0(sample_id, "_metadata.csv"))
if (!file.exists(meta_path)) {
  stop("[rctd_reference_build.R] missing metadata: ", meta_path)
}

# ----- Load the 10X mtx bundle ---------------------------------------
cat(sprintf("[rctd_reference_build.R] reading mtx = %s\n", mtx_path))
mtx <- ReadMtx(
  mtx            = mtx_path,
  features       = feats_path,
  cells          = barcodes_path,
  feature.column = 1,   # single-column features file (see export_mtx.py)
  cell.column    = 1
)
cat(sprintf("[rctd_reference_build.R]   mtx shape: %d genes x %d cells\n",
            nrow(mtx), ncol(mtx)))

# ----- Load + join metadata ------------------------------------------
cat(sprintf("[rctd_reference_build.R] reading metadata = %s\n", meta_path))
meta <- read.csv(meta_path, row.names = 1, check.names = FALSE,
                 stringsAsFactors = FALSE)

if (!(celltype_col %in% colnames(meta))) {
  stop("[rctd_reference_build.R] celltype column '", celltype_col,
       "' not found in metadata; columns present: ",
       paste(colnames(meta), collapse = ", "))
}

# Align the metadata row order to the mtx column order.
common <- intersect(colnames(mtx), rownames(meta))
if (length(common) != ncol(mtx)) {
  stop("[rctd_reference_build.R] cell-id mismatch between mtx barcodes and ",
       "metadata rownames: ", length(common), " common vs ", ncol(mtx),
       " mtx columns.")
}
meta  <- meta[colnames(mtx), , drop = FALSE]
labels_raw <- meta[[celltype_col]]

# ----- Sanitise cell-type labels -------------------------------------
# spacexr factor levels reject '/' — replace with the configured
# replacement (default `_`). Example: B/Plasma_T/NK_rbc -> B_Plasma_T_NK_rbc.
labels <- gsub("/", slash_repl, labels_raw, fixed = TRUE)
if (any(labels != labels_raw)) {
  changed <- which(labels != labels_raw)
  n_changed <- length(unique(labels_raw[changed]))
  cat(sprintf("[rctd_reference_build.R] sanitised %d unique labels containing '/'\n",
              n_changed))
}

ref_labels <- factor(labels)
names(ref_labels) <- colnames(mtx)
cat(sprintf("[rctd_reference_build.R]   celltype set (%d): %s\n",
            nlevels(ref_labels),
            paste(levels(ref_labels), collapse = ", ")))

# ----- Build the spacexr Reference -----------------------------------
# Ported from ref-build-summary-v3.md lines 264-272 and the user's rebuild
# scripts (an internal source path, an internal source path). We
# construct the counts matrix straight from `mtx` (already a dgCMatrix
# from ReadMtx) — no need to detour through a Seurat object first, and
# it saves a couple GB on large primaries.
cat("[rctd_reference_build.R] calling spacexr::Reference...\n")
ref <- Reference(
  counts      = mtx,
  cell_types  = ref_labels,
  min_UMI     = min_umi,
  require_int = require_int
)

# ----- Save ----------------------------------------------------------
dir.create(dirname(out_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(ref, file = out_rds)

cat(sprintf("[rctd_reference_build.R] wrote %s\n", out_rds))
