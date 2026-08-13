#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# split_purify.R -- Stage 2 of the rctd-split pipeline.
#
# Reads rctd_results.rds (from stage 1) + the spatial test object RDS
# (step 1 xenium-preprocess output), applies SPLIT::run_post_process_RCTD
# to produce the "unpurified" variant (all cells + RCTD calls, before
# purification), and SPLIT::purify to produce the "purified" variant
# (confident subset).
#
# Ported verbatim from an internal SPLIT/Proseg workflow
#   the internal SPLIT/Proseg workflow
# lines 494-523.
#
# Invocation (via rctd_split.stages.split_purify):
#   Rscript --vanilla split_purify.R \
#       --test-object=<path> \
#       --rctd-results-rds=<path> \
#       --out-unpurified-rds=<path> \
#       --out-purified-rds=<path> \
#       --assay-name=Proseg \
#       --do-purify-singlets=TRUE
# ---------------------------------------------------------------------

# ----- Prepend the user-local R lib so spacexr + SPLIT are found -----
user_lib <- file.path(Sys.getenv("HOME"), ".claude", "r_libs", "4.4.1")
if (dir.exists(user_lib)) {
  .libPaths(c(user_lib, .libPaths()))
  cat(sprintf("[split_purify.R] prepended user lib: %s\n", user_lib))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(spacexr)
  library(SPLIT)
})

# Source the meta-preserve helper from the same directory as this script.
# (Rscript sets --file=<path> on commandArgs so we can locate ourselves.)
this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(trailingOnly = FALSE),
                      value = TRUE))
if (length(this_file) == 0L) {
  stop("[split_purify.R] cannot determine script path from commandArgs()")
}
source(file.path(dirname(normalizePath(this_file)), "_preserve_meta.R"))

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
  "test-object", "rctd-results-rds",
  "out-unpurified-rds", "out-purified-rds",
  "assay-name", "do-purify-singlets"
)
missing <- setdiff(need, names(kv))
if (length(missing)) {
  stop("[split_purify.R] missing required args: ",
       paste(missing, collapse = ", "))
}

test_object_path   <- kv[["test-object"]]
rctd_results_rds   <- kv[["rctd-results-rds"]]
out_unpurified_rds <- kv[["out-unpurified-rds"]]
out_purified_rds   <- kv[["out-purified-rds"]]
assay_name         <- kv[["assay-name"]]
DO_purify_singlets <- toupper(kv[["do-purify-singlets"]]) %in% c("TRUE", "T", "YES", "1")

cat(sprintf("[split_purify.R] test_object        = %s\n", test_object_path))
cat(sprintf("[split_purify.R] rctd_results_rds   = %s\n", rctd_results_rds))
cat(sprintf("[split_purify.R] out_unpurified_rds = %s\n", out_unpurified_rds))
cat(sprintf("[split_purify.R] out_purified_rds   = %s\n", out_purified_rds))
cat(sprintf("[split_purify.R] assay_name         = %s\n", assay_name))
cat(sprintf("[split_purify.R] DO_purify_singlets = %s\n",
            ifelse(DO_purify_singlets, "TRUE", "FALSE")))

# ----- Load inputs ---------------------------------------------------
cat("[split_purify.R] reading test object...\n")
xe <- readRDS(test_object_path)
cat(sprintf("[split_purify.R]   test object: cells=%d, features=%d\n",
            ncol(xe), nrow(xe)))

cat("[split_purify.R] reading RCTD results...\n")
rctd_result <- readRDS(rctd_results_rds)

# ----- SPLIT::run_post_process_RCTD (Rmd lines 495-497) --------------
# The Rmd:
#   rctd_result <- SPLIT::run_post_process_RCTD(rctd_result)
#   xe <- AddMetaData(xe, rctd_result@results$results_df)
cat("[split_purify.R] SPLIT::run_post_process_RCTD...\n")
rctd_result <- SPLIT::run_post_process_RCTD(rctd_result)

results_df <- rctd_result@results$results_df
cat(sprintf("[split_purify.R]   results_df: %d rows x %d cols\n",
            nrow(results_df), ncol(results_df)))
cat(sprintf("[split_purify.R]   columns: %s\n",
            paste(colnames(results_df), collapse = ", ")))

# Attach the per-cell RCTD calls to the spatial Seurat's metadata.
# This is the "unpurified" variant — all cells that survived the
# spatial-side QC + RCTD's UMI_min gate, each now carrying first_type,
# second_type, spot_class, and the SPLIT-emitted columns
# (purification_status, w1_larger_w2, same_class).
xe_unpurified <- AddMetaData(xe, results_df)
cat(sprintf("[split_purify.R]   unpurified variant: cells=%d, features=%d\n",
            ncol(xe_unpurified), nrow(xe_unpurified)))

# ----- Save unpurified.rds ------------------------------------------
dir.create(dirname(out_unpurified_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(xe_unpurified, file = out_unpurified_rds)
cat(sprintf("[split_purify.R] wrote %s\n", out_unpurified_rds))

# ----- SPLIT::purify (Rmd lines 509-513) -----------------------------
# The Rmd:
#   res_split <- SPLIT::purify(
#     counts = GetAssayData(xe, assay="Xenium", layer="counts"),
#     rctd = rctd_result,
#     DO_purify_singlets = TRUE
#   )
# We swap the assay name (Xenium -> Proseg by default) via --assay-name.
cat("[split_purify.R] SPLIT::purify...\n")
res_split <- SPLIT::purify(
  counts             = GetAssayData(xe_unpurified, assay = assay_name, layer = "counts"),
  rctd               = rctd_result,
  DO_purify_singlets = DO_purify_singlets
)

# ----- Build purified Seurat (Rmd lines 517-522) ---------------------
# The Rmd:
#   xe_purified <- CreateSeuratObject(
#     counts    = res_split$purified_counts,
#     meta.data = res_split$cell_meta,
#     assay     = "Xenium"
#   )
xe_purified <- CreateSeuratObject(
  counts    = res_split$purified_counts,
  meta.data = res_split$cell_meta,
  assay     = assay_name
)

# Back-copy pre-purification metadata from xe_unpurified onto the
# purified object (proseg originals, xenium-preprocess Step 1
# enrichments, x/y, any caller-attached columns). RCTD-added columns
# and Seurat-managed nCount_*/nFeature_*/orig.ident are excluded; the
# helper is strictly additive (never overwrites columns SPLIT::purify
# or CreateSeuratObject already put on xe_purified). See
# _preserve_meta.R for the exclusion rules.
xe_purified <- preserve_meta_from_unpurified(
  pur_seurat      = xe_purified,
  unp_seurat      = xe_unpurified,
  rctd_added_cols = colnames(results_df)
)

cat(sprintf("[split_purify.R]   purified variant: cells=%d, features=%d\n",
            ncol(xe_purified), nrow(xe_purified)))

# ----- Save purified.rds --------------------------------------------
dir.create(dirname(out_purified_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(xe_purified, file = out_purified_rds)
cat(sprintf("[split_purify.R] wrote %s\n", out_purified_rds))
