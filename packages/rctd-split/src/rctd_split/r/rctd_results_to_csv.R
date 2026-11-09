#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# rctd_results_to_csv.R -- v0.2.9
#
# Extract per-cell RCTD results from a spacexr RCTD .rds (+ SPLIT's
# post-process step) into a flat CSV, ready for Python-side obs joins
# onto `xenium_ranger.h5ad`.
#
# Mirrors the first few lines of `split_purify.R`:
#   rctd_result <- readRDS(rctd_results_rds)
#   rctd_result <- SPLIT::run_post_process_RCTD(rctd_result)
#   results_df  <- rctd_result@results$results_df
#
# The post_process step adds the SPLIT-emitted columns
# (``purification_status``, ``w1_larger_w2``, ``same_class``) that the
# bare RCTD object doesn't carry. Matches exactly what `split_purify.R`
# consumes, so the CSV is a 1:1 Python-readable view of what the
# pipeline already uses on the R side.
#
# Invocation:
#   Rscript --vanilla rctd_results_to_csv.R \
#       --rctd-results-rds=<path>  \
#       --out-csv=<path>          \
#       [--r-lib-paths=<lib1:lib2:...>]
#
# Output CSV columns (first column is `cell_id` — the RCTD rowname):
#   cell_id, spot_class, first_type, second_type, first_class,
#   second_class, min_score, singlet_score, conv_all, conv_doublet,
#   purification_status, w1_larger_w2, same_class
#
# Columns that aren't present in a particular rctd_result's results_df
# (e.g. older RCTD runs) are silently omitted from the CSV. Missing
# values are serialised as empty strings for categoricals and as `NA`
# for numerics.
# ---------------------------------------------------------------------

options(warn = 1, stringsAsFactors = FALSE)

# --- parse --key=value args -----------------------------------------
parse_kv <- function(args) {
  kv <- list()
  for (a in args) {
    m <- regmatches(a, regexec("^--([^=]+)=(.*)$", a))[[1]]
    if (length(m) == 3L) kv[[m[2]]] <- m[3]
  }
  kv
}

args <- commandArgs(trailingOnly = TRUE)
kv   <- parse_kv(args)

rctd_results_rds <- kv[["rctd-results-rds"]]
out_csv          <- kv[["out-csv"]]
r_lib_paths      <- kv[["r-lib-paths"]]

for (nm in c("rctd-results-rds", "out-csv")) {
  if (is.null(kv[[nm]]) || !nzchar(kv[[nm]]))
    stop(sprintf("missing required --%s=<path>", nm))
}

if (!is.null(r_lib_paths) && nzchar(r_lib_paths)) {
  paths <- strsplit(r_lib_paths, ":", fixed = TRUE)[[1]]
  paths <- paths[nzchar(paths)]
  .libPaths(c(paths, .libPaths()))
}

cat(sprintf("[rctd_results_to_csv.R] rctd_results_rds = %s\n", rctd_results_rds))
cat(sprintf("[rctd_results_to_csv.R] out_csv          = %s\n", out_csv))

suppressPackageStartupMessages({
  library(spacexr)
  library(SPLIT)
})

# --- load + post-process --------------------------------------------
cat("[rctd_results_to_csv.R] reading RCTD results...\n")
rctd_result <- readRDS(rctd_results_rds)

cat("[rctd_results_to_csv.R] SPLIT::run_post_process_RCTD...\n")
rctd_result <- SPLIT::run_post_process_RCTD(rctd_result)

results_df <- rctd_result@results$results_df
cat(sprintf("[rctd_results_to_csv.R]   results_df: %d rows x %d cols\n",
            nrow(results_df), ncol(results_df)))
cat(sprintf("[rctd_results_to_csv.R]   columns: %s\n",
            paste(colnames(results_df), collapse = ", ")))

# --- select the columns the Python side joins onto obs ---------------
# Any column requested here that's absent from results_df is dropped
# silently — supports older RCTD results that may not carry every
# SPLIT-emitted column.
want_cols <- c(
  "spot_class", "first_type", "second_type",
  "first_class", "second_class",
  "min_score", "singlet_score",
  "conv_all", "conv_doublet",
  "purification_status", "w1_larger_w2", "same_class"
)
have_cols <- intersect(want_cols, colnames(results_df))
if (length(have_cols) == 0L)
  stop("rctd_results.rds has none of the expected results_df columns.")

cat(sprintf("[rctd_results_to_csv.R]   emitting cols: %s\n",
            paste(have_cols, collapse = ", ")))

out <- data.frame(
  cell_id = rownames(results_df),
  results_df[, have_cols, drop = FALSE],
  check.names = FALSE,
  row.names   = NULL
)

# --- write CSV (atomic via tmp + file.rename) ------------------------
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
tmp <- paste0(out_csv, ".tmp")
write.csv(out, tmp, row.names = FALSE, quote = TRUE, na = "")
file.rename(tmp, out_csv)

cat(sprintf("[rctd_results_to_csv.R] wrote %d rows x %d cols -> %s\n",
            nrow(out), ncol(out) - 1L, out_csv))
