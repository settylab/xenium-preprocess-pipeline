#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# rctd_run.R -- Stage 1 of the rctd-split pipeline (Stage D of the
# broader SPLIT / RCTD workflow).
#
# Reads the spatial test object RDS (from step 1 xenium-preprocess's
# rctd_prep stage) and the scRNA reference RDS (from step 3 ref-build's
# rctd_reference_build stage), runs create.RCTD + run.RCTD, saves the
# result.
#
# Ported verbatim from an internal SPLIT/Proseg workflow
#   the internal SPLIT/Proseg workflow
# lines 305-363 -- every parameter below traces to a line in the Rmd.
#
# Invocation (via rctd_split.stages.rctd_run):
#   Rscript --vanilla rctd_run.R \
#       --test-object=<path> \
#       --reference-rds=<path> \
#       --out-rds=<path> \
#       --assay-name=Proseg \
#       --spatial-reduction=spatial \
#       --umi-min=10 \
#       --counts-min=10 \
#       --umi-min-sigma=100 \
#       --max-cores=4 \
#       --cell-min-instance=20 \
#       --doublet-mode=doublet \
#       --label-slash-replacement=_
# ---------------------------------------------------------------------

# ----- Prepend the user-local R lib so spacexr is found ---------------
# The default fhR/4.4.1-foss-2023b module carries Seurat + Matrix but
# NOT spacexr. See docs/install.md for the one-time user-local install
# at ~/.claude/r_libs/4.4.1. Adding it to .libPaths() here is
# idempotent: if the dir doesn't exist, we skip it.
user_lib <- file.path(Sys.getenv("HOME"), ".claude", "r_libs", "4.4.1")
if (dir.exists(user_lib)) {
  .libPaths(c(user_lib, .libPaths()))
  cat(sprintf("[rctd_run.R] prepended user lib: %s\n", user_lib))
}

# ----- Propagate .libPaths() to R_LIBS_USER for PSOCK workers ---------
# spacexr's process_beads_batch / decompose_batch call
# parallel::makeCluster(numCores) which spawns fresh Rscript workers by
# fork+exec — PSOCK workers inherit the parent's environment vars, and
# each worker's .libPaths() is derived at startup from R_LIBS_USER. The
# parent-only .libPaths() mutation above is in-memory ONLY; it does NOT
# reach the workers. Result: workers try to loadNamespace("spacexr")
# during closure deserialisation and die with
#   "there is no package called 'spacexr'"
# whenever spacexr lives on a path (user_lib, an renv library, ...) that
# isn't in the R_LIBS_USER env value the parent inherited from its
# caller. We fix this by writing the parent's full .libPaths() back into
# R_LIBS_USER so every PSOCK worker sees the same library search path.
# Pre-registering our own cluster with clusterCall(cl, .libPaths, ...)
# would NOT work: spacexr always calls makeCluster() itself and its
# doParallel::registerDoParallel(cl) call overwrites any pre-registered
# backend.
Sys.setenv(R_LIBS_USER = paste(.libPaths(), collapse = .Platform$path.sep))
cat(sprintf("[rctd_run.R] R_LIBS_USER set to parent .libPaths() (%d entries) for PSOCK worker inheritance\n",
            length(.libPaths())))

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
  "test-object", "reference-rds", "out-rds",
  "assay-name", "spatial-reduction",
  "umi-min", "counts-min", "umi-min-sigma",
  "max-cores", "cell-min-instance",
  "doublet-mode", "label-slash-replacement"
)
missing <- setdiff(need, names(kv))
if (length(missing)) {
  stop("[rctd_run.R] missing required args: ",
       paste(missing, collapse = ", "))
}

test_object_path  <- kv[["test-object"]]
reference_rds     <- kv[["reference-rds"]]
out_rds           <- kv[["out-rds"]]
assay_name        <- kv[["assay-name"]]
spatial_reduction <- kv[["spatial-reduction"]]
UMI_min           <- as.integer(kv[["umi-min"]])
counts_MIN        <- as.integer(kv[["counts-min"]])
UMI_min_sigma     <- as.integer(kv[["umi-min-sigma"]])
max_cores         <- as.integer(kv[["max-cores"]])
CELL_MIN_INSTANCE <- as.integer(kv[["cell-min-instance"]])
doublet_mode      <- kv[["doublet-mode"]]
slash_repl        <- kv[["label-slash-replacement"]]

cat(sprintf("[rctd_run.R] test_object       = %s\n", test_object_path))
cat(sprintf("[rctd_run.R] reference_rds     = %s\n", reference_rds))
cat(sprintf("[rctd_run.R] out_rds           = %s\n", out_rds))
cat(sprintf("[rctd_run.R] assay_name        = %s\n", assay_name))
cat(sprintf("[rctd_run.R] spatial_reduction = %s\n", spatial_reduction))
cat(sprintf("[rctd_run.R] UMI_min           = %d\n", UMI_min))
cat(sprintf("[rctd_run.R] counts_MIN        = %d\n", counts_MIN))
cat(sprintf("[rctd_run.R] UMI_min_sigma     = %d\n", UMI_min_sigma))
cat(sprintf("[rctd_run.R] max_cores         = %d\n", max_cores))
cat(sprintf("[rctd_run.R] CELL_MIN_INSTANCE = %d\n", CELL_MIN_INSTANCE))
cat(sprintf("[rctd_run.R] doublet_mode      = %s\n", doublet_mode))
cat(sprintf("[rctd_run.R] slash_repl        = %s\n", slash_repl))

# ----- Load inputs ---------------------------------------------------
cat(sprintf("[rctd_run.R] reading test object = %s\n", test_object_path))
xe <- readRDS(test_object_path)
cat(sprintf("[rctd_run.R]   test object: class=%s, cells=%d, features=%d\n",
            class(xe)[1], ncol(xe), nrow(xe)))

cat(sprintf("[rctd_run.R] reading reference   = %s\n", reference_rds))
ref.obj <- readRDS(reference_rds)
cat(sprintf("[rctd_run.R]   reference: class=%s, cells=%d, features=%d\n",
            class(ref.obj)[1],
            length(ref.obj@cell_types),
            nrow(ref.obj@counts)))

# ----- Intersect gene panels (Rmd line 311) --------------------------
common_genes <- intersect(rownames(xe), rownames(ref.obj@counts))
if (length(common_genes) == 0) {
  stop("[rctd_run.R] zero common genes between spatial (", nrow(xe),
       ") and reference (", nrow(ref.obj@counts),
       "); check row-name conventions / gene symbols.")
}
cat(sprintf("[rctd_run.R]   common genes: %d\n", length(common_genes)))

# ----- Build SpatialRNA test.obj (Rmd lines 324-327) -----------------
# The Rmd pulls coords from the "spatial" reduction's cell.embeddings.
# Step 1's rctd_prep sets this reduction with columns ST_1, ST_2; we
# preserve that but only take the first two dims for the coords frame.
if (!(spatial_reduction %in% names(xe@reductions))) {
  stop("[rctd_run.R] test object has no reduction named '", spatial_reduction,
       "'; available: ",
       paste(names(xe@reductions), collapse = ", "))
}
coords <- as.data.frame(xe@reductions[[spatial_reduction]]@cell.embeddings)
if (ncol(coords) < 2) {
  stop("[rctd_run.R] spatial reduction '", spatial_reduction,
       "' has < 2 columns; RCTD expects (x, y).")
}
coords <- coords[, 1:2, drop = FALSE]

# Sanitise reference cell-type labels — spacexr factor levels reject '/'
# (matches the Rmd's gsub("/", "_", ...) on line 306). The step-3
# reference already applies this on build, but the guard is cheap and
# lets this stage stay standalone.
ref_labels_raw <- as.character(ref.obj@cell_types)
ref_labels <- gsub("/", slash_repl, ref_labels_raw, fixed = TRUE)
if (any(ref_labels != ref_labels_raw)) {
  n_changed <- length(unique(ref_labels_raw[ref_labels != ref_labels_raw]))
  cat(sprintf("[rctd_run.R]   sanitised %d unique reference labels containing '/'\n",
              n_changed))
  cell_names <- names(ref.obj@cell_types)
  if (is.null(cell_names)) cell_names <- colnames(ref.obj@counts)
  new_ct <- factor(ref_labels)
  names(new_ct) <- cell_names
  ref.obj@cell_types <- new_ct
}

test.obj <- SpatialRNA(
  coords      = coords,
  counts      = GetAssayData(xe, assay = assay_name, layer = "counts")[common_genes, ],
  require_int = TRUE
)
cat(sprintf("[rctd_run.R]   test.obj built: %d spots, %d genes\n",
            ncol(test.obj@counts), nrow(test.obj@counts)))

# ----- Restrict reference to the common genes -----------------------
ref.obj@counts <- ref.obj@counts[common_genes, , drop = FALSE]

# ----- create.RCTD (Rmd lines 344-353) ------------------------------
# Every keyword arg below is verbatim from the Rmd; class_df=NULL
# mirrors the Rmd's `class_df = class_df` where `class_df` is set to
# NULL on line 330-331.
cat("[rctd_run.R] calling create.RCTD...\n")
RCTD <- create.RCTD(
  test.obj,
  ref.obj,
  UMI_min           = UMI_min,
  counts_MIN        = counts_MIN,
  UMI_min_sigma     = UMI_min_sigma,
  max_cores         = max_cores,
  CELL_MIN_INSTANCE = CELL_MIN_INSTANCE,
  class_df          = NULL
)

# ----- run.RCTD (Rmd line 362) --------------------------------------
cat(sprintf("[rctd_run.R] calling run.RCTD (doublet_mode = %s)...\n", doublet_mode))
RCTD <- run.RCTD(RCTD, doublet_mode = doublet_mode)

# ----- Save ----------------------------------------------------------
dir.create(dirname(out_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(RCTD, file = out_rds)

n_results <- if (!is.null(RCTD@results$results_df))
  nrow(RCTD@results$results_df) else NA_integer_
cat(sprintf("[rctd_run.R] wrote %s (results_df rows: %s)\n",
            out_rds,
            ifelse(is.na(n_results), "NA", as.character(n_results))))
