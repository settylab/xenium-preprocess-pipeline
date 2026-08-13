#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# verify_parallel_fix.R -- live end-to-end verification that the
# libPaths-propagation fix in ../src/rctd_split/r/rctd_run.R lets
# spacexr's parallel PSOCK workers find `spacexr` when it's only on a
# user-local .libPaths() entry that is NOT already on R_LIBS_USER.
#
# Reproduces Tracy's failure mode:
#   * caller sets R_LIBS_USER to a path that does NOT contain spacexr
#     (Tracy's SPLIT renv path).
#   * rctd_run.R prepends ~/.claude/r_libs/4.4.1 (where spacexr lives)
#     to .libPaths() but WITHOUT the fix, doesn't propagate to
#     R_LIBS_USER, so PSOCK workers spawned by spacexr's
#     process_beads_batch %dopar% loop can't loadNamespace("spacexr").
# ---------------------------------------------------------------------

# --- The fix (mirror of the block added to r/rctd_run.R) --------------
user_lib <- file.path(Sys.getenv("HOME"), ".claude", "r_libs", "4.4.1")
if (dir.exists(user_lib)) {
  .libPaths(c(user_lib, .libPaths()))
  cat(sprintf("[verify] prepended user lib: %s\n", user_lib))
}
Sys.setenv(R_LIBS_USER = paste(.libPaths(), collapse = .Platform$path.sep))
cat(sprintf("[verify] R_LIBS_USER set to parent .libPaths() (%d entries) for PSOCK inheritance\n",
            length(.libPaths())))

suppressPackageStartupMessages({
  library(spacexr)
  library(Matrix)
})

# --- Load the Merfish example data (small, ships with spacexr) --------
ext <- system.file("extdata", package = "spacexr")

ref_counts <- read.csv(file.path(ext, "Reference/Merfish_Ref/counts.csv"),
                       row.names = 1, check.names = FALSE)
ref_ct <- read.csv(file.path(ext, "Reference/Merfish_Ref/cell_types.csv"),
                   row.names = 1)
ref_numi <- read.csv(file.path(ext, "Reference/Merfish_Ref/nUMI.csv"),
                     row.names = 1)
ct_vec <- setNames(as.factor(ref_ct[[1]]), rownames(ref_ct))
numi_vec <- setNames(as.numeric(ref_numi[[1]]), rownames(ref_numi))
ref <- Reference(as.matrix(ref_counts), ct_vec, numi_vec)

sp_counts <- read.csv(file.path(ext, "SpatialRNA/MerfishVignette/counts.csv"),
                      row.names = 1, check.names = FALSE)
sp_coords <- read.csv(file.path(ext, "SpatialRNA/MerfishVignette/coords.csv"),
                      row.names = 1)
sp_numi <- read.csv(file.path(ext, "SpatialRNA/MerfishVignette/nUMI.csv"),
                    row.names = 1)
sp_numi_vec <- setNames(as.numeric(sp_numi[[1]]), rownames(sp_numi))
puck <- SpatialRNA(sp_coords, as.matrix(sp_counts), sp_numi_vec)

cat(sprintf("[verify] reference: %d cells, %d genes\n", ncol(ref@counts), nrow(ref@counts)))
cat(sprintf("[verify] spatial:   %d spots, %d genes\n", ncol(puck@counts), nrow(puck@counts)))

# --- Full RCTD run with max_cores=4 -----------------------------------
# CELL_MIN_INSTANCE=1 because a couple of Merfish cell types have < 25
# cells; that's not what we're testing, so relax it.
cat("[verify] create.RCTD (max_cores=4)...\n")
myRCTD <- create.RCTD(puck, ref,
                     max_cores = 4,
                     CELL_MIN_INSTANCE = 1,
                     UMI_min = 10,
                     counts_MIN = 10,
                     UMI_min_sigma = 100)

cat("[verify] run.RCTD (doublet_mode = doublet)...\n")
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")

n <- if (!is.null(myRCTD@results$results_df)) nrow(myRCTD@results$results_df) else NA
cat(sprintf("[verify] SUCCESS: results_df rows = %s\n",
            ifelse(is.na(n), "NA", as.character(n))))
