#!/usr/bin/env Rscript
# ---------------------------------------------------------------------
# test_preserve_meta.R -- unit test for _preserve_meta.R.
#
# Fabricates a pair of Seurat objects (unpurified: 100 cells with
# proseg-origin + Step 1 enrichment + x/y + RCTD-added metadata;
# purified: 50-cell subset with only SPLIT::purify-shaped cell_meta)
# and asserts that preserve_meta_from_unpurified():
#
#   1. Copies every pre-purification column to the purified object —
#      proseg originals (cell_area, my_custom_col), Step 1 enrichments
#      (original_cell_id, xenium_cell_id_nn, xenium_id_nn_distance,
#      xenium_id_nn_note, xenium_id_match), and x/y coordinates.
#   2. Preserves the correct per-cell VALUES (indexed to the purified
#      cell subset, not misaligned).
#   3. Does NOT overwrite RCTD-added columns already on the purified
#      object from SPLIT::purify's cell_meta.
#   4. Does NOT overwrite the Seurat-managed nCount_*/nFeature_*/
#      orig.ident columns (those must reflect the PURIFIED counts,
#      not the unpurified ones).
#
# Exits 0 on success; stopifnot() fires with exit 1 on any failure.
# Driven by tests/test_preserve_meta_R.py (which handles the Rscript
# launch, module load, and pytest skip when Rscript is unavailable).
# ---------------------------------------------------------------------

user_lib <- file.path(Sys.getenv("HOME"), ".claude", "r_libs", "4.4.1")
if (dir.exists(user_lib)) {
  .libPaths(c(user_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# Locate this script -> package r/ dir; source the helper under test.
this_file <- sub("^--file=", "",
                 grep("^--file=", commandArgs(trailingOnly = FALSE),
                      value = TRUE))
tests_dir <- dirname(normalizePath(this_file))
r_dir <- normalizePath(file.path(tests_dir, "..", "src", "rctd_split", "r"))
source(file.path(r_dir, "_preserve_meta.R"))

# ----- Fabricate xe_unpurified --------------------------------------
set.seed(1234)
n_cells <- 100
n_genes <- 50
counts_unp <- matrix(rpois(n_genes * n_cells, lambda = 2),
                     nrow = n_genes, ncol = n_cells)
rownames(counts_unp) <- paste0("gene", seq_len(n_genes))
colnames(counts_unp) <- paste0("cell", seq_len(n_cells))

meta_unp <- data.frame(
  # Proseg-origin columns (Tracy's ask):
  cell_area          = runif(n_cells, 10, 500),
  my_custom_col      = seq_len(n_cells),
  # Xenium-preprocess Step 1 enrichment columns:
  original_cell_id   = paste0("orig_", seq_len(n_cells)),
  xenium_cell_id_nn  = paste0("xen_", seq_len(n_cells)),
  xenium_id_nn_distance = runif(n_cells, 0, 10),
  xenium_id_nn_note  = sample(c("exact", "nearest", "unmatched"),
                              n_cells, replace = TRUE),
  xenium_id_match    = sample(c(TRUE, FALSE), n_cells, replace = TRUE),
  # Spatial coords:
  x = runif(n_cells, 0, 1000),
  y = runif(n_cells, 0, 1000),
  # RCTD-added columns (these come from AddMetaData(xe, results_df) in
  # the real pipeline; simulate them so the exclusion rule can be
  # exercised):
  first_type         = sample(c("tumor", "stroma", "immune"),
                              n_cells, replace = TRUE),
  second_type        = sample(c("tumor", "stroma", "immune"),
                              n_cells, replace = TRUE),
  spot_class         = sample(c("singlet", "doublet_certain"),
                              n_cells, replace = TRUE),
  purification_status = sample(c("singlet", "doublet"),
                               n_cells, replace = TRUE),
  w1_larger_w2       = sample(c(TRUE, FALSE), n_cells, replace = TRUE),
  same_class         = sample(c(TRUE, FALSE), n_cells, replace = TRUE),
  stringsAsFactors   = FALSE,
  row.names          = colnames(counts_unp)
)
xe_unpurified <- CreateSeuratObject(counts    = counts_unp,
                                    meta.data = meta_unp,
                                    assay     = "Proseg")

# The simulated RCTD-added columns (mirrors real results_df cols):
rctd_added_cols <- c("first_type", "second_type", "spot_class",
                     "purification_status", "w1_larger_w2", "same_class")

# ----- Fabricate xe_purified (subset of 50 cells) -------------------
# Emulate what CreateSeuratObject(..., meta.data = res_split$cell_meta)
# yields: only the RCTD-shaped columns from SPLIT::purify's cell_meta
# land on purified; proseg / enrichment / x/y columns do NOT.
pur_cells <- colnames(counts_unp)[1:50]
counts_pur <- counts_unp[, pur_cells]
cell_meta_pur <- meta_unp[pur_cells, rctd_added_cols, drop = FALSE]
xe_purified <- CreateSeuratObject(counts    = counts_pur,
                                  meta.data = cell_meta_pur,
                                  assay     = "Proseg")

# Sanity: verify the fixture matches the real pipeline shape.
pre_cols <- colnames(xe_purified@meta.data)
stopifnot("orig.ident" %in% pre_cols)
stopifnot("nCount_Proseg" %in% pre_cols)
stopifnot("nFeature_Proseg" %in% pre_cols)
stopifnot(all(rctd_added_cols %in% pre_cols))
# The proseg/enrichment cols must NOT be on purified yet — that's the
# bug this helper fixes.
stopifnot(!("cell_area" %in% pre_cols))
stopifnot(!("original_cell_id" %in% pre_cols))
stopifnot(!("x" %in% pre_cols))

# ----- Run the helper -----------------------------------------------
xe_purified <- preserve_meta_from_unpurified(
  pur_seurat      = xe_purified,
  unp_seurat      = xe_unpurified,
  rctd_added_cols = rctd_added_cols
)

# ----- Assertion 1: every expected column now on purified -----------
expect_now_present <- c(
  "cell_area", "my_custom_col",
  "original_cell_id", "xenium_cell_id_nn", "xenium_id_nn_distance",
  "xenium_id_nn_note", "xenium_id_match",
  "x", "y"
)
post_cols <- colnames(xe_purified@meta.data)
missing_cols <- setdiff(expect_now_present, post_cols)
if (length(missing_cols) > 0) {
  stop("[test_preserve_meta] expected columns not preserved: ",
       paste(missing_cols, collapse = ", "))
}
cat(sprintf("[test_preserve_meta] OK: all %d expected columns preserved.\n",
            length(expect_now_present)))

# ----- Assertion 2: per-cell VALUES match unpurified ----------------
for (col in expect_now_present) {
  got  <- xe_purified@meta.data[pur_cells, col]
  want <- meta_unp[pur_cells, col]
  if (!isTRUE(all.equal(got, want))) {
    stop(sprintf("[test_preserve_meta] value mismatch on column %s:\n  got:  %s\n  want: %s",
                 col,
                 paste(head(as.character(got), 5), collapse = ","),
                 paste(head(as.character(want), 5), collapse = ",")))
  }
}
cat("[test_preserve_meta] OK: per-cell values match unpurified for all preserved columns.\n")

# ----- Assertion 3: RCTD-added columns not overwritten --------------
# xe_purified had cell_meta_pur values pre-helper; ensure they didn't
# get clobbered. (Fabricated: pre-helper values == meta_unp values on
# the purified cells since we built cell_meta_pur that way. So this is
# really testing "no double-column and no NA-overwrite.")
for (col in rctd_added_cols) {
  n_matches <- sum(colnames(xe_purified@meta.data) == col)
  if (n_matches != 1) {
    stop(sprintf("[test_preserve_meta] RCTD column %s appears %d times (expected 1)",
                 col, n_matches))
  }
  got  <- xe_purified@meta.data[pur_cells, col]
  want <- cell_meta_pur[pur_cells, col]
  if (!isTRUE(all.equal(got, want))) {
    stop(sprintf("[test_preserve_meta] RCTD column %s was mutated by preserve helper",
                 col))
  }
}
cat("[test_preserve_meta] OK: RCTD-added columns not duplicated or mutated.\n")

# ----- Assertion 4: Seurat-managed columns still reflect purified counts
# nCount_Proseg on the purified object must equal colSums of the
# PURIFIED counts (which are identical to unp counts here because we
# didn't drop any UMIs, but a wrong copy from unp would still be equal
# by coincidence — test that they equal the purified counts colSums,
# and also that copying didn't produce duplicate columns).
n_seurat_matches <- length(grep("^nCount_Proseg$", colnames(xe_purified@meta.data)))
stopifnot(n_seurat_matches == 1)
n_seurat_matches <- length(grep("^nFeature_Proseg$", colnames(xe_purified@meta.data)))
stopifnot(n_seurat_matches == 1)
n_seurat_matches <- length(grep("^orig\\.ident$", colnames(xe_purified@meta.data)))
stopifnot(n_seurat_matches == 1)

expected_ncount <- Matrix::colSums(counts_pur)
got_ncount <- xe_purified@meta.data[pur_cells, "nCount_Proseg"]
stopifnot(all(got_ncount == expected_ncount))
cat("[test_preserve_meta] OK: Seurat-managed columns intact.\n")

# ----- Assertion 5: no purified cells lost ---------------------------
stopifnot(ncol(xe_purified) == length(pur_cells))
stopifnot(all(colnames(xe_purified) == pur_cells))
cat("[test_preserve_meta] OK: purified cell set unchanged.\n")

# ----- Assertion 6: NA-cell handling --------------------------------
# If xe_purified contained a cell NOT present in xe_unpurified (should
# never happen in practice but test the boundary), the helper's
# intersect() must silently drop it rather than error, and no NA rows
# should propagate. Re-run with a fabricated extra cell.
extra_cell <- "cellNEW"
counts_pur2 <- cbind(counts_pur,
                     matrix(rpois(n_genes, 2), nrow = n_genes, ncol = 1,
                            dimnames = list(NULL, extra_cell)))
cell_meta_pur2 <- rbind(
  cell_meta_pur,
  data.frame(
    first_type = "tumor", second_type = "tumor", spot_class = "singlet",
    purification_status = "singlet", w1_larger_w2 = TRUE, same_class = TRUE,
    stringsAsFactors = FALSE, row.names = extra_cell
  )
)
xe_pur2 <- CreateSeuratObject(counts    = counts_pur2,
                              meta.data = cell_meta_pur2,
                              assay     = "Proseg")
xe_pur2 <- preserve_meta_from_unpurified(xe_pur2, xe_unpurified,
                                         rctd_added_cols = rctd_added_cols)
stopifnot("cell_area" %in% colnames(xe_pur2@meta.data))
# The 50 known cells keep their real values:
stopifnot(all(xe_pur2@meta.data[pur_cells, "cell_area"] ==
              meta_unp[pur_cells, "cell_area"]))
# The extra cell gets NA (not in unp):
stopifnot(is.na(xe_pur2@meta.data[extra_cell, "cell_area"]))
cat("[test_preserve_meta] OK: cells not present in unpurified get NA (no error).\n")

cat("[test_preserve_meta] ALL ASSERTIONS PASSED.\n")
