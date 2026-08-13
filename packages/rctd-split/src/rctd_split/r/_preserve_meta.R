# ---------------------------------------------------------------------
# _preserve_meta.R -- helper for split_purify.R (Stage 2).
#
# Additively back-copies pre-purification meta.data columns from
# xe_unpurified onto xe_purified. Fix for (internal issue review): before
# this helper, only x/y were carried over; now the full pre-SPLIT
# metadata (proseg originals + xenium-preprocess Step 1 enrichments +
# any caller-attached columns) survives on the purified variant.
#
# Deliberately excluded:
#   - RCTD-added columns (already on xe_purified via SPLIT::purify's
#     cell_meta; copying again would either error on duplicate or
#     silently overwrite with the same value).
#   - Seurat-managed per-assay columns matching ^n(Count|Feature)_
#     plus orig.ident (recomputed by CreateSeuratObject against the
#     purified counts; copying from unpurified would corrupt them).
#   - Any column already present on xe_purified (strictly additive —
#     never overwrite what SPLIT::purify or CreateSeuratObject
#     produced for the purified cells).
#
# Covered by tests/test_preserve_meta_R.py (which drives an Rscript
# invocation of test_preserve_meta.R).
# ---------------------------------------------------------------------

preserve_meta_from_unpurified <- function(pur_seurat, unp_seurat,
                                          rctd_added_cols) {
  seurat_managed_cols <- c(
    "orig.ident",
    grep("^n(Count|Feature)_", colnames(unp_seurat@meta.data), value = TRUE)
  )
  already_on_purified <- colnames(pur_seurat@meta.data)
  preserve_cols <- setdiff(
    colnames(unp_seurat@meta.data),
    c(rctd_added_cols, seurat_managed_cols, already_on_purified)
  )
  purified_cells <- colnames(pur_seurat)
  carryover_cells <- intersect(purified_cells, colnames(unp_seurat))
  if (length(preserve_cols) == 0 || length(carryover_cells) == 0) {
    cat(sprintf(
      "[preserve_meta] nothing to copy (candidate cols=%d, carryover cells=%d)\n",
      length(preserve_cols), length(carryover_cells)
    ))
    return(pur_seurat)
  }
  extra_meta <- unp_seurat@meta.data[carryover_cells, preserve_cols,
                                     drop = FALSE]
  pur_seurat <- Seurat::AddMetaData(pur_seurat, extra_meta)
  cat(sprintf(
    "[preserve_meta] preserved %d meta.data column(s) on purified (%d/%d cells): %s\n",
    length(preserve_cols), length(carryover_cells), length(purified_cells),
    paste(preserve_cols, collapse = ", ")
  ))
  pur_seurat
}
