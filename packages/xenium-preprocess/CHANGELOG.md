# Changelog

All notable changes to xenium-preprocess will be documented here. Follows
[Keep a Changelog](https://keepachangelog.com/) shape.

## [Unreleased]

### Changed (2026-08-18)
- Top-level driver-side rename: `scripts/submit_step1.sbatch` →
  `scripts/submit_xenium-preprocess.sbatch`. The Python package
  surface is unaffected; only affects direct-sbatch callers (the
  documented flow is `submit_workflow.sh`). See top-level
  `CHANGELOG.md` §0.2.0 for the sibling renames.

### Changed (2026-08-12, the user clarification (internal issue review)
- `rctd_prep` RESTORED to `DEFAULT_STAGES`. The 2026-08-11 change
  above dropped `split_prep` + `rctd_prep` together, reading the caller's
  "remove all 10x bundle" ask as targeting the downstream RDS as
  well. Internal review clarifies she still needs
  `<sample>_test_object.rds` — only the persisted
  mtx/features/barcodes bundle under `split_prep/` was meant to
  disappear. `split_prep` stays out of defaults.
- `pipeline.run` now auto-produces the `split_prep/` mtx bundle as a
  TRANSIENT prerequisite for `rctd_prep` when the user did not
  explicitly request `split_prep` in `--stages`, and deletes
  `<run_dir>/split_prep/` after `rctd_prep` writes the RDS. Users
  who want to keep the bundle on disk can add `split_prep` to
  `--stages`; the earlier explicit block still writes it there and
  the transient-cleanup skips (opt-in is honoured).
- CLI `--stages` help text updated to describe the new behaviour:
  only `preprocess` and `split_prep` are omitted by default; the
  RCTD RDS is produced.
- `tests/test_config.py::test_default_stages_excludes_legacy_10x_bundle`
  replaced by `test_default_stages_split_out_rctd_in` +
  `test_default_stages_rctd_runs_after_prereqs`, pinning the
  new invariant.

### Changed (2026-08-11, user request (internal issue review)
- `enrich_xenium_id` direction REVERSED. Previously: for each proseg
  cell → find nearest xenium cell → write `xenium_cell_id_nn` +
  `xenium_id_match` (plus distance/note columns) onto
  `<sample>_proseg_raw.h5ad`. Now: for each xenium cell → find
  nearest proseg cell → write `proseg_cell_id_nn` +
  `proseg_id_nn_distance` + `proseg_id_nn_note` onto
  `<sample>_xenium_ranger.h5ad`. The proseg h5ad is READ-ONLY under
  the new stage. The `original_id_col` / `match_flag_col` params
  and the "audit peer" comparison are dropped (there is no
  original proseg id on the xenium side to compare against).
- CLI flag `--xenium-cells` removed — the xenium cells file is
  no longer an input to `enrich_xenium_id`; the stage reads
  centroids directly from `<sample>_xenium_ranger.h5ad`.
- Stage dispatch order updated so `enrich_xenium_id` runs AFTER
  `xenium_ranger_to_anndata` (needs both source h5ads on disk).
- `DEFAULT_STAGES` now: `proseg_to_anndata → qc_filter →
  xenium_ranger_to_anndata → enrich_xenium_id`. `split_prep` and
  `rctd_prep` dropped from defaults — they still live in
  `VALID_STAGES` for opt-in via `--stages`, but the 10x-style
  bundle + downstream RDS are legacy artifacts of the pre-reorg
  architecture and no longer emitted by default (user request — "remove
  all 10x bundle" output).
- `config/default.yaml` `enrich_xenium_id` section rewritten to
  reflect the new column names; obsolete
  `xenium_cells_path` / `xenium_id_col` / `xenium_x_col` /
  `xenium_y_col` / `original_id_col` / `match_flag_col` keys
  removed.

### Added
- Dual-matrix support (user request 2026-07-10, (internal issue review):
  - `proseg_to_anndata` reads BOTH the `expected-counts` and
    `maxpost_counts` proseg outputs and stores them as
    `adata.layers['expected_counts']` (mirrored to `.X`) and
    `adata.layers['maxpost_counts']`. Validates shape + gene-column
    order across both matrices; fails loud on mismatch.
    New config: `proseg_to_anndata.maxpost_matrix_glob` (default
    `"maxpost_counts*"`; set null to disable). New CLI flag:
    `--maxpost-matrix`.
  - `preprocess` runs TWO independent passes when
    `preprocess.dual_matrix_mode: true` (default) — one on the
    `expected_counts` layer, one on the `maxpost_counts` layer.
    PCA/KNN/UMAP/Leiden and (optional) cell-type annotation are
    layer-scoped via `_expected` / `_maxpost` suffixes; the two
    label sets may disagree, and that divergence is inspectable.
    New CLI flag: `--dual-matrix-mode`.
  - `split_prep.layer` default changed from `counts` to
    `maxpost_counts` (the integer proseg output that RCTD needs).
    Threaded through by `pipeline.py` from the NEW config
    `rctd_prep.source_layer: maxpost_counts` — one truth-source
    for which layer feeds the RCTD test object.
- `tests/test_dual_matrix.py` — new tests covering: both layers
  populated on a synthetic proseg tree; shape-mismatch failure mode;
  `maxpost_matrix_glob=None` back-compat; dual-mode + single-mode
  preprocess round-trips.

### Changed
- `preprocess.daniel_approach_hvg_pca` gains `source_layer` +
  `suffix` params so two dual-matrix passes coexist without
  overwriting each other's `X_pca` / `X_umap` / `X_clipped` /
  KNN / PC-loadings / variance-ratio outputs.
- `preprocess.leiden_clustering_via_knn_graph` gains a `suffix`
  param aligned with the same pattern (was previously suffix-aware
  but the calling driver never used it).
- `preprocess.calculate_module_score`, `find_potential_cell_types`,
  `detect_low_count`, `map_celltypes` now take an explicit
  `leiden_key` (rather than computing `f"leiden_{resolution}"`
  internally), so the dual-mode driver can point each pass at its
  own `leiden_{res}_expected` / `leiden_{res}_maxpost` partition.
- `docs/methods.md`, `docs/usage.md`, and `README.md` updated to
  document the dual-matrix outputs and the `rctd_prep.source_layer`
  truth-source.

### Backward compat
- Setting `preprocess.dual_matrix_mode: false` AND
  `proseg_to_anndata.maxpost_matrix_glob: null` reverts to exactly
  the pre-2026-07-10 single-matrix, single-pass pipeline (unsuffixed
  `X_pca` / `X_umap` / `leiden_{res}` / `global_level1_celltype`
  keys, `umap_leiden.png` / `umap_global_level1_celltype.png` /
  `umap_tumor_markers_<TumorType>.png` plot names).

## [0.1.0] — 2026-07-09

### Added
- Initial package scaffold (xenium-preprocess of the internal SPLIT/Proseg workflow spatial-data
  preprocessing pipeline).
- Four stage modules under `xenium_preprocess.stages`
  (`proseg_to_anndata`, `preprocess`, `split_prep`, `rctd_prep`),
  each a verbatim port of one of the user's reference scripts /
  notebooks.
- `xenium_preprocess/r/rctd_prep.R` — the R script the `rctd_prep`
  stage shells out to. Builds a Seurat spatial object with the
  Proseg assay from the SPLIT mtx/features/barcodes triple + the
  metadata + spatial-coords sidecars, saves as `test_object.rds`.
- Config machinery split into `xenium_preprocess.config` (YAML load
  + deep merge + validation).
- CLI split into `xenium_preprocess.cli` (argparse) +
  `xenium_preprocess.pipeline` (stage-orchestration loop with
  sentinel-file resume and `--stages` filter).
- `pyproject.toml` + `environment.yml` describing the `xeniumPreprocess`
  conda env.
- `scripts/submit.slurm.sh` — sbatch wrapper carrying the same
  Slurm-tee/pipefail deadlock fix that landed in the H&E pipeline.
- `CITATION.cff`, MIT `LICENSE`, `docs/` markdown pages, four smoke
  tests under `tests/`.
