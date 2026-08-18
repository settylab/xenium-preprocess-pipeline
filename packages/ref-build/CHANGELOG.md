# Changelog

All notable changes to ref-build will be documented here. Follows
[Keep a Changelog](https://keepachangelog.com/) shape.

## [0.1.0] — 2026-07-09

### Added
- Initial package scaffold (ref-build of the internal SPLIT/Proseg workflow spatial-data
  preprocessing pipeline).
- Five stage modules under `ref_build.stages`
  (`load_primary_and_donors`, `census`, `assemble`, `export_mtx`,
  `rctd_reference_build`), ported from the internal SPLIT/Proseg workflow's
  the reference-build notebook notebook family
  under the internal SPLIT/Proseg notebook family.
- the updated per-celltype migration rules from an internal issue
  (internal issue review):
  - **Rule 1** — if the primary's cell count for a celltype exceeds
    `rule1_threshold` (default 100), skip donor migration for that
    celltype (`primary_only`).
  - **Rule 2** — if the primary is missing a celltype declared in the
    marker JSON, borrow from same-primary-tumor-type donors that
    ACTUALLY HAVE that celltype (capped at `donor_borrow_cap`, default
    100). If no donor has it, emit a warning and skip.
  - Intermediate case (primary has 1 to `rule1_threshold` cells for
    a celltype) supplements the primary via the classical
    `donor_balanced_sample_by_reference` recipe from the ref-build
    summary.
- `primary_only_celltypes` (default `[tumor, liver]`) guards the
  summary's Caveat §1 union-not-intersection leak — those celltypes
  always come exclusively from the primary regardless of the census
  decision.
- `census.csv` — per-celltype audit trail written by the `census`
  stage: `celltype`, `primary_count`, `per_donor_counts` (JSON),
  `decision`, `source_samples`. The operator-visible "why did we
  borrow X from Y" record.
- `ref_build/r/rctd_reference_build.R` — Stage C R script. Reads the
  10X-style bundle written by `export_mtx` via
  `Seurat::ReadMtx(feature.column=1, cell.column=1)`, joins the
  metadata CSV, replaces `/` in cell-type labels with `_` (spacexr
  factor-level constraint), calls `spacexr::Reference(counts,
  cell_types, min_UMI=10, require_int=TRUE)`, and `saveRDS`s the
  result to `<output_root>/<sample_id>/rctd_reference/<sample_id>_scRNA_ref.rds`.
- `noGeneFilter` invariant preserved — no `sc.pp.filter_genes(min_cells=20)`
  anywhere on the reference-build path.
- Raw integer counts on export — the `assemble` stage refuses to cast
  a non-integer `.layers["counts"]` (or `.X` fallback); `export_mtx`
  writes `int32` MatrixMarket so `require_int=TRUE` on the R side
  passes.
- Config machinery split into `ref_build.config` (YAML load + deep
  merge + validation).
- CLI split into `ref_build.cli` (argparse) + `ref_build.pipeline`
  (stage-orchestration loop with sentinel-file resume and `--stages`
  filter).
- `pyproject.toml` + `environment.yml` describing the `refBuild`
  conda env.
- `scripts/submit.slurm.sh` — sbatch wrapper carrying the same
  Slurm-tee/pipefail deadlock fix that landed in the H&E pipeline
  and xenium-preprocess.
- `CITATION.cff`, MIT `LICENSE`, `docs/` markdown pages, smoke tests
  + rule-1 / rule-2 / rule-2-edge unit tests under `tests/`.
