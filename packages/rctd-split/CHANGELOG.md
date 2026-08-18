# Changelog

All notable changes to rctd-split will be documented here. Follows
[Keep a Changelog](https://keepachangelog.com/) shape.

## [Unreleased]

### Added
- `qc_report`: optional `extra_reports` list appends external HTML
  reports to the bottom of `summary_report.html` (Tracy's ask on
  `settylab/TracyY123-nexus#26` comment 5322126093, item 4).
  Each entry is a `{path, name}` dict; the report renders each as
  an `<h3>` heading + click-through `<a>` link + sandboxed inline
  iframe. When the linked file lives under `summary/`, the path
  is rewritten as relative so the folder stays portable; otherwise
  the absolute path is used. Missing files render a visible
  "missing at render time" note rather than a broken link.
  Interfaces:
  - Config YAML: `qc_report.extra_reports: [{path: ..., name: ...}, ...]`.
    Wire in via `--config` on `rctd-split run` or `--step4-config`
    on the workflow driver.
  - CLI: `rctd-split run --extra-report PATH,NAME` (repeatable;
    first-comma split lets the display name contain commas).
- `qc_report`: fixed color maps for the `purification_status`
  and `first_type` categorical UMAPs so colors stay consistent
  across samples (Tracy's ask on `settylab/TracyY123-nexus#26`
  comment 5322126093, item 3).
  - `purification_status`: high-contrast Wong-derived palette;
    `singlet` / `purified` → dark green, `reject` / `discarded`
    → dark magenta, doublet variants get distinct hues, unknowns
    fall through to a fallback cycle.
  - `first_type` (SPLIT-inferred celltype): deterministic
    per-name SHA-1 hash into the `tab20` palette. The SAME
    celltype name gets the SAME hex color across every sample
    the pipeline runs — no shared global celltype list required.
  - Full color mapping is persisted to
    `summary/<sample>_color_map.json` and linked from the two
    UMAP captions in the HTML report for reproducibility.
  - Palette lives in `rctd_split._internal.palette`;
    `_scatter_umap` grew an optional `color_map` kwarg.

### Changed
- `qc_report` HTML output filename renamed
  `<sample>_qc_report.html` → `<sample>_summary_report.html`
  (title + `<h1>` updated to "Summary report" to match). The
  `summary/` folder is where the HTML lives; the filename now
  matches. Referenced by
  `_QC_BASENAMES["html_report"]` in `_internal/layout.py`;
  no downstream reader relies on the old name.
  Requested via
  `settylab/TracyY123-nexus#26` comment 5322126093 (item 1).

### Added
- `qc_report` now extracts an RCTD spot_class summary + a first_type
  breakdown of RCTD-rejected cells from `raw.obs` (`spot_class` +
  `first_type`, folded from `step4_unpurified.obs` by
  `writeback_to_step1_raw` Part 2; originally emitted by
  `SPLIT::run_post_process_RCTD` from RCTD's `results_df`). Renders a
  two-part table in the HTML report:
    1. RCTD spot_class counts + % of RCTD-categorized (singlet,
       doublet_certain, doublet_uncertain, reject; any unexpected
       values are appended in sorted order rather than dropped).
    2. `first_type` breakdown within `spot_class == 'reject'` cells.
  Also writes `qc/<sample_id>_rctd_summary.csv` as a machine-parseable
  sidecar. Fail-loud when `raw.obs['spot_class']` or
  `raw.obs['first_type']` is absent — points at
  `writeback_to_step1_raw` as the responsible upstream stage.
  Requested via internal review comment (internal).

### Changed
- `qc_report` histograms now plot the distribution of ALL cells in
  each source h5ad rather than the previous positive-only slice. A
  dashed vertical line is drawn at the upstream filter threshold read
  from the merged `config.yaml`:
    - `proseg_raw` histogram: `step1.qc_filter.min_counts_cell`
    - `xenium_ranger` histogram: no line (step-1 has no min-counts
      gate on xenium_ranger)
    - `proseg_purified` histogram: `step4.postprocess.qc.min_counts`
  Thresholds are strictly data-driven — if `config.yaml` is
  absent or a key is missing, the caption falls back to "no dashed
  threshold line drawn" (never a hard-coded value). Cells with counts
  ≤ 0 are drawn as a separate hatched leftmost bar so no cell is
  dropped from the visual. Fixes the user's "why is the html showing
  only non-positive cells" complaint on
  (internal issue review) comment (internal).
- `rctd_run.CELL_MIN_INSTANCE` default lowered from `25` (the Rmd
  value / spacexr default) to `20`. Motivation: detailed-annotation
  references (e.g. `Pt32_scRNA_ref_new.rds`) have smaller per-type
  cell counts, and the detailed-annotation rerun already
  used `20` for the same reason (an internal source).
  Requested via internal review (2026-07-28). Existing
  invocations that pass `--cell-min-instance N` or set
  `rctd_run.CELL_MIN_INSTANCE` in a YAML override continue to work
  unchanged; only the fall-through default moves.

### Fixed
- `qc_report.raw_layer` default changed from `maxpost` → `maxpost_counts`
  to match the canonical layer name proseg's step-1 export actually
  emits (`expected_counts`, `maxpost_counts`). users hit
  `SystemExit: proseg_raw.h5ad has no layers['maxpost']` at runtime
  on (internal issue review) comment (internal). Touches:
  `config/default.yaml`, `pipeline.py`, `cli.py`, `stages/qc_report.py`
  (docstring + error message now points at `maxpost_counts` as the
  canonical name), and `tests/test_qc_report.py`. Existing runs that
  set `--qc-raw-layer maxpost` or `qc_report.raw_layer: maxpost` in a
  YAML override continue to work only if the h5ad actually has that
  layer; production h5ads use `maxpost_counts`.
- `rctd_run.R` now propagates its full `.libPaths()` back into
  `R_LIBS_USER` before invoking `run.RCTD`, so the PSOCK workers
  `spacexr::process_beads_batch` / `decompose_batch` spawn inherit the
  same library search path. Previously the parent-only in-memory
  `.libPaths()` prepend was invisible to workers, and workers died with
  `there is no package called 'spacexr'` during closure deserialisation
  whenever spacexr lived on a path (user-local install, an renv library)
  that wasn't already on the inherited `R_LIBS_USER`. Pre-registering a
  cluster with `doParallel::registerDoParallel(cl)` after
  `parallel::clusterCall(cl, .libPaths, .libPaths())` was tried and
  ruled out: spacexr unconditionally calls `parallel::makeCluster` and
  re-registers its own backend. See `docs/methods.md` §"`max_cores` and
  PSOCK worker library paths". Config default `rctd_run.max_cores`
  reverted to `4` (Rmd value).

## [0.1.0] — 2026-07-09

### Added
- Initial package scaffold (step 4 of the internal SPLIT/Proseg workflow spatial-data
  preprocessing pipeline).
- Four stage modules under `rctd_split.stages`:
  - `rctd_run` — Python wrapper → `r/rctd_run.R`: `create.RCTD` +
    `run.RCTD(doublet_mode="doublet")`, saves `rctd_results.rds`.
  - `split_purify` — Python wrapper → `r/split_purify.R`:
    `SPLIT::run_post_process_RCTD` (unpurified variant) then
    `SPLIT::purify(..., DO_purify_singlets=TRUE)` (purified variant);
    saves `unpurified.rds` and `purified.rds`.
  - `export_mtx` — Python wrapper → `r/export_mtx.R`: reads both RDS
    variants and writes matched 10X-style bundles
    (`{counts.mtx.gz, features.tsv.gz, barcodes.tsv.gz, metadata.csv,
    spatial_coords.csv.gz}`).
  - `mtx_to_h5ad` — pure Python: reads each bundle, builds AnnData with
    `.var` = features, `.obs_names` = barcodes, merges `metadata.csv`
    into `.obs`, attaches spatial coords to `.obsm["spatial"]`; writes
    `<sample>_unpurified.h5ad` and `<sample>_purified.h5ad`.
- Config machinery `rctd_split.config` (YAML load + deep merge +
  validation).
- CLI split into `rctd_split.cli` (argparse) + `rctd_split.pipeline`
  (stage-orchestration loop with sentinel-file resume and `--stages`
  filter).
- `pyproject.toml` + `environment.yml` describing the `rctdSplit`
  conda env.
- `scripts/submit.slurm.sh` — sbatch wrapper carrying the same
  Slurm-tee/pipefail deadlock fix that landed in step 1
  (the xenium-preprocess submit script).
- `CITATION.cff`, MIT `LICENSE`, `docs/` markdown pages, smoke tests
  under `tests/` plus a synthetic round-trip test for `mtx_to_h5ad`.

### Ported from
- the internal SPLIT/Proseg workflow —
  the canonical Rmd driving RCTD + SPLIT for the MH sample family.
  Every parameter (`UMI_min=10`, `counts_MIN=10`, `UMI_min_sigma=100`,
  `max_cores=4`, `CELL_MIN_INSTANCE=25`, `doublet_mode="doublet"`,
  `DO_purify_singlets=TRUE`) is inline-commented back to a line in the
  Rmd.
- `the internal reference summary` Stage D (lines 282-298) — the
  workflow narrative + the metadata columns downstream (`first_type`,
  `purification_status`, `w1_larger_w2`, `same_class`, `nCount_Proseg`,
  spatial `x`/`y`) that this pipeline propagates into the h5ad's
  `.obs`.
