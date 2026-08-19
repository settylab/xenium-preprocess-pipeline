# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] — semantic-name migration

### Changed

- Renamed the per-step Slurm stubs to semantic names, matching the
  package CLIs (`xenium-preprocess` / `ref-build` / `rctd-split`):
  - `scripts/submit_step1.sbatch` → `scripts/submit_xenium-preprocess.sbatch`
  - `scripts/submit_step3.sbatch` → `scripts/submit_ref-build.sbatch`
  - `scripts/submit_step4.sbatch` → `scripts/submit_rctd-split.sbatch`
  `submit_workflow.sh` invokes the new names; all `--job-name`, log-tag,
  and internal-tmp-file identifiers inside the stubs use the semantic
  names too. Direct-`sbatch` callers (rare — the documented flow is via
  `submit_workflow.sh`) need to update the filename.

## [0.1.0] — Initial migration

Initial standalone-repo release, migrated from an internal SPLIT/Proseg workflow tree.
Bundles three previously ad-hoc Python packages plus a shell driver so a new
Xenium sample runs end-to-end from a single command.

### Added

- `packages/xenium-preprocess/` — proseg → adata preprocessing (5 stages).
  Ships `xenium-preprocess` CLI.
- `packages/ref-build/` — celltype-marker-driven reference build. Ships
  `ref-build` CLI.
- `packages/rctd-split/` — RCTD + SPLIT typing (4 stages). Ships `rctd-split`
  CLI.
- `scripts/submit_workflow.sh` — Slurm driver chaining
  xenium-preprocess → ref-build → rctd-split with `--dependency=afterok:`.
- `scripts/submit_{xenium-preprocess,ref-build,rctd-split}.sbatch` — per-step
  Slurm templates. (Renamed from the original `submit_step{1,3,4}.sbatch` in
  the `v0.2.0` semantic-migration; see the top-level CHANGELOG entry below.)
- `docs/` — installation, quickstart, and pipeline-overview docs.
- `environments/xenium.yml` — micromamba environment spec for the Python side.
- Per-package `pyproject.toml`, `README.md`, `CHANGELOG.md`, `CITATION.cff`,
  `LICENSE`, and `tests/`.
