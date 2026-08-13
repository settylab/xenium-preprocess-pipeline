# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
- `scripts/submit_workflow.sh` — Slurm driver chaining steps 1 → 3 → 4 with
  `--dependency=afterok:`.
- `scripts/submit_step{1,3,4}.sbatch` — per-step Slurm templates.
- `docs/` — installation, quickstart, and pipeline-overview docs.
- `environments/xenium.yml` — micromamba environment spec for the Python side.
- Per-package `pyproject.toml`, `README.md`, `CHANGELOG.md`, `CITATION.cff`,
  `LICENSE`, and `tests/`.
