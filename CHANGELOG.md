# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.7] — 2026-10-01

### Added

- **New sibling package `packages/exploration/`** (Tracy `#31`) — a home
  for downstream / optional analyses that run AFTER the canonical
  `xenium-preprocess` output. Peer of `xenium-preprocess`, `rctd-split`,
  and `ref-build` under `packages/`, same conventions (own
  `pyproject.toml`, `src/` layout, `tests/`, `environment.yml`,
  `README.md`, `CHANGELOG.md`, `LICENSE`, `CITATION.cff`). See
  `packages/exploration/README.md` for the full conceptual boundary
  (preprocess = stable canonical producer; exploration = downstream
  question-specific analyses). Initial version `0.1.0`.
- `exploration.annotate.annotate_regions_and_tumor_distance` — region +
  signed tumor-distance annotator. Takes a preprocessed AnnData (`SPLIT`
  output with `obs['first_type']` + `obsm['spatial']` in µm) and a
  user-authored region CSV (`Selection`, `X`, `Y`, optional `Class`), and
  returns a NEW AnnData with three strictly additive `.obs` columns:
  `selected_region` (all cells), `tumor_distance_bin` (hepatocytes only;
  `intratumoral` / `0-100um` / `100-300um` / `>300um`), and
  `tumor_boundary_distance_um` (hepatocytes only; signed Euclidean
  distance in µm, negative inside the tumor union, positive outside).
  Tumor-selection resolution uses the F1-patched gate (name-based tumor
  + Class-majority for ambiguous selections; `non-fibrotic tumor`
  correctly stays tumor, `non-tumor` is excluded). Refuses to overwrite
  any of the three target columns if already present. The input adata
  is never mutated. Also exported at package root as
  `exploration.annotate_regions_and_tumor_distance`.
- `packages/exploration/tests/test_annotate.py` — 10 unit tests pinning
  the overlap rule, bin edges, hep-only distances, F1 classifier,
  preservation invariants, and parity against the AUDITED reference
  algorithm in `scratch/tracy-i31-mh-all-metx-workflow/
  run_sample.py:load_distance_bins`.
- `shapely>=2` added to the NEW `exploration` package's `environment.yml`
  + `pyproject.toml` (NOT to `xenium-preprocess` — see "Removed" below).

### Removed

- `xenium_preprocess.annotate.annotate_regions_and_tumor_distance` — the
  function previously added to `xenium-preprocess` in the earlier
  v0.2.7-draft commit `9a8530c` has been **moved out** of
  `xenium-preprocess` into the new `exploration` package. Rationale:
  region/distance annotation is downstream, question-specific, and
  should not block the canonical preprocessing pipeline when no region
  CSV is provided. `xenium-preprocess` now has **zero** dependency on
  the region CSV; the canonical processed `.h5ad` is produced
  identically to v0.2.6. `shapely>=2` is removed from
  `xenium-preprocess`'s deps.
- **Import break (no shims)**: callers that previously used
  `from xenium_preprocess import annotate_regions_and_tumor_distance`
  must switch to
  `from exploration import annotate_regions_and_tumor_distance`.

### Validation (MH6_BR_2328 full sample)

See `validation_mh6.log` under the branch worktree for the full output.
The MH6 validation was run against the pre-move implementation (same
code, same algorithm, same byte-for-byte parity vs the AUDITED
reference); after the move the unit tests (10, now under
`packages/exploration/tests/`) re-pass unchanged.

## [0.2.6] — 2026-08-26

Patch release: defensive env-config resolution + longer default runtime
for large RCTD split jobs. `scripts/*` only — no `packages/*/` code
changes, no CLI-surface changes.

### Fixed

- `scripts/create-env.sh` (`e3108de`): prefer `$MAMBA_EXE` over
  `command -v micromamba` when `--micromamba-bin` isn't given. A
  shadowed / wrong-arch `micromamba` earlier on `PATH` (e.g. a stale
  `~/bin/micromamba` copied from another machine) passes `command -v`
  cleanly, then fails with `Exec format error` at exec — well after
  the guard is happy. `$MAMBA_EXE` is exported by micromamba's
  shell-hook to the exact binary the hook sourced, so it's
  known-working by construction. Cache-isolation semantics preserved
  (`MAMBA_ROOT_PREFIX` pinning, `CONDA_PKGS_DIRS` export,
  `~/.mamba/pkgs` untouched assertion, `record_env_prefix` all
  unchanged). Ref `settylab/TracyY123-nexus#26` comment 5417538704.
- `scripts/write-env-config.sh` (`2cc9179`): companion fix — same
  three-level resolution (`--micromamba-bin` → `$MAMBA_EXE` →
  `command -v micromamba`) applied to the sibling script, catching
  the same shadowed / wrong-arch shape.

### Changed

- `scripts/submit_rctd-split.sbatch` (`a970b97`): bump default
  `#SBATCH --time` from `24:00:00` to `6-00:00:00` (D-HH:MM:SS
  form). 24h is not enough for most TMA-scale RCTD split runs.
  Mirrors the equivalent bump on main (`9b188c0`, pre-semantic-rename
  filename). Ref `settylab/TracyY123-nexus#26` comment 5418040437.

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
