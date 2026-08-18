# ref-build

**ref-build of the Xenium spatial-data preprocessing pipeline: primary + donor-pool preprocessed scRNA AnnData → per-celltype adaptive migration → 10X-style mtx/features/barcodes bundle → spacexr Reference `.rds` for RCTD.**

At a glance:

- **`load_primary_and_donors`** — reads the primary sample's preprocessed scRNA h5ad (originally produced by step 2 `flex-preprocess`; step 2 was **[deprecated 2026-08-10](internal issue review)** and the primary h5ad now comes from manual scRNA preprocessing, meeting the input contract documented under [Configuration](#configuration)) plus a list of same-primary-tumor-type donor files, concatenates with `anndata.concat(..., join="outer", fill_value=0)`, and attaches a `sample_ID` column.
- **`census`** — reads the celltype-marker JSON (same shape as xenium-preprocess's `--global-non-tumor-json`), extracts the expected-celltype set, tallies per-sample per-celltype counts, decides the per-celltype migration policy, and writes `census.csv` — the operator-visible "why did we borrow X from Y" audit trail.
- **`assemble`** — applies the census's per-celltype policy to build the merged reference AnnData. Preserves raw integer counts (no normalization) so the R side's `require_int=TRUE` passes.
- **`export_mtx`** — writes a 10X-style `{counts.mtx.gz, features.tsv.gz, barcodes.tsv.gz, metadata.csv}` bundle. Single-column `features.tsv` so `Seurat::ReadMtx(feature.column=1)` on the R side lines up.
- **`rctd_reference_build`** — shells out to `Rscript` (via the `fhR/4.4.1-foss-2023b` R module (Lmod) by default) to build a `spacexr::Reference` object with `min_UMI=10, require_int=TRUE`, `saveRDS`'d as `<sample_id>_scRNA_ref.rds`. That RDS is what RCTD's `create.RCTD(spatial, reference)` takes as the reference side.
- Every stage is idempotent — sentinel-file resume — and every run snapshots its resolved configuration to disk.

Table of contents: [What it does](#what-it-does) · [Installation](#installation) · [Quickstart](#quickstart) · [Outputs](#outputs) · [Configuration](#configuration) · [HPC (Slurm)](#hpc-usage-slurm) · [Citation](#citation) · [License](#license) · [Acknowledgements](#acknowledgements)

## What it does

Five stages, run in order (`load_primary_and_donors → census → assemble → export_mtx → rctd_reference_build`). Restrict a run to a subset with `--stages`; each stage is guarded by a sentinel file, so re-running with the same `--sample-id` and `--output-root` picks up where the last run left off. Nuke sentinels with `--force-rerun`.

### the per-celltype migration rules (internal issue review)

For each celltype declared in the marker JSON:

- **Rule 1** — if the primary has more than `rule1_threshold` (default 100) cells of that celltype, use ONLY the primary's cells for that celltype. Do not migrate cells from other donors.
- **Rule 2** — if the primary has zero cells of that celltype, borrow cells from same-primary-tumor-type donors that HAVE THAT CELLTYPE, capped at `donor_borrow_cap` (default 100). If no donor has the celltype, emit a warning and continue without.
- **Intermediate case** — 1 to `rule1_threshold` cells in the primary: keep every primary cell, supplement with donor cells via the classical `donor_balanced_sample_by_reference` recipe (per-donor cap = primary's per-celltype count, deterministic under `random_state=1`).

Two additional guards:

- **`primary_only_celltypes`** (default `[tumor, liver]`) — these celltypes ALWAYS come exclusively from the primary sample, regardless of the census decision. Guards against Caveat §1 in the ref-build summary (the union-not-intersection recombine).
- **`noGeneFilter`** — no `sc.pp.filter_genes(min_cells=20)` anywhere on the reference-build path. This preserves the full Xenium panel per the pipeline directive.

### Census variants — `original` vs `balanced` (internal issue review)

Under the classical `original` variant, the intermediate case supplements each donor up to `primary_count` cells, so the **total** donor supplement scales linearly with the number of donors (2 donors ⇒ up to `2 × primary_count`, 4 donors ⇒ up to `4 × primary_count`). the user asked for a variant that keeps the total donor contribution invariant with respect to donor count, so donor-count changes don't inflate donor bias.

Selectable per-run via `--census-variant original|balanced` (or `census.variant` in a config YAML). Same `census.csv` schema so runs can be diffed cell-for-cell.

| Variant | Intermediate case (primary in `[1, rule1_threshold]`) | Primary = 0, donors have it |
|---|---|---|
| `original` (default) | `decision=balanced`; each donor contributes up to `primary_count` cells. Total donor ≤ `primary_count × n_donors`. | `decision=borrowed`; `donor_borrow_cap` total, round-robin. |
| `balanced` | `decision=balanced_borrow`; total donor budget = `primary_count`, split as `per_donor_target = primary_count // n_donors`. Each donor contributes `min(per_donor_target, donor_count)` — **shortfall is not redistributed** (decision (A)). | `decision=borrowed`; hybrid fallback to the same `donor_borrow_cap` round-robin sampler (decision (C)). |

Worked example — primary has 60 fibroblasts, 3 donors have 30, 20, and 5:

- `original`: donor supplement is up to `min(30, 60) + min(20, 60) + min(5, 60) = 55` cells (each donor capped by `primary_count=60`).
- `balanced`: `per_donor_target = 60 // 3 = 20`. Supplement is `min(20, 30) + min(20, 20) + min(20, 5) = 20 + 20 + 5 = 45` cells. Donor 3's shortfall of 15 is NOT redistributed to donor 1 or 2.

### Rule 5 — celltype completeness + fallback donors (internal issue review)

Flex data commonly misses B/Plasma, T/NK, and rbc across an entire sample set — no primary or regular donor has them. Rule 5 opts in to an external "fallback donor" h5ad that's consulted ONLY when a target-list celltype has 0 cells in primary + all regular donors. Fallback donors are IGNORED by Rules 1-4 (so a big immune atlas can't drown out a real primary), and each rescue is capped per-celltype at `fallback_borrow_cap` (default 100).

Two knobs enable it:

- `--fallback-donor-h5ad <path>` (repeatable) — the external h5ad(s). No fallback donor supplied ⇒ Rule 5 is inert (backwards compatible; the default).
- `--celltype-target-list <path>` (optional) — declares the CANONICAL celltype list to enforce completeness over. Defaults to reusing `--celltype-marker-json`'s keys. Add `--celltype-target-key <str>` if the JSON is nested `{tissue: {celltype: markers}}`.

Under Rule 5, `census.csv` gets a new `per_fallback_counts` column (JSON) and a new decision `fallback_borrow` whose `source_samples` names the fallback (`<fallback>fallback:<basename>`). Assemble uses the round-robin borrow sampler across fallback donors that have the celltype, capped at `fallback_borrow_cap` total.

### `load_primary_and_donors` — primary + donor pool → concat AnnData

Adapted from the shared Stage-B backbone in `the internal reference summary` (Stage B, steps 1-3). Reads the primary sample's `_preprocessed_scRNA.h5ad` and each donor's `_preprocessed_scRNA.h5ad` (paths supplied on the CLI as `--primary-h5ad` + one or more `--donor-h5ad`), concatenates them with `anndata.concat(objs, join="outer", fill_value=0, keys=[sample_ids])`, and attaches a `sample_ID` column indicating which sample each cell came from.

**Writes:** `<output_root>/<sample_id>/loaded/concat.h5ad`.

### `census` — expected-celltype set → per-celltype migration decision

Reads the marker JSON (default: `--celltype-marker-json /data/markers/markers.json`), extracts the expected celltype set as `{k.rsplit("_marker", 1)[0] for k in markers}` (same shape as xenium-preprocess's `--global-non-tumor-json`), tallies per-sample per-celltype cell counts, and decides the per-celltype migration policy under the rules:

| Primary count | Any donor has it? | Any fallback has it? | Decision |
|---|---|---|---|
| > `rule1_threshold` (100) | ignored | ignored | `primary_only` |
| in `[1, rule1_threshold]` | any | ignored | `balanced` (use primary + donor supplementation) |
| 0 | yes | ignored | `borrowed` (borrow from donors that have it, capped by `donor_borrow_cap`) |
| 0 | no | yes | `fallback_borrow` — Rule 5; borrow from fallback donor(s), capped by `fallback_borrow_cap` |
| 0 | no | no | `missing_no_donor` (warn + continue without) |

`primary_only_celltypes` (default `[tumor, liver]`) overrides the decision to `primary_only` regardless of the count.

**Writes:** `<output_root>/<sample_id>/census/census.csv` — one row per celltype, columns `celltype, primary_count, per_donor_counts (JSON), per_fallback_counts (JSON), decision, source_samples, note, matched_labels`.

### `assemble` — apply census → merged reference AnnData

Applies the census's per-celltype policy to build the merged reference AnnData. For each celltype:

- `primary_only` — every primary cell of that celltype passes through.
- `balanced` — every primary cell passes through; donor cells are added via `donor_balanced_sample_by_reference` (per-donor cap = primary's per-celltype count, `random_state=1`).
- `balanced_borrow` — every primary cell passes through; donor cells are added via `donor_balanced_borrow_sample` (per_donor_target = primary_count // n_donors, no shortfall redistribution).
- `borrowed` — donor cells are added up to `donor_borrow_cap` total; primary contributes zero (by definition — primary has 0 of that celltype).
- `fallback_borrow` — fallback-donor cells are added up to `fallback_borrow_cap` total, round-robin across fallbacks that have the celltype (Rule 5).
- `missing_no_donor` — celltype is skipped, warning logged.

Preserves raw integer counts. If `.layers["counts"]` is not integer-valued, the stage refuses to proceed rather than silently truncating (the R side's `require_int=TRUE` needs integer counts).

**Writes:** `<output_root>/<sample_id>/assembled/reference.h5ad`.

### `export_mtx` — merged reference → 10X-style bundle

Writes the four 10X-style files that both the Stage C R script and any downstream RCTD workflow can read:

- `<sample>_counts.mtx.gz` — MatrixMarket, gzipped, cast to `int32`, transposed to `genes × cells`.
- `<sample>_features.tsv.gz` — one gene name per line (single-column form; matches `Seurat::ReadMtx(feature.column=1)`).
- `<sample>_barcodes.tsv.gz` — one cell id per line.
- `<sample>_metadata.csv` — `adata.obs`, cell-indexed, plain (uncompressed) CSV.

**Writes:** `<output_root>/<sample_id>/mtx_bundle/`.

### `rctd_reference_build` — 10X bundle → spacexr Reference `.rds` (Stage C)

Python (`ref_build.stages.rctd_reference_build`) shells out to `Rscript src/ref_build/r/rctd_reference_build.R`. The R script:

1. `Seurat::ReadMtx(mtx=…, features=…, cells=…, feature.column=1, cell.column=1)`.
2. Reads the metadata CSV, joins by cell id, extracts the celltype column (default `Final_level1_celltype_annotation`).
3. Replaces `/` characters in cell-type labels with `_` — spacexr factor levels don't accept slashes (e.g. `B/Plasma_T/NK_rbc` → `B_Plasma_T_NK_rbc`).
4. `ref <- spacexr::Reference(counts=GetAssayData(seu, "RNA", "counts"), cell_types=labels, min_UMI=10, require_int=TRUE)`.
5. `saveRDS(ref, "<sample_id>_scRNA_ref.rds")`.

The resulting RDS is the RCTD "reference object" — downstream `spacexr::create.RCTD(spatial, reference)` takes it as the reference side of the deconvolution.

**Writes:** `<output_root>/<sample_id>/rctd_reference/<sample_id>_scRNA_ref.rds`.

## Installation

`ref-build` targets Python 3.10+ (validated on 3.11). Stage C shells out to `Rscript` and needs `Seurat`, `Matrix`, `readr`, and `spacexr`.

```bash
# 1. Install micromamba (skip if you already have it)
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# 2. Create the env
micromamba env create -f environment.yml
micromamba activate refBuild

# 3. Editable install
pip install -e /path/to/ref-build

# 4. Verify
ref-build --help
ref-build run --help
```

For the R side, see [`docs/install.md`](docs/install.md) — the recommended path is `ml fhR/4.4.1-foss-2023b` + a one-time `spacexr` install into `~/.claude/r_libs/4.4.1`.

## Quickstart

```bash
micromamba activate refBuild
ml fhR/4.4.1-foss-2023b

ref-build run \
    --sample-id            SAMPLE1 \
    --primary-h5ad         /data/SAMPLE1/scRNA/SAMPLE1_preprocessed_scRNA.h5ad \
    --donor-h5ad           /data/SAMPLE2/scRNA/SAMPLE2_preprocessed_scRNA.h5ad \
    --donor-h5ad           /data/SAMPLE3/scRNA/SAMPLE3_preprocessed_scRNA.h5ad \
    --celltype-marker-json /data/markers/markers_global_no_tumor_cells_level1.json \
    --tumor-type           Bladder \
    --output-root          /data/ref_build_runs
```

Skip a stage set with `--stages`:

```bash
# Regenerate only the .rds (assumes the 10X bundle exists).
ref-build run --sample-id SAMPLE1 --primary-h5ad … --donor-h5ad … \
    --celltype-marker-json … --output-root … \
    --stages rctd_reference_build
```

## Outputs

Under `<output_root>/<sample_id>/`:

```
loaded/concat.h5ad                            — stage 1
census/census.csv                             — stage 2 (the audit trail)
assembled/reference.h5ad                      — stage 3
mtx_bundle/<sample>_counts.mtx.gz             — stage 4
mtx_bundle/<sample>_features.tsv.gz           — stage 4
mtx_bundle/<sample>_barcodes.tsv.gz           — stage 4
mtx_bundle/<sample>_metadata.csv              — stage 4
rctd_reference/<sample>_scRNA_ref.rds         — stage 5 (the RCTD-ready reference)
config.yaml                          — snapshot of the resolved config
```

## Configuration

Every knob lives in `config/default.yaml`. Override precedence: CLI flag > user YAML (`--config …`) > default YAML.

Required-at-runtime keys (all `null` in the default YAML — supplied via CLI or user YAML):

| Key | Meaning |
|---|---|
| `sample_id` | Primary sample id (e.g. `SAMPLE1`). Used to name the output subdirectory + the .rds. |
| `primary_h5ad` | Path to the primary's preprocessed scRNA h5ad. Historically produced by step 2 (`flex-preprocess`, [deprecated 2026-08-10](internal issue review); now typically manual scRNA prep. Must meet the input contract: `.obs[celltype_col]` populated (default `Final_level1_celltype_annotation`); counts source `layers["counts"]` → `layers["raw_count"]` → `.X` (first present wins) and MUST be integer (`spacexr::Reference(..., require_int=TRUE)`). |
| `output_root` | Root output directory. |
| `celltype_marker_json` | Marker-gene JSON declaring the expected celltype set. |

Optional (with sensible defaults):

| Key | Default | Meaning |
|---|---|---|
| `donor_h5ads` | `[]` | List of donor preprocessed-scRNA h5ad paths (same primary tumor type). Same input contract as `primary_h5ad`. |
| `tumor_type` | `null` | Cosmetic label written into census provenance + resolved-config snapshot. |
| `celltype_col` | `Final_level1_celltype_annotation` | Per-cell celltype column in `.obs`. |
| `rule1_threshold` | `100` | Primary count above which we don't migrate donors (the rule 1). |
| `donor_borrow_cap` | `100` | Max donor cells to borrow per missing celltype (the rule 2). |
| `primary_only_celltypes` | `[tumor, liver]` | Always primary-only, regardless of census decision. |
| `random_state` | `1` | Seed for `donor_balanced_sample_by_reference`. |
| `census.variant` | `original` | `original` (per-donor cap, classic) or `balanced` (revised rule: total donor = primary_count, split per-donor, no shortfall redistribution). |
| `min_UMI` | `10` | Stage C — cells below this UMI count are dropped by `spacexr::Reference`. |
| `require_int` | `true` | Stage C — `spacexr::Reference` refuses non-integer counts. |

Full CLI reference: `ref-build run --help`. Semantics of each knob: `docs/usage.md`.

## HPC usage (Slurm)

```bash
sbatch scripts/submit.slurm.sh SAMPLE1 \
    /data/SAMPLE1/scRNA/SAMPLE1_preprocessed_scRNA.h5ad \
    --donor-h5ad /data/SAMPLE2/scRNA/SAMPLE2_preprocessed_scRNA.h5ad \
    --donor-h5ad /data/SAMPLE3/scRNA/SAMPLE3_preprocessed_scRNA.h5ad \
    --celltype-marker-json /data/markers/markers_global_no_tumor_cells_level1.json \
    --tumor-type Bladder \
    --output-root /data/ref_build_runs
```

The submit wrapper carries the Slurm-tee/pipefail deadlock fix (see `scripts/submit.slurm.sh` for the comment block).

## Citation

If you use `ref-build` in published work, cite it via the metadata in [`CITATION.cff`](CITATION.cff). Also cite the upstream tools this pipeline glues together:

- **scanpy** — Wolf, F.A. *et al.* Genome Biology 19, 15 (2018).
- **anndata** — Virshup *et al.* bioRxiv 2021.12.16.473007 (2021).
- **Seurat** — Hao *et al.* Cell 184, 3573–3587.e29 (2021).
- **RCTD / spacexr** — Cable *et al.* Nat. Biotechnol. 40, 517–526 (2022).

## License

MIT — see [`LICENSE`](LICENSE). See the "Third-party attribution" block in `LICENSE` for the internal SPLIT/Proseg workflow this package is a port of.

## Acknowledgements

Adapted from an internal SPLIT/Proseg workflow with updated per-celltype migration rules from issue (internal issue review). See `docs/methods.md` for a per-stage attribution to the source notebooks and scripts.
