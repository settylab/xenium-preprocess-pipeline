# Methods

`ref-build` is a stage-orchestrated wrapper around the internal SPLIT/Proseg workflow's per-sample scRNA reference-dataset build, with the per-celltype migration rules updated per issue (internal issue review). Behaviour is a port of the reference notebook family and the summary document; this file points each stage back at its source.

## Stage 1: `load_primary_and_donors`

**Source:** `the internal reference summary` Stage B, steps 1-3 (lines 76-79).

Reads the primary sample's `_preprocessed_scRNA.h5ad` and each donor's `_preprocessed_scRNA.h5ad`, concatenates them with `anndata.concat(objs, join="outer", fill_value=0, keys=[sample_ids], label="sample_ID", index_unique=None)`. Outer-join preserves the union of genes so mismatched panels get 0-filled (the internal SPLIT/Proseg workflow summary documents outer-join keeps the union). The `keys` list drives the `sample_ID` obs column that downstream stages read.

The `noGeneFilter` invariant (the pipeline directive; summary lines 30, 80, 237) is preserved from load onward — no `sc.pp.filter_genes(min_cells=20)` anywhere.

## Stage 2: `census`

**Source:** An internal reference recipe (Stage B), but with **updated per-celltype migration rules** from an internal issue (internal issue review):

- **Rule 1** — if the primary's cell count for a celltype exceeds `rule1_threshold` (default 100), decision is `primary_only` — no donor migration for that celltype.
- **Rule 2** — if the primary has zero cells for a celltype:
    - **AND** at least one donor has the celltype → decision is `borrowed` (borrow from donors that HAVE THAT celltype, capped by `donor_borrow_cap`, default 100).
    - **OR** no donor has the celltype → decision is `missing_no_donor` — WARN log line + skip.
- **Intermediate case** (1 to `rule1_threshold` primary cells) → decision is `balanced` — keep every primary cell, supplement with donor cells via `donor_balanced_sample_by_reference` (per-donor cap = primary's per-celltype count; `random_state=1` — matches summary line 89-91).

Two additional guards:

- **`primary_only_celltypes`** (default `[tumor, liver]`) — always `primary_only` regardless of the count. Guards summary Caveat §1 (lines 358-369): the summary's Stage B step 10 recombine is union-not-intersection on `celltype.isin({"tumor","liver"})`; if a donor tumor cell survived the earlier subset, it would silently sweep in. `primary_only_celltypes` makes the guard explicit rather than relying on an earlier subset having removed the donor-tumor cells.
- **Expected-celltype set** — the marker JSON's keys with `_marker` stripped (same shape as step 1's `--global-non-tumor-json`; the expected set is fixed by the JSON, per the user's clarification).

### Fuzzy celltype-label matching

From user request 2026-07-10 (internal issue review) — the census stage no longer requires an exact `obs[celltype_col] == X` match to count a cell toward expected celltype `X`. Instead, three match rules apply (any-of):

1. **Exact** — `L == X`. Always active.
2. **Unknown-maybe** — `L == f"unknown_maybe_{X}"`. Active when `census.include_unknown_maybe: true` (default). Rationale: cells labelled `unknown_maybe_Fibroblast` are Fibroblast candidates whose annotation carries some uncertainty; including them lets downstream spatial deconvolution reduce that uncertainty against a larger spatial panel.
3. **Composite delimited** — `X` appears in `L` as a token surrounded by `_` / `/` / string boundary. Regex: ``(?:^|[_/])re.escape(X)(?:$|[_/])``. Active when `census.fuzzy_matching: true` (default). Rationale: preprocessing (step 3 of the summary) produces composite labels like `B/Plasma_T/NK_rbc` for cells whose marker signal maps to multiple lineages; these should count toward every constituent celltype, not be silently dropped.

The `re.escape(X)` is essential because expected celltypes may themselves contain `/` (e.g. `T/NK`).

A cell may match multiple `X`'s (composite labels do so by design). The **assemble** stage dedupes selected cells by cell id at the end, so a composite-labelled cell appears **once** in the final reference, with its **original label preserved** (no relabeling — RCTD/`spacexr::Reference` consume any label set).

For **backward compatibility**, set both `census.fuzzy_matching: false` and `census.include_unknown_maybe: false` to revert to the pre-2026-07-10 strict-equality behavior.

Two extra columns land in `census.csv`:

- **`matched_labels`** — comma-separated list of actual on-cell labels that were rolled up under this expected celltype. Empty for `missing_no_donor` rows. Example: for `celltype=T/NK`, `matched_labels` might be `T/NK,unknown_maybe_T/NK,B/Plasma_T/NK_rbc`. The audit trail lets an operator answer "which composite/unknown_maybe labels drove this count?" without re-running.

`census.csv` is the operator-visible audit trail: `celltype, primary_count, per_donor_counts (JSON), decision, source_samples, note, matched_labels`. One row per expected celltype. Human-readable — inspect it to answer "why did we borrow X from Y" without re-running.

## Stage 3: `assemble`

**Source:** `the internal reference summary` Stage B, steps 9-10 (lines 85-87), with the recombine step rewritten to use the census decisions instead of the summary's `mask = obs_names.isin(keep_ids) | celltype.isin({"tumor","liver"})` union.

For each celltype's census row, applies:

- `primary_only` — keep every primary cell of that celltype.
- `balanced` — keep every primary cell; supplement with donor cells via a helper that faithfully ports `donor_balanced_sample_by_reference` (per-donor cap = primary's per-celltype count, `random_state=1`).
- `borrowed` — distribute the `donor_borrow_cap` across donors that HAVE the celltype using seed-deterministic round-robin (ceiling division per donor, trim from the tail if we overshoot).
- `missing_no_donor` — skip.

Cell selection uses the **same fuzzy match rules** as the census stage (mirrored from `census.fuzzy_matching` + `census.include_unknown_maybe`) — a cell counted in the census MUST be pickable here, or the two would drift. Composite-labelled cells that match multiple expected celltypes are added to the kept-ids list once per matching row; the final **dedupe by cell id** (already the last step of the assembly loop) keeps them exactly once in the merged reference, with their original label untouched.

Preserves raw integer counts. `.layers["counts"]` (or the fallback layer) must be integer — the stage refuses to proceed on a non-integer counts layer rather than silently truncating in the `int32` cast at export time. The Stage-C R side runs `require_int=TRUE`; a float counts layer would fail there anyway, but failing loudly here surfaces the upstream preprocessing bug sooner.

The `noGeneFilter` invariant continues — no gene filter here.

## Stage 4: `export_mtx`

**Source:** `the internal reference summary` Stage B, step 11 (line 87), and the file naming convention from step 1's `split_prep` module.

Writes:

- `<sample>_counts.mtx.gz` — MatrixMarket, gzipped, cast to `int32`, transposed to `genes × cells` (matches `Seurat::ReadMtx` conventions).
- `<sample>_features.tsv.gz` — one gene name per line (single-column form; matches `Seurat::ReadMtx(feature.column=1)`).
- `<sample>_barcodes.tsv.gz` — one cell id per line.
- `<sample>_metadata.csv` — `adata.obs`, cell-indexed, plain (uncompressed) CSV.

The Stage-C R side reads all four files.

## Stage 5: `rctd_reference_build`

**Source:** `the internal reference summary` Stage C (lines 252-281), and the user's rebuild scripts under an internal source + an internal source.

Python (`ref_build.stages.rctd_reference_build`) shells out to `Rscript src/ref_build/r/rctd_reference_build.R`. The R script:

1. Prepends `~/.claude/r_libs/4.4.1` to `.libPaths()` so the user-local `spacexr` install is found (idempotent — the prepend is skipped if the dir doesn't exist).
2. Loads `Seurat`, `Matrix`, `spacexr`.
3. `Seurat::ReadMtx(mtx=…, features=…, cells=…, feature.column=1, cell.column=1)` on the mtx bundle.
4. Reads the metadata CSV, joins by cell id, extracts the celltype column (default `Final_level1_celltype_annotation`).
5. Replaces `/` in celltype labels with `_` (default; configurable) — `spacexr::Reference` factor levels reject slashes (e.g. `B/Plasma_T/NK_rbc` → `B_Plasma_T_NK_rbc`, per summary line 262).
6. Calls
   ```r
   ref <- Reference(
       counts      = mtx,
       cell_types  = ref_labels,
       min_UMI     = 10,
       require_int = TRUE
   )
   ```
7. `saveRDS(ref, "<sample_id>_scRNA_ref.rds")`.

The resulting `.rds` is the RCTD "reference object" — downstream `spacexr::create.RCTD(spatial, reference)` calls take it as the reference side of the deconvolution.

## Reproducibility

- Every run writes a `config.yaml` at the run root, capturing the exact merged config (default YAML + user YAML + CLI overrides).
- The one RNG path — `donor_balanced_sample_by_reference` in the balanced case, and the round-robin cap in the borrowed case — is seeded from `census.random_state` (default 1, matching the notebook family). With the same seed + same inputs, the build is deterministic.
- The `assemble` stage's O(cells) allocation and O(celltypes) census loop mean the runtime scales linearly in both. No per-cell branching costs; the census does the heavy per-celltype decision work once.
- The `.rds` filename `<sample_id>_scRNA_ref.rds` matches the shipped-reference naming in summary Table 1 (lines 340-352) so downstream Stage-D driver code can pick up the reference without a rename.
