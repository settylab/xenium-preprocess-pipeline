# Usage

## Command-line reference

```
ref-build run [OPTIONS]
```

The four required-at-runtime flags:

| Flag | Meaning |
|---|---|
| `--sample-id`            | Primary sample identifier (e.g. `SAMPLE1`). Names the output subdirectory + the `.rds`. |
| `--primary-h5ad`         | Path to the primary's `_preprocessed_scRNA.h5ad` (from step 2 `flex-preprocess`). |
| `--celltype-marker-json` | Marker-gene JSON declaring the expected celltype set. Same shape as xenium-preprocess's `--global-non-tumor-json`. |
| `--output-root`          | Root output directory. Per-sample results land at `<output_root>/<sample_id>/`. |

Every other flag has a default in `config/default.yaml`. Override precedence: CLI flag > user YAML (`--config user.yaml`) > default YAML.

### Donor pool

Pass one or more donor `_preprocessed_scRNA.h5ad` paths via repeated `--donor-h5ad`:

```bash
ref-build run \
    --sample-id       SAMPLE1 \
    --primary-h5ad    /data/SAMPLE1/scRNA/SAMPLE1_preprocessed_scRNA.h5ad \
    --donor-h5ad      /data/SAMPLE2/scRNA/SAMPLE2_preprocessed_scRNA.h5ad \
    --donor-h5ad      /data/SAMPLE3/scRNA/SAMPLE3_preprocessed_scRNA.h5ad \
    --celltype-marker-json /data/markers/markers_global_no_tumor_cells_level1.json \
    --tumor-type      Bladder \
    --output-root     /data/ref_build_runs
```

If no `--donor-h5ad` is supplied (or `donor_h5ads: []` in a user YAML), the build is a **sample-only** reference — per pipeline directives (single-sample references). The census stage will emit `primary_only` for every celltype that the primary has and `missing_no_donor` for the rest.

### Stage selection + resume

```bash
ref-build run --stages load_primary_and_donors census   # partial pipeline
ref-build run --stages rctd_reference_build              # only stage 5
ref-build run --force-rerun                              # nuke sentinels + re-run
```

Each stage has a sentinel file at its expected output path. Re-running with the same `--sample-id` and `--output-root` picks up where the last run left off. `--force-rerun` overrides.

### Rule-set flags

the rules (internal issue review) are configurable:

- `--rule1-threshold 100` — primary count above which the celltype is primary-only (rule 1).
- `--donor-borrow-cap 100` — max donor cells to borrow per missing celltype (rule 2).
- `--primary-only-celltype tumor --primary-only-celltype liver` — celltypes always sourced exclusively from the primary. Pass repeatedly. Guards summary Caveat §1.
- `--random-state 1` — seed for `donor_balanced_sample_by_reference` in the balanced (intermediate) case.
- `--celltype-col Final_level1_celltype_annotation` — per-cell celltype column in `.obs`. Same column across primary + all donors.

### Stage-C (R side) flags

- `--rscript-bin Rscript` — path to the R interpreter. Default `Rscript` (PATH lookup); pass an absolute path to skip PATH.
- `--min-umi 10` — `spacexr::Reference` `min_UMI` (drops cells below this).
- `--require-int true` — `spacexr::Reference` `require_int`.
- `--export-layer counts` — which anndata layer to export as the mtx (default: `counts`, fallback: `raw_count`, final fallback: `.X`).

## User YAML

If you have a lot of overrides, park them in a user YAML and load it in one flag:

```yaml
# my_bladder_run.yaml
donor_h5ads:
  - /data/SAMPLE2/scRNA/SAMPLE2_preprocessed_scRNA.h5ad
  - /data/SAMPLE3/scRNA/SAMPLE3_preprocessed_scRNA.h5ad
tumor_type: Bladder
census:
  rule1_threshold: 100          # the user rule 1
  donor_borrow_cap: 100         # the user rule 2
  primary_only_celltypes: [tumor, liver]
```

```bash
ref-build run \
    --sample-id  SAMPLE1 \
    --primary-h5ad /data/SAMPLE1/scRNA/SAMPLE1_preprocessed_scRNA.h5ad \
    --celltype-marker-json /data/markers/markers_global_no_tumor_cells_level1.json \
    --output-root /data/ref_build_runs \
    --config my_bladder_run.yaml
```

CLI flags always win over the user YAML.

## Reading the census.csv

`<output_root>/<sample_id>/census/census.csv` is the operator-facing audit trail. Columns:

| Column | Meaning |
|---|---|
| `celltype` | The expected celltype (one row per key in the marker JSON). |
| `primary_count` | Number of cells the primary has for this celltype. |
| `per_donor_counts` | JSON dict `{donor_id: count}` — how many cells each donor has. |
| `decision` | One of `primary_only`, `balanced`, `borrowed`, `missing_no_donor`. |
| `source_samples` | Which samples contribute cells (`<primary>` + donor ids, semicolon-separated). |
| `note` | One-line rationale for the decision (why rule 1 fired, cap value, etc.). |

If a celltype ends up as `missing_no_donor`, the pipeline logs a WARN line to stdout AND records the row in the CSV — you can act on either.

## Output layout

See the [Outputs](../README.md#outputs) section of the README. Every run also writes a `config.yaml` at the run root, capturing the exact merged config that the run used — useful when reproducing a specific result months later.
