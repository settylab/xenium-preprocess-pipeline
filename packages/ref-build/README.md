# ref-build

**Function.** Build an RCTD-compatible cell-type reference from Flex scRNA-seq
datasets — a primary sample plus an optional donor pool — with a marker-JSON-driven
completeness policy so under-represented celltypes get borrowed from donors.

Part of the [`xenium-preprocess-pipeline`](../../README.md) chain. Runs after
`xenium-preprocess` and before `rctd-split` under a shared `--run-id`.

## Sub-stages (run in order)

| Stage | Purpose |
|---|---|
| `load_primary_and_donors` | Primary + donor h5ads → concatenated AnnData (`anndata.concat(join="outer")`) with a `sample_ID` provenance column. |
| `census`                  | Read the celltype-marker JSON, tally per-sample per-celltype counts, decide a per-celltype migration policy (`primary_only` / `balanced` / `borrowed` / `fallback_borrow` / `missing_no_donor`). Writes `census.csv`. |
| `assemble`                | Apply the census policy → merged reference AnnData with raw integer counts preserved (RCTD's `require_int=TRUE`). |
| `export_mtx`              | Merged AnnData → 10X-style bundle (`counts.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`, `metadata.csv`) that Seurat + spacexr can read. |
| `rctd_reference_build`    | R shell-out — Seurat reads the bundle, `spacexr::Reference(min_UMI=10, require_int=TRUE)`, `saveRDS` to `<sample>_reference.rds`. |

Migration policy (from `census`):

| Primary count | Donor has it? | Fallback has it? | Decision |
|---|---|---|---|
| > `rule1_threshold` (100) | ignored | ignored | `primary_only` |
| `[1, rule1_threshold]`    | any     | ignored | `balanced` (primary + donor supplementation) |
| `0`                       | yes     | ignored | `borrowed` (donors that have it, capped by `donor_borrow_cap`) |
| `0`                       | no      | yes     | `fallback_borrow` (fallback donor(s), capped by `fallback_borrow_cap`) |
| `0`                       | no      | no      | `missing_no_donor` (warn + skip) |

`primary_only_celltypes` (default `[tumor, liver]`) always forces `primary_only`,
regardless of count. Census variants (`original` / `balanced`) control the
intermediate-case donor-supplement math — see `docs/methods.md`.

Every stage is sentinel-gated (idempotent resume); pass `--force-rerun` to redo.

## Quickstart

```bash
micromamba activate xenium
ml fhR/4.4.1-foss-2023b     # for rctd_reference_build

ref-build run \
    --sample-id            SAMPLE1 \
    --run-id               demo_v1 \
    --output-root          /data/workflow_runs \
    --primary-h5ad         /data/SAMPLE1/scRNA/SAMPLE1_preprocessed_scRNA.h5ad \
    --donor-h5ad           /data/SAMPLE2/scRNA/SAMPLE2_preprocessed_scRNA.h5ad \
    --celltype-marker-json /data/markers/markers.json
```

Run a subset with `--stages`:

```bash
ref-build run --sample-id SAMPLE1 --run-id demo_v1 --output-root /data/workflow_runs \
    --primary-h5ad /… --celltype-marker-json /… \
    --stages census assemble
```

## Outputs

Under `<output_root>/<sample_id>/<sample_id>_<run_id>/`:

```
intermediate/loaded/concat.h5ad                        load_primary_and_donors
intermediate/census/census.csv                         census — the audit trail
intermediate/assembled/reference.h5ad                  assemble
intermediate/mtx_bundle/<sample>_{counts,features,barcodes}.*   export_mtx
rctd/<sample>_reference.rds                            rctd_reference_build — RCTD-ready reference
config.yaml                                            resolved config (this stage writes the ref_build: key)
logs/slurm-<jobid>-ref-build.log                       slurm-captured stdout+stderr for the stage
```

## Configuration

Every knob lives in `config/default.yaml`. Override precedence:
**CLI flag > user YAML (`--config …`) > default YAML**.

Required-at-runtime: `sample_id`, `run_id`, `output_root`, `primary_h5ad`,
`celltype_marker_json`. `primary_h5ad` and every `donor_h5ad` must have integer
counts in `layers["counts"]` (or `layers["raw_count"]` / `.X`) and a `.obs`
celltype column (default `celltypes`; override via `--celltype-col`).

Key optional knobs (defaults in parentheses):

| Key | Default | Meaning |
|---|---|---|
| `donor_borrow_cap`         | `100`                          | Max donor cells to borrow per missing celltype. |
| `rule1_threshold`          | `100`                          | Primary count above which donors are ignored. |
| `primary_only_celltypes`   | `[tumor, liver]`               | Always primary-only regardless of census. |
| `census.variant`           | `original`                     | `original` (per-donor cap) or `balanced` (total donor = primary_count). |
| `random_seed`              | `42`                           | Seed for donor sampling. `random_state` accepted as a deprecated alias. |
| `min_UMI` / `require_int`  | `10` / `true`                  | Passed to `spacexr::Reference`. |

Full CLI reference: `ref-build run --help`.

## Installation

```bash
../../scripts/create-env.sh -n xenium -f environment.yml   # once
micromamba activate xenium
uv pip install -e .
```

For the R stage, load an R module with Seurat + spacexr + Matrix — on Fred Hutch:

```bash
ml fhR/4.4.1-foss-2023b
```

`spacexr` installs from GitHub, not the module or CRAN — see
[`docs/install.md`](docs/install.md) § R side. No GitHub credentials
required (it's a public repo), but see that doc for GitHub's anonymous
rate-limit caveat.

## Citation

`CITATION.cff`. Upstream tools to cite alongside:

- **scanpy** — Wolf, F.A. *et al.* Genome Biology 19, 15 (2018).
- **anndata** — Virshup *et al.* bioRxiv 2021.12.16.473007 (2021).
- **Seurat** — Hao *et al.* Cell 184, 3573–3587.e29 (2021).
- **RCTD / spacexr** — Cable *et al.* Nat. Biotechnol. 40, 517–526 (2022).

## License

MIT — see [`LICENSE`](LICENSE).
