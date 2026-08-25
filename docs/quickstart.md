# Quickstart

Once installed (see [`installation.md`](installation.md)), the fastest
path to a first result is one command:

```bash
./scripts/submit_workflow.sh \
    --sample-id             SAMPLE1 \
    --output-root           /data/workflow_runs \
    --run-id                demo_v1 \
    --flex-h5ad             /data/SAMPLE1/scRNA/SAMPLE1_flex.h5ad \
    --celltype-marker-json  /data/markers/markers.json \
    --proseg-dir            /data/SAMPLE1/proseg \
    --xenium-ranger-dir     /data/SAMPLE1/xenium_ranger \
    --donor-h5ad            /data/SAMPLE2/scRNA/SAMPLE2_flex.h5ad \
    --donor-h5ad            /data/SAMPLE3/scRNA/SAMPLE3_flex.h5ad
```

This submits three Slurm jobs chained with `--dependency=afterok:` for
`xenium-preprocess`, `ref-build`, and `rctd-split`. The driver prints
the job ids and the dependency chain, and writes an authoritative
record to `<run-dir>/logs/workflow-submit.log`.

## Inputs

| Flag                       | What it is                                           | Required |
|----------------------------|------------------------------------------------------|----------|
| `--sample-id`              | Short sample identifier, used in all output filenames | Yes      |
| `--output-root`            | Directory under which run folders are created         | Yes      |
| `--run-id`                 | Identifier for this specific run of the sample        | Yes      |
| `--flex-h5ad`              | 10x Flex scRNA h5ad for the reference build           | Yes      |
| `--celltype-marker-json`   | Marker gene JSON for the reference build              | Yes      |
| `--proseg-dir`             | Directory of proseg output for `xenium-preprocess`    | Yes      |
| `--xenium-ranger-dir`      | Xenium Ranger output bundle for `xenium-preprocess`   | Yes      |
| `--donor-h5ad` (repeatable)| Donor h5ads mixed into the reference                  | Optional |
| `--fallback-donor-h5ad`    | Fallback donor if a specific donor is missing         | Optional |
| `--start-step`             | Resume at `ref-build` or `rctd-split` in an existing run folder | Optional |
| `--reuse-run-dir`          | Proceed against an existing run folder                | Optional |
| `--force`                  | `rm -rf` the run folder and start from `xenium-preprocess` | Optional |
| `--max-cores`              | Cap parallelism for `rctd-split`                      | Optional |
| `--dry-run`                | Print sbatch calls without submitting                 | Optional |

Full flag reference: `./scripts/submit_workflow.sh --help`, or see the
per-package `--help` output for the CLI flags of an individual step.

## Outputs

All three steps write into one run folder:

```
<output-root>/<sample-id>/<sample-id>_<run-id>/
├── spatial_adata/
├── rctd/
├── config.yaml
└── logs/
```

See the [README](../README.md#outputs) for the full tree.

## Resume a run

If `xenium-preprocess` already ran, resume at `ref-build`:

```bash
./scripts/submit_workflow.sh \
    --sample-id             SAMPLE1 \
    --output-root           /data/workflow_runs \
    --run-id                demo_v1 \
    --flex-h5ad             /data/SAMPLE1/scRNA/SAMPLE1_flex.h5ad \
    --celltype-marker-json  /data/markers/markers.json \
    --start-step            ref-build
```

`--start-step` past `xenium-preprocess` implies `--reuse-run-dir`.
