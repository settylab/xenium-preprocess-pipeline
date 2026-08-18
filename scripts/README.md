# workflow-driver

Top-level sbatch driver + step-N stubs for the user's step-1 → step-3 →
step-4 spatial-genomics workflow.

Specification: (internal issue review) comment
(internal issue review).

## Files

- `submit_workflow.sh` — the driver. Parses CLI args, enforces the
  `--force` guard, submits step 1, then chains step 3 (`afterok:JOB1`)
  and step 4 (`afterok:JOB3`) with a shared `RUN_ID` threaded via
  `sbatch --export=ALL,…,RUN_ID=…`.
- `submit_step1.sbatch` — invokes `xenium-preprocess run`.
- `submit_step3.sbatch` — invokes `ref-build run --flex-h5ad …`.
- `submit_step4.sbatch` — invokes `rctd-split run`.
- `tests/test_submit_workflow.py` — integration test (18 cases) using
  mock `sbatch` + mock package CLIs on `PATH`.

## Usage

```bash
./submit_workflow.sh --sample-id SAMPLE1 \
                     --flex-h5ad /data/SAMPLE1/scRNA/SAMPLE1_flex.h5ad \
                     --celltype-marker-json /data/markers/markers.json
```

Options:

| Flag | Purpose |
|---|---|
| `--sample-id <S>` | Sample id (required). |
| `--flex-h5ad <path>` | Flex scRNA h5ad — recorded verbatim under `step3.flex_h5ad_path`, no copy, no symlink (required). |
| `--celltype-marker-json <path>` | Marker-gene JSON declaring the expected celltype set — threaded to `ref-build run --celltype-marker-json` (required by step 3). |
| `--output-root <dir>` | Root output directory. Required unless the `OUTPUT_ROOT` env var is set. |
| `--run-id <id>` | Explicit run identifier. Precedence: `--run-id` > `$RUN_ID` env > JOB1's `SLURM_JOB_ID`. |
| `--proseg-dir <dir>` | Override step-1 config: proseg output dir. |
| `--xenium-cells <path>` | Override step-1 config: xenium `cells.parquet`. |
| `--xenium-ranger-dir <dir>` | Override step-1 config: xenium-ranger dir. |
| `--donor-h5ad <path>` | Step-3 donor scRNA h5ad. **Repeatable** for multiple donors. Threaded to `ref-build run --donor-h5ad …`. Default: empty (no donor supplementation). |
| `--fallback-donor-h5ad <path>` | Step-3 Rule-5 fallback donor h5ad. **Repeatable**. Threaded to `ref-build run --fallback-donor-h5ad …`. Default: empty (Rule 5 skipped). |
| `--start-step <1\|3\|4>` | Skip earlier steps and start submission at step N. Default `1` (full chain). `3` submits step 3 (no dependency) then step 4; `4` submits step 4 only. Requires `--run-id`. Implies `--reuse-run-dir`. Fails loud if the resumed step's prior inputs are missing. |
| `--reuse-run-dir` | Proceed even if `<run-dir>` exists, **keeping its contents**. Each step overwrites the files it writes; other files preserved. Recommended for resume flows. Mutually exclusive with `--force`. |
| `--force` | `rm -rf <run-dir>` then proceed. Destructive; intended for from-scratch re-run under an already-used run-id. Rejected with `--start-step 3/4` (would wipe the prerequisites). |
| `--dry-run` | Print sbatch commands without submitting. |

### Multi-donor invocation (donors + fallback)

```bash
./submit_workflow.sh --sample-id SAMPLE1 \
                     --flex-h5ad          /data/SAMPLE1/scRNA/SAMPLE1_flex.h5ad \
                     --celltype-marker-json /data/markers/markers.json \
                     --donor-h5ad         /data/SAMPLE2/scRNA/SAMPLE2_flex.h5ad \
                     --donor-h5ad         /data/SAMPLE3/scRNA/SAMPLE3_flex.h5ad \
                     --fallback-donor-h5ad /data/fallback/fallback_donor.h5ad
```

Donor / fallback paths change per run, so they live on the CLI, not
in `config/default.yaml` (per the user's comment
(internal issue review).
Both flags are safe to omit — `ref-build`'s Rule 5 fallback is
skipped when no fallback donors are supplied and no supplementation
happens when no donors are supplied.

## Output layout (locked)

```
<output_root>/<S>/<S>_<run_id>/
├── spatial_adata/
│   ├── <S>_xenium_ranger.h5ad
│   ├── <S>_proseg_raw.h5ad
│   └── <S>_proseg_purified.h5ad
├── rctd/
│   ├── <S>_reference_post_rules.h5ad
│   ├── <S>_test_object.rds
│   ├── <S>_reference.rds
│   └── <S>_rctd_results.rds
├── config.yaml
└── logs/{step1,step3,step4}.log
```

## Tests

```bash
uv pip install pytest
python -m pytest scripts/tests/ -v
```

Expected: 36 passed (9 pre-existing failures in the mocked-chain tests
require a real `micromamba`/`conda` `xenium` env — unrelated to the
driver).

## Resume flows (`--start-step` + `--reuse-run-dir`)

`--start-step <N>` skips earlier steps of the chain when their outputs
already exist in the run folder — the user's case from
(internal issue review).

Resume from step 3 (step 1 already produced its outputs, keep them):

```bash
./submit_workflow.sh \
    --sample-id           SAMPLE1 \
    --output-root         /data/workflow_runs_test \
    --run-id              demo_v1 \
    --start-step          3 \
    --reuse-run-dir \
    --flex-h5ad           /data/SAMPLE1/scRNA/SAMPLE1_flex.h5ad \
    --celltype-marker-json /data/markers/markers.json
```

Preflight checks (all fail-loud, before any sbatch submission):

- `--start-step 3` requires `spatial_adata/<S>_proseg_raw.h5ad`,
  `spatial_adata/<S>_xenium_ranger.h5ad`, and `rctd/<S>_test_object.rds`.
- `--start-step 4` additionally requires `rctd/<S>_reference.rds`.
- `--start-step > 1` requires `--run-id` (or `$RUN_ID` env).
- `--start-step 3/4` cannot combine with `--force` (would delete the
  prerequisites).

`--start-step 3/4` implies `--reuse-run-dir`; passing it explicitly is
fine (redundant but harmless).

## Existing-folder policy (`--force` vs `--reuse-run-dir`)

Two edge cases per comment (internal) §2 overwrite audit require an
explicit user opt-in to fire:

- **A**: user reuses `--run-id` deliberately — driver refuses by default.
- **B**: step 1 re-run after step 4 already wrote back — same refusal.

Three ways to proceed when `<run-dir>` exists:

- Default: **refuse** with exit 3.
- `--reuse-run-dir`: **keep** the folder. Each step overwrites the files
  it writes; other files preserved. Recommended for resumes.
- `--force`: **`rm -rf` then proceed**. Destructive; use for
  from-scratch re-run under an already-used run-id. Mutually exclusive
  with `--reuse-run-dir` and rejected with `--start-step 3/4`.

The guard only fires when `RUN_ID` is bound BEFORE submission (via
`--run-id` or `$RUN_ID`). A fresh submission (no `RUN_ID` set) uses
JOB1's `SLURM_JOB_ID`, which by construction names a never-before-seen
folder.
