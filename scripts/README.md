# workflow-driver

Top-level sbatch driver + per-stage stubs for the user's
`xenium-preprocess` → `ref-build` → `rctd-split` spatial-genomics
workflow.

Specification: (internal issue review) comment
(internal issue review).

## Files

- `submit_workflow.sh` — the driver. Parses CLI args, enforces the
  `--force` guard, submits `xenium-preprocess`, then chains `ref-build`
  (`afterok:JOB1`) and `rctd-split` (`afterok:JOB3`) with a shared
  `RUN_ID` threaded via `sbatch --export=ALL,…,RUN_ID=…`.
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
| `--flex-h5ad <path>` | Flex scRNA h5ad — recorded verbatim under `ref_build.flex_h5ad_path`, no copy, no symlink (required). |
| `--celltype-marker-json <path>` | Marker-gene JSON declaring the expected celltype set — threaded to `ref-build run --celltype-marker-json` (required by `ref-build`). |
| `--output-root <dir>` | Root output directory. Required unless the `OUTPUT_ROOT` env var is set. |
| `--run-id <id>` | Explicit run identifier. Precedence: `--run-id` > `$RUN_ID` env > JOB1's `SLURM_JOB_ID`. |
| `--proseg-dir <dir>` | Override `xenium-preprocess` config: proseg output dir. |
| `--xenium-cells <path>` | Override `xenium-preprocess` config: xenium `cells.parquet`. |
| `--xenium-ranger-dir <dir>` | Override `xenium-preprocess` config: xenium-ranger dir. |
| `--donor-h5ad <path>` | `ref-build` donor scRNA h5ad. **Repeatable** for multiple donors. Threaded to `ref-build run --donor-h5ad …`. Default: empty (no donor supplementation). |
| `--fallback-donor-h5ad <path>` | `ref-build` Rule-5 fallback donor h5ad. **Repeatable**. Threaded to `ref-build run --fallback-donor-h5ad …`. Default: empty (Rule 5 skipped). |
| `--start-step <name>` | Skip earlier steps and start submission at the named step. Accepts `xenium-preprocess` (default; full chain), `ref-build` (submits `ref-build` then `rctd-split`), or `rctd-split` (submits `rctd-split` only). Numeric aliases `1`/`3`/`4` are accepted for one release for backwards compat. Requires `--run-id`. Implies `--reuse-run-dir`. Fails loud if the resumed step's prior inputs are missing. |
| `--reuse-run-dir` | Proceed even if `<run-dir>` exists, **keeping its contents**. Each stage overwrites the files it writes; other files preserved. Recommended for resume flows. Mutually exclusive with `--force`. |
| `--force` | `rm -rf <run-dir>` then proceed. Destructive; intended for from-scratch re-run under an already-used run-id. Rejected with `--start-step ref-build`/`rctd-split` (would wipe the prerequisites). |
| `--env-name <name>` | Convenience: set the conda env name for **all three** stages at once. Equivalent to passing `--xenium-preprocess-env`/`--ref-build-env`/`--rctd-split-env` with the same value. Default: each per-stage sbatch script's own fallback (`xenium`). |
| `--xenium-preprocess-env <name>` | Conda env for `submit_step1.sbatch` (threaded via `XENIUM_PREPROCESS_ENV`). Default: `xenium`. |
| `--ref-build-env <name>` | Conda env for `submit_step3.sbatch` (threaded via `REF_BUILD_ENV`). Default: `xenium`. |
| `--rctd-split-env <name>` | Conda env for `submit_step4.sbatch` (threaded via `RCTD_SPLIT_ENV`). Default: `xenium`. |
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

### Per-step conda envs

By default every step activates the `xenium` conda env. Override
per-step (e.g. one env for `xenium-preprocess`/`ref-build`, another
for `rctd-split`):

```bash
./submit_workflow.sh --sample-id SAMPLE1 \
                     --flex-h5ad          /data/SAMPLE1/scRNA/SAMPLE1_flex.h5ad \
                     --celltype-marker-json /data/markers/markers.json \
                     --xenium-preprocess-env xenium-test \
                     --ref-build-env         xenium-test \
                     --rctd-split-env        xenium-alt
```

For the common case of "one shared env for the whole workflow":

```bash
./submit_workflow.sh --sample-id SAMPLE1 \
                     --flex-h5ad          /data/SAMPLE1/scRNA/SAMPLE1_flex.h5ad \
                     --celltype-marker-json /data/markers/markers.json \
                     --env-name xenium-test
```

`--env-name` sets all three at once; per-step flags **later on the
CLI** override the bulk assignment (`--env-name shared
--ref-build-env alt` → `{shared, alt, shared}`).

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
└── logs/{xenium-preprocess,ref-build,rctd-split}.log
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

`--start-step <name>` skips earlier steps of the chain when their
outputs already exist in the run folder — the user's case from
(internal issue review).

Resume from `ref-build` (`xenium-preprocess` already produced its
outputs, keep them):

```bash
./submit_workflow.sh \
    --sample-id           SAMPLE1 \
    --output-root         /data/workflow_runs_test \
    --run-id              demo_v1 \
    --start-step          ref-build \
    --reuse-run-dir \
    --flex-h5ad           /data/SAMPLE1/scRNA/SAMPLE1_flex.h5ad \
    --celltype-marker-json /data/markers/markers.json
```

Preflight checks (all fail-loud, before any sbatch submission):

- `--start-step ref-build` requires `spatial_adata/<S>_proseg_raw.h5ad`,
  `spatial_adata/<S>_xenium_ranger.h5ad`, and `rctd/<S>_test_object.rds`.
- `--start-step rctd-split` additionally requires `rctd/<S>_reference.rds`.
- `--start-step` past `xenium-preprocess` requires `--run-id` (or
  `$RUN_ID` env).
- `--start-step ref-build`/`rctd-split` cannot combine with `--force`
  (would delete the prerequisites).

`--start-step ref-build`/`rctd-split` implies `--reuse-run-dir`;
passing it explicitly is fine (redundant but harmless).

## Existing-folder policy (`--force` vs `--reuse-run-dir`)

Two edge cases per comment (internal) §2 overwrite audit require an
explicit user opt-in to fire:

- **A**: user reuses `--run-id` deliberately — driver refuses by default.
- **B**: `xenium-preprocess` re-run after `rctd-split` already wrote
  back — same refusal.

Three ways to proceed when `<run-dir>` exists:

- Default: **refuse** with exit 3.
- `--reuse-run-dir`: **keep** the folder. Each stage overwrites the files
  it writes; other files preserved. Recommended for resumes.
- `--force`: **`rm -rf` then proceed**. Destructive; use for
  from-scratch re-run under an already-used run-id. Mutually exclusive
  with `--reuse-run-dir` and rejected with `--start-step ref-build`/`rctd-split`.

The guard only fires when `RUN_ID` is bound BEFORE submission (via
`--run-id` or `$RUN_ID`). A fresh submission (no `RUN_ID` set) uses
JOB1's `SLURM_JOB_ID`, which by construction names a never-before-seen
folder.
