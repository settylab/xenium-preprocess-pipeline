#!/usr/bin/env bash
# submit_workflow.sh — top-level driver for the step-1 → step-3 → step-4
# spatial-genomics workflow (settylab/TracyY123-nexus#26 comment 5251080220).
#
# Submits up to three sbatch jobs and chains them with --dependency=afterok:
#     JOB1 (step 1, xenium-preprocess)
#       └── JOB3 (step 3, ref-build)          afterok:JOB1
#             └── JOB4 (step 4, rctd-split)   afterok:JOB3
#
# --start-step lets an operator resume the chain when earlier steps have
# already produced their outputs in the run folder (settylab/TracyY123-nexus#26
# comment 5260289249). --start-step 3 skips JOB1 (step 3 submits with no
# dependency); --start-step 4 skips JOB1 + JOB3.
#
# All submitted jobs share a single RUN_ID, propagated via --export=ALL,RUN_ID=…
# so outputs land in the run-scoped folder:
#     <output_root>/<sample_id>/<sample_id>_<run_id>/
#         ├── spatial_adata/   ← step 1 + step 4 augmentations
#         ├── rctd/            ← step 1 + step 3 + step 4
#         ├── config.yaml  ← merged across steps
#         └── logs/            ← {xenium-preprocess,ref-build,rctd-split}.log
#
# RUN_ID precedence:
#     --run-id <id>        > $RUN_ID env             > JOB1's SLURM_JOB_ID
# When neither --run-id nor $RUN_ID is set, the driver submits step 1
# first (with no explicit RUN_ID) and uses its --parsable job id as the
# shared RUN_ID for steps 3 and 4. Each step-N sbatch script then defaults
# `RUN_ID=${RUN_ID:-$SLURM_JOB_ID}` internally, so JOB1 always names its own
# folder correctly.
#
# --start-step > 1 REQUIRES an explicit RUN_ID (via --run-id or $RUN_ID) —
# there's no way to resume into an existing run folder without knowing which
# run to resume.
#
# Existing-folder policy:
#     Default: if <output_root>/<S>/<S>_<run_id>/ already exists, REFUSE
#     with exit 3 (safety catch for accidental re-use of a bound run-id).
#     --reuse-run-dir: proceed and KEEP the folder. Each step overwrites the
#         specific files it writes; other files in the folder are preserved
#         (e.g. step-1 outputs survive an isolated step-3+4 resume). This is
#         the resume path — implied by --start-step 3 / --start-step 4.
#     --force: `rm -rf` the run folder, then proceed. Destructive: intended
#         for a from-scratch re-run under an already-used run-id. Mutually
#         exclusive with --reuse-run-dir. Rejected when combined with
#         --start-step 3/4 (would wipe the prerequisites we're resuming from).
#
# Usage:
#     ./submit_workflow.sh --sample-id MH10 --flex-h5ad /path/flex.h5ad \
#                          --celltype-marker-json /path/markers.json
#     # Resume from step 3 (step 1 already ran, keep its outputs):
#     ./submit_workflow.sh --sample-id MH10 --flex-h5ad /path/flex.h5ad \
#                          --celltype-marker-json /path/markers.json \
#                          --run-id my_experiment_v2 \
#                          --start-step 3 --reuse-run-dir
#     # From-scratch re-run under an existing run-id (destructive):
#     ./submit_workflow.sh --sample-id MH10 --flex-h5ad /path/flex.h5ad \
#                          --celltype-marker-json /path/markers.json \
#                          --run-id my_experiment_v2 --force
#     RUN_ID=my_id ./submit_workflow.sh --sample-id MH10 \
#                          --flex-h5ad /path/flex.h5ad \
#                          --celltype-marker-json /path/markers.json
#     ./submit_workflow.sh --sample-id MH10 --flex-h5ad /path/flex.h5ad \
#                          --celltype-marker-json /path/markers.json \
#                          --dry-run       # print sbatch commands, don't submit
#
# --donor-h5ad / --fallback-donor-h5ad:
#     Both flags are OPTIONAL and REPEATABLE. Each occurrence threads a
#     path down to step 3 (ref-build) via numbered env vars:
#       DONOR_H5AD_COUNT=N, DONOR_H5AD_1=…, DONOR_H5AD_2=…, …
#       FALLBACK_H5AD_COUNT=M, FALLBACK_H5AD_1=…, FALLBACK_H5AD_2=…, …
#     submit_step3.sbatch expands them into --donor-h5ad / --fallback-donor-h5ad
#     args on `ref-build run`. Empty (default) is safe: ref-build's Rule 5
#     fallback path is skipped when no fallback donors are passed and no
#     supplementation happens when no donors are passed.
#
# --test-object / --reference-rds
# (settylab/TracyY123-nexus#26 comment 5278875595):
#     Explicit step-4 inputs that bypass the default run-folder
#     auto-discovery (`<run>/rctd/<sample>_test_object.rds` +
#     `<run>/rctd/<sample>_reference.rds`). Use this to mix a
#     test object from one sample with a reference from another
#     (e.g. `--sample-id MH3 --test-object …MH3….rds
#           --reference-rds …MH2_scRNA_ref.rds`).
#     Both flags are OPTIONAL, but MUTUAL — passing only one is
#     a hard error. When both are set, the driver threads the
#     paths to submit_step4.sbatch (STEP4_TEST_OBJECT /
#     STEP4_REFERENCE_RDS), which materialises them as
#     `rctd-split run --test-object … --reference-rds …`; the
#     CLI marks them as `source: config` and does NOT touch the
#     run-folder layout for these two artifacts.
#
# Per-step parameter exposure
# (settylab/TracyY123-nexus#26 comment 5277569725):
#     Every step's config surface is regulable from this driver via three
#     layers, applied in this precedence (last wins):
#         (1) step's own config/default.yaml
#         (2) --stepN-config <path>              (full user YAML for step N)
#         (3) --override stepN.<dotted.key>=<yaml-val>   (repeatable)
#         (4) --stepN-<param> <value>            (named driver flag)
#     Named flags cover the ~15 params Tracy has historically tuned. The
#     --override form is a catch-all for any nested config key not covered
#     by a named flag; values are YAML-parsed (so lists / bools / ints
#     round-trip via `[0.5, 0.7]` / `true` / `42`).
#     Threading: each step gets a base64-encoded YAML env var
#     (STEPN_OVERRIDES_B64) built from layers 2 + 3, decoded in the sbatch
#     script and passed via `--config`. Named flags (layer 4) are threaded
#     as env vars (STEPN_<PARAM>=<value>) and materialized as `--flag value`
#     on the step CLI, whose own precedence order lets them override the
#     yaml naturally. Base64 encoding sidesteps the comma-in-value trap
#     with slurm's `--export=ALL,K=V,K=V` payload format.

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults + arg parsing
# ---------------------------------------------------------------------------

# Script dir — location of submit_stepN.sbatch stubs.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

DEFAULT_OUTPUT_ROOT=${OUTPUT_ROOT:-/fh/fast/setty_m/user/ryang/workflow_runs}

SAMPLE_ID=""
FLEX_H5AD=""
CELLTYPE_MARKER_JSON=""
OUTPUT_ROOT="$DEFAULT_OUTPUT_ROOT"
RUN_ID_OVERRIDE="${RUN_ID:-}"
FORCE=0
REUSE_RUN_DIR=0
DRY_RUN=0
START_STEP=1
# Optional extra step-1 inputs (all resolved by config/default.yaml if
# omitted, but Tracy's normal flow needs them named).
PROSEG_DIR=""
XENIUM_CELLS=""
XENIUM_RANGER_DIR=""
# Optional step-3 donor pool (both repeatable, both default empty).
DONOR_H5ADS=()
FALLBACK_H5ADS=()
# Optional step-3 override: which .obs column ref-build reads for the
# per-cell celltype label (threaded to `ref-build run --celltype-col`).
# Empty ⇒ ref-build's config default (Final_level1_celltype_annotation).
CELLTYPE_COL_FOR_REF_BUILD=""
# Optional step-4 override: RCTD parallelism. Threaded to submit_step4.sbatch
# as MAX_CORES → `rctd-split run --max-cores N`. Also honors $MAX_CORES env.
# Empty ⇒ submit_step4.sbatch's own default (12). NOTE: only meaningful up to
# step-4's sbatch alloc (--cpus-per-task=16); larger values will oversubscribe.
MAX_CORES_OVERRIDE="${MAX_CORES:-}"
# Optional step-4 explicit inputs — bypass the run-folder layout auto-discovery
# of test_object.rds + reference.rds. Both are OPTIONAL individually but
# MUTUAL: passing only one is rejected below. When both are set, they're
# threaded as STEP4_TEST_OBJECT / STEP4_REFERENCE_RDS and materialised as
# `rctd-split run --test-object … --reference-rds …`.
TEST_OBJECT=""
REFERENCE_RDS=""

# ---------------------------------------------------------------------------
# Per-step named parameter overrides. Empty ⇒ step CLI's config default.
# Every entry here maps 1:1 onto a flag the step's own CLI already accepts,
# so we can thread as env var → materialize as `--flag value` in the sbatch
# script. Precedence: named flag wins over --override and --stepN-config
# (matches each step CLI's own precedence).
# ---------------------------------------------------------------------------
# Step 1 (xenium-preprocess)
STEP1_X_SOURCE=""                # --x-source           (maxpost_counts|expected_counts)
STEP1_QC_MIN_COUNTS_CELL=""      # --qc-min-counts-cell (int)
STEP1_GEX_ONLY=""                # --gex-only           (bool)
STEP1_FORCE_RERUN=0              # --force-rerun        (flag)
# Step 3 (ref-build)
STEP3_DONOR_BORROW_CAP=""        # --donor-borrow-cap        (int)
STEP3_CELL_MIN_INSTANCE=""       # --cell-min-instance       (int)
STEP3_MIN_UMI=""                 # --min-umi                 (int)
STEP3_RANDOM_SEED=""             # --random-seed             (int)
STEP3_CELLTYPE_TARGET_LIST=""    # --celltype-target-list    (path)
# Step 4 (rctd-split)
STEP4_UMI_MIN=""                 # --umi-min                    (int)
STEP4_COUNTS_MIN=""              # --counts-min                 (int)
STEP4_CELL_MIN_INSTANCE=""       # --cell-min-instance          (int)
STEP4_DOUBLET_MODE=""            # --doublet-mode               (str)
STEP4_POSTPROCESS_MIN_COUNTS=""  # --postprocess-min-counts     (int)
STEP4_KEEP_INTERMEDIATE=0        # --keep-intermediate          (flag)

# ---------------------------------------------------------------------------
# Per-step config file (layer 2) and dotted-key --override list (layer 3).
# Empty by default ⇒ step falls through to its own default.yaml. See the
# module docstring's "Per-step parameter exposure" block for precedence.
# ---------------------------------------------------------------------------
STEP1_CONFIG=""
STEP3_CONFIG=""
STEP4_CONFIG=""
STEP1_OVERRIDES=()
STEP3_OVERRIDES=()
STEP4_OVERRIDES=()

# ---------------------------------------------------------------------------
# --config <path>: driver-level config file (settylab/TracyY123-nexus#26
# comment 5320970944 item 6, design under
# reports/no-shared-folder-config-design_*). One file per run under
# runs/<sample>_<run_id>/config.yaml carries the values the operator would
# otherwise pass as CLI flags to submit_workflow.sh: driver top-level fields
# (sample_id, run_id, output_root, flex_h5ad, celltype_marker_json,
# test_object, reference_rds, proseg_dir, xenium_cells, xenium_ranger_dir,
# celltype_col_for_ref_build, max_cores) plus per-step scalars nested under
# step1: / step3: / step4:.
#
# Precedence is "CLI wins over YAML": we pre-scan argv for --config and
# apply its values as DEFAULTS here, before the arg-parse loop runs; the
# loop then naturally overwrites anything the operator also passed on the
# CLI. This mirrors how DEFAULT_OUTPUT_ROOT etc. seed defaults above.
# ---------------------------------------------------------------------------

WORKFLOW_CONFIG=""
for ((_i=1; _i<=$#; _i++)); do
    if [[ "${!_i}" == "--config" ]]; then
        _j=$((_i+1))
        if [[ $_j -gt $# ]]; then
            echo "error: --config requires a path argument" >&2
            exit 2
        fi
        WORKFLOW_CONFIG="${!_j}"
        break
    fi
done
unset _i _j

if [[ -n "$WORKFLOW_CONFIG" ]]; then
    if [[ ! -f "$WORKFLOW_CONFIG" ]]; then
        echo "error: --config file not found: $WORKFLOW_CONFIG" >&2
        exit 2
    fi
    # Parse the YAML into a bash `eval`-able block. Every emitted line is
    # `NAME=<shlex-quoted-value>`, matching the driver's existing env-var
    # names 1:1. Keys missing from the YAML emit nothing, so the driver's
    # own scalar defaults stay in place. YAML `true`/`false` for the
    # gex-only / *-force-rerun / keep-intermediate flag surfaces round-trip
    # as the strings "1"/"0" the sbatch scripts already look for.
    _CONFIG_ASSIGNS=$(WORKFLOW_CONFIG="$WORKFLOW_CONFIG" python3 - <<'PY'
import os, shlex, sys, yaml

path = os.environ["WORKFLOW_CONFIG"]
with open(path) as f:
    cfg = yaml.safe_load(f) or {}
if not isinstance(cfg, dict):
    sys.exit(f"error: --config {path}: top-level must be a mapping, got {type(cfg).__name__}")

def emit(name, value):
    if value is None:
        return
    print(f"{name}={shlex.quote(str(value))}")

def emit_flag(name, value):
    # Bool-ish YAML values → the "1"/"0" the sbatch scripts already look for.
    if value is None:
        return
    if isinstance(value, bool):
        v = "1" if value else "0"
    else:
        v = str(value)
    print(f"{name}={shlex.quote(v)}")

# Driver top-level scalars. Accept `sample` as alias for `sample_id`
# (matches the memo Tracy uses in issue-thread discussion).
emit("SAMPLE_ID",              cfg.get("sample_id") or cfg.get("sample"))
emit("RUN_ID_OVERRIDE",        cfg.get("run_id"))
emit("OUTPUT_ROOT",            cfg.get("output_root"))
emit("FLEX_H5AD",              cfg.get("flex_h5ad"))
emit("CELLTYPE_MARKER_JSON",   cfg.get("celltype_marker_json"))
emit("TEST_OBJECT",            cfg.get("test_object"))
emit("REFERENCE_RDS",          cfg.get("reference_rds"))
emit("PROSEG_DIR",             cfg.get("proseg_dir"))
emit("XENIUM_CELLS",           cfg.get("xenium_cells"))
emit("XENIUM_RANGER_DIR",      cfg.get("xenium_ranger_dir"))
emit("CELLTYPE_COL_FOR_REF_BUILD", cfg.get("celltype_col_for_ref_build"))
emit("MAX_CORES_OVERRIDE",     cfg.get("max_cores"))

# Per-step scalars. Each of these has a matching --stepN-<param> CLI flag
# in the driver; the YAML key mirrors the flag name with underscores.
step1 = cfg.get("step1") or {}
if isinstance(step1, dict):
    emit(     "STEP1_X_SOURCE",           step1.get("x_source"))
    emit(     "STEP1_QC_MIN_COUNTS_CELL", step1.get("qc_min_counts_cell"))
    emit(     "STEP1_GEX_ONLY",           step1.get("gex_only"))
    emit_flag("STEP1_FORCE_RERUN",        step1.get("force_rerun"))

step3 = cfg.get("step3") or {}
if isinstance(step3, dict):
    emit("STEP3_DONOR_BORROW_CAP",     step3.get("donor_borrow_cap"))
    emit("STEP3_CELL_MIN_INSTANCE",    step3.get("cell_min_instance"))
    emit("STEP3_MIN_UMI",              step3.get("min_umi"))
    emit("STEP3_RANDOM_SEED",          step3.get("random_seed"))
    emit("STEP3_CELLTYPE_TARGET_LIST", step3.get("celltype_target_list"))

step4 = cfg.get("step4") or {}
if isinstance(step4, dict):
    emit(     "STEP4_UMI_MIN",                step4.get("umi_min"))
    emit(     "STEP4_COUNTS_MIN",             step4.get("counts_min"))
    emit(     "STEP4_CELL_MIN_INSTANCE",      step4.get("cell_min_instance"))
    emit(     "STEP4_DOUBLET_MODE",           step4.get("doublet_mode"))
    emit(     "STEP4_POSTPROCESS_MIN_COUNTS", step4.get("postprocess_min_counts"))
    emit_flag("STEP4_KEEP_INTERMEDIATE",      step4.get("keep_intermediate"))
PY
    )
    if [[ -n "$_CONFIG_ASSIGNS" ]]; then
        eval "$_CONFIG_ASSIGNS"
    fi
    unset _CONFIG_ASSIGNS
fi

usage() {
    cat <<'EOF'
Usage: submit_workflow.sh --sample-id <S> --flex-h5ad <path>
                          --celltype-marker-json <path> [OPTIONS]

Required:
  --sample-id <S>              Sample identifier (e.g. MH10).
  --flex-h5ad <path>           Flex scRNA h5ad (recorded verbatim under
                               step3.flex_h5ad_path in config.yaml;
                               no copy / no symlink).
  --celltype-marker-json <p>   Marker-gene JSON declaring the expected
                               celltype set (step-3 required input;
                               threaded to ref-build run --celltype-marker-json).

Optional:
  --config <path>              Driver-level config YAML. Top-level keys
                               `sample_id`, `run_id`, `output_root`,
                               `flex_h5ad`, `celltype_marker_json`,
                               `test_object`, `reference_rds`,
                               `proseg_dir`, `xenium_cells`,
                               `xenium_ranger_dir`,
                               `celltype_col_for_ref_build`, `max_cores`
                               plus per-step scalars nested under
                               `step1:` / `step3:` / `step4:` (each key
                               mirrors the matching --stepN-<param> CLI
                               flag with underscores). CLI flags win over
                               the YAML — the file provides defaults.
                               Recommended location:
                               <output_root>/<sample>/<sample>_<run_id>/config.yaml.
  --output-root <dir>          Root output directory (default: env
                               OUTPUT_ROOT, else
                               /fh/fast/setty_m/user/ryang/workflow_runs).
  --run-id <id>                Explicit run identifier. Precedence:
                               --run-id > $RUN_ID env > JOB1 SLURM_JOB_ID.
                               Required when --start-step > 1.
  --start-step <name>          Skip earlier steps and start submission at
                               the named step. Accepts either the numeric
                               name (1, 3, 4) or the semantic name
                               (xenium-preprocess, ref-build, rctd-split);
                               they alias 1:1. Default xenium-preprocess
                               (full chain). --start-step ref-build submits
                               step 3 (no dep) then step 4; assumes step-1
                               outputs already exist in the run folder.
                               --start-step rctd-split submits only step 4;
                               assumes step-1 + step-3 outputs exist. Both
                               imply --reuse-run-dir. Fails loud if the
                               required prior outputs are missing. The
                               numeric aliases are kept for one release for
                               backwards compat.
  --reuse-run-dir              Proceed even if <run-dir> already exists,
                               KEEPING its contents. Each step overwrites
                               the files it writes; other files preserved.
                               Recommended for resume flows. Mutually
                               exclusive with --force.
  --force                      rm -rf <run-dir> then proceed. Destructive;
                               intended for from-scratch re-run under an
                               already-used run-id. Rejected with
                               --start-step 3/4. Mutually exclusive with
                               --reuse-run-dir.
  --proseg-dir <dir>           Overrides step-1 config: proseg output dir.
  --xenium-cells <path>        Overrides step-1 config: xenium cells.parquet.
  --xenium-ranger-dir <dir>    Overrides step-1 config: xenium-ranger dir.
  --donor-h5ad <path>          Step-3 donor scRNA h5ad. REPEATABLE for
                               multiple donors. Default: empty (no
                               donor supplementation).
  --fallback-donor-h5ad <path> Step-3 Rule-5 fallback donor h5ad.
                               REPEATABLE. Default: empty (Rule 5
                               skipped).
  --celltype-col-for-ref-build <col-name>
                               Step-3 override: which .obs column
                               ref-build reads for the per-cell
                               celltype label (threaded to
                               `ref-build run --celltype-col`).
                               Default: ref-build's config default
                               (Final_level1_celltype_annotation).
  --max-cores <N>              Step-4 override: RCTD parallelism.
                               Threaded to submit_step4.sbatch as
                               MAX_CORES → `rctd-split run --max-cores N`.
                               Default: submit_step4.sbatch default (12).
                               Only meaningful up to step-4's sbatch
                               alloc (--cpus-per-task=16); larger
                               values oversubscribe the R workers.
  --test-object <path>         Step-4 explicit test_object.rds path.
                               Bypasses the run-folder layout
                               auto-discovery. MUTUAL with
                               --reference-rds — passing only one
                               is a hard error. Use to mix a test
                               object from one sample with a
                               reference from another.
                               When combined with --start-step 4,
                               also bypasses the run-folder
                               existence check (the driver mkdir
                               -p's a fresh <run-dir>/logs/); use
                               this to run step 4 from external
                               rds files without having produced
                               step 1 / step 3 outputs in-tree.
  --reference-rds <path>       Step-4 explicit reference.rds path.
                               Bypasses the run-folder layout
                               auto-discovery. MUTUAL with
                               --test-object (see above). Same
                               --start-step 4 bypass semantics.

Per-step named parameters
(TracyY123-nexus#26 comment 5277569725; every one maps to an existing
step CLI flag; empty ⇒ step CLI's default.yaml value):

  Step 1 (xenium-preprocess):
    --step1-x-source <src>            proseg_to_anndata.x_source
                                      (maxpost_counts|expected_counts).
                                      Default: maxpost_counts.
    --step1-qc-min-counts-cell <N>    qc_filter.min_counts_cell (int).
                                      Default: 10.
    --step1-gex-only <bool>           xenium_ranger_to_anndata.gex_only.
                                      Default: true.
    --step1-force-rerun               Nuke step-1's per-stage sentinels
                                      and re-run every stage.

  Step 3 (ref-build):
    --step3-donor-borrow-cap <N>      census.donor_borrow_cap (int).
                                      Default: 100.
    --step3-cell-min-instance <N>     census.cell_min_instance (int).
                                      Default: 20.
    --step3-min-umi <N>               rctd_reference_build.min_UMI (int).
                                      Default: 10.
    --step3-random-seed <N>           census.random_seed (int).
                                      Default: 42.
    --step3-celltype-target-list <p>  census.celltype_target_list (path).
                                      Default: null (marker JSON keys).

  Step 4 (rctd-split):
    --step4-umi-min <N>               rctd_run.UMI_min (int). Default: 10.
    --step4-counts-min <N>            rctd_run.counts_MIN (int).
                                      Default: 10.
    --step4-cell-min-instance <N>     rctd_run.CELL_MIN_INSTANCE (int).
                                      Default: 20.
    --step4-doublet-mode <mode>       rctd_run.doublet_mode (str).
                                      Default: doublet.
    --step4-postprocess-min-counts <N>
                                      postprocess.qc.min_counts (int).
                                      Default: 50.
    --step4-keep-intermediate         Keep <run_dir>/intermediate/ after
                                      step 4 completes. Default: drop.

Per-step catch-all overrides
(for any nested config key NOT covered by the named flags above):

  --stepN-config <path>        Full YAML for step N (N in 1|3|4). Threaded
                               as the step CLI's --config <path>. See each
                               package's config/default.yaml for keys.
  --override stepN.<key>=<val>
                               Repeatable. `stepN` in {step1, step3, step4};
                               `<key>` is a dotted path into that step's
                               YAML (e.g. `postprocess.leiden.resolutions`);
                               `<val>` is YAML-parsed (`42`, `true`,
                               `[0.5, 0.7]`, `some_string`). Merged on top
                               of --stepN-config (if any); named flags
                               above still win over --override.

  --dry-run                    Print sbatch commands but don't submit.
  -h, --help                   Show this message.

Env vars honored:
  OUTPUT_ROOT, RUN_ID,         Same as their --flags.
  MAX_CORES
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample-id)             SAMPLE_ID="$2"; shift 2 ;;
        --flex-h5ad)             FLEX_H5AD="$2"; shift 2 ;;
        --celltype-marker-json)  CELLTYPE_MARKER_JSON="$2"; shift 2 ;;
        --output-root)           OUTPUT_ROOT="$2"; shift 2 ;;
        --run-id)                RUN_ID_OVERRIDE="$2"; shift 2 ;;
        --start-step)            START_STEP="$2"; shift 2 ;;
        --reuse-run-dir)         REUSE_RUN_DIR=1; shift ;;
        --proseg-dir)            PROSEG_DIR="$2"; shift 2 ;;
        --xenium-cells)          XENIUM_CELLS="$2"; shift 2 ;;
        --xenium-ranger-dir)     XENIUM_RANGER_DIR="$2"; shift 2 ;;
        --donor-h5ad)            DONOR_H5ADS+=("$2"); shift 2 ;;
        --fallback-donor-h5ad)   FALLBACK_H5ADS+=("$2"); shift 2 ;;
        --celltype-col-for-ref-build)
                                 CELLTYPE_COL_FOR_REF_BUILD="$2"; shift 2 ;;
        --max-cores)             MAX_CORES_OVERRIDE="$2"; shift 2 ;;
        --test-object)           TEST_OBJECT="$2"; shift 2 ;;
        --reference-rds)         REFERENCE_RDS="$2"; shift 2 ;;
        # Step 1 named params
        --step1-x-source)                STEP1_X_SOURCE="$2"; shift 2 ;;
        --step1-qc-min-counts-cell)      STEP1_QC_MIN_COUNTS_CELL="$2"; shift 2 ;;
        --step1-gex-only)                STEP1_GEX_ONLY="$2"; shift 2 ;;
        --step1-force-rerun)             STEP1_FORCE_RERUN=1; shift ;;
        # Step 3 named params
        --step3-donor-borrow-cap)        STEP3_DONOR_BORROW_CAP="$2"; shift 2 ;;
        --step3-cell-min-instance)       STEP3_CELL_MIN_INSTANCE="$2"; shift 2 ;;
        --step3-min-umi)                 STEP3_MIN_UMI="$2"; shift 2 ;;
        --step3-random-seed)             STEP3_RANDOM_SEED="$2"; shift 2 ;;
        --step3-celltype-target-list)    STEP3_CELLTYPE_TARGET_LIST="$2"; shift 2 ;;
        # Step 4 named params
        --step4-umi-min)                 STEP4_UMI_MIN="$2"; shift 2 ;;
        --step4-counts-min)              STEP4_COUNTS_MIN="$2"; shift 2 ;;
        --step4-cell-min-instance)       STEP4_CELL_MIN_INSTANCE="$2"; shift 2 ;;
        --step4-doublet-mode)            STEP4_DOUBLET_MODE="$2"; shift 2 ;;
        --step4-postprocess-min-counts)  STEP4_POSTPROCESS_MIN_COUNTS="$2"; shift 2 ;;
        --step4-keep-intermediate)       STEP4_KEEP_INTERMEDIATE=1; shift ;;
        # Driver-level workflow config (already consumed by the pre-scan
        # above; skipped here without an error so the loop stays uniform
        # and later CLI flags — which by policy WIN over YAML — parse
        # normally).
        --config)                shift 2 ;;
        # Per-step config files (yaml passthrough)
        --step1-config)          STEP1_CONFIG="$2"; shift 2 ;;
        --step3-config)          STEP3_CONFIG="$2"; shift 2 ;;
        --step4-config)          STEP4_CONFIG="$2"; shift 2 ;;
        # Per-step --override <stepN.dotted.key=yaml-value>
        --override)
            _ov="$2"
            case "$_ov" in
                step1.*) STEP1_OVERRIDES+=("${_ov#step1.}") ;;
                step3.*) STEP3_OVERRIDES+=("${_ov#step3.}") ;;
                step4.*) STEP4_OVERRIDES+=("${_ov#step4.}") ;;
                *)
                    echo "error: --override must be prefixed with step1./step3./step4. (got: $_ov)" >&2
                    exit 2
                    ;;
            esac
            shift 2
            ;;
        --force)                 FORCE=1; shift ;;
        --dry-run)               DRY_RUN=1; shift ;;
        -h|--help)               usage; exit 0 ;;
        *)                       echo "unknown arg: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -z "$SAMPLE_ID" ]]; then
    echo "error: --sample-id is required" >&2
    usage >&2
    exit 2
fi
if [[ -z "$FLEX_H5AD" ]]; then
    echo "error: --flex-h5ad is required" >&2
    usage >&2
    exit 2
fi
if [[ -z "$CELLTYPE_MARKER_JSON" ]]; then
    echo "error: --celltype-marker-json is required" >&2
    usage >&2
    exit 2
fi

# --test-object / --reference-rds: both-or-neither. Passing only one is a
# hard error — mixing an explicit path with layout-derived discovery for
# its sibling would be silently wrong most of the time (different sample
# → wrong reference / wrong test object).
if [[ -n "$TEST_OBJECT" && -z "$REFERENCE_RDS" ]]; then
    echo "error: --test-object requires --reference-rds (both-or-neither)." >&2
    echo "       Passing one explicit input while auto-discovering the other" >&2
    echo "       from the run folder is silently mismatch-prone. Pass both," >&2
    echo "       or omit both and let the layout-derived defaults apply." >&2
    exit 2
fi
if [[ -n "$REFERENCE_RDS" && -z "$TEST_OBJECT" ]]; then
    echo "error: --reference-rds requires --test-object (both-or-neither)." >&2
    echo "       Passing one explicit input while auto-discovering the other" >&2
    echo "       from the run folder is silently mismatch-prone. Pass both," >&2
    echo "       or omit both and let the layout-derived defaults apply." >&2
    exit 2
fi

# Announce the explicit-rds mode once for the caller (visible in --dry-run
# too), so the log makes it obvious step 4 will bypass the run-folder
# layout auto-discovery for test_object.rds / reference.rds.
if [[ -n "$TEST_OBJECT" ]]; then
    echo "info: --test-object/--reference-rds → using explicit rds paths (skipping auto-discovery)" >&2
    echo "      test-object:   $TEST_OBJECT" >&2
    echo "      reference-rds: $REFERENCE_RDS" >&2
fi

# ---------------------------------------------------------------------------
# --start-step validation + alias normalization. Legal values are the numeric
# step names (1, 3, 4) or their semantic aliases (xenium-preprocess, ref-build,
# rctd-split); each pair aliases 1:1. Semantic aliases are normalized to
# numeric here so the rest of the script keeps its integer-step logic. > 1
# needs a pre-bound RUN_ID (nothing to resume without one) and implies
# --reuse-run-dir (the whole point is to keep the earlier steps' outputs).
# Numeric aliases are kept for one release for backwards compat.
# ---------------------------------------------------------------------------

case "$START_STEP" in
    1|xenium-preprocess) START_STEP=1 ;;
    3|ref-build)         START_STEP=3 ;;
    4|rctd-split)        START_STEP=4 ;;
    *)
        echo "error: --start-step must be one of 1|xenium-preprocess, 3|ref-build, 4|rctd-split (got: $START_STEP)" >&2
        exit 2
        ;;
esac

if [[ "$START_STEP" -ne 1 ]]; then
    if [[ -z "$RUN_ID_OVERRIDE" ]]; then
        echo "error: --start-step $START_STEP requires --run-id (or \$RUN_ID env)." >&2
        echo "       Without a bound run-id there is no run folder to resume from." >&2
        exit 2
    fi
    if [[ "$FORCE" -eq 1 ]]; then
        echo "error: --force is incompatible with --start-step $START_STEP." >&2
        echo "       --force wipes the run folder, which would delete the step-1" >&2
        echo "       (and step-3) outputs that --start-step $START_STEP resumes from." >&2
        exit 2
    fi
    # Resume semantics require reusing the existing folder — set the flag
    # implicitly so downstream logic stays uniform.
    REUSE_RUN_DIR=1
fi

if [[ "$FORCE" -eq 1 && "$REUSE_RUN_DIR" -eq 1 ]]; then
    echo "error: --force and --reuse-run-dir are mutually exclusive." >&2
    echo "       --force removes the run folder; --reuse-run-dir preserves it." >&2
    exit 2
fi

# ---------------------------------------------------------------------------
# Existing-folder policy — refuse by default, --force wipes, --reuse-run-dir
# keeps. Only meaningful when RUN_ID is bound BEFORE JOB1 submits (either
# --run-id or $RUN_ID). When RUN_ID = JOB1's SLURM_JOB_ID (Slurm-assigned),
# the folder is by construction fresh.
# ---------------------------------------------------------------------------

RUN_DIR=""
if [[ -n "$RUN_ID_OVERRIDE" ]]; then
    RUN_DIR="$OUTPUT_ROOT/$SAMPLE_ID/${SAMPLE_ID}_${RUN_ID_OVERRIDE}"
    if [[ -e "$RUN_DIR" ]]; then
        if [[ "$FORCE" -eq 1 ]]; then
            echo "warning: --force → removing existing run folder:" >&2
            echo "         $RUN_DIR" >&2
            if [[ "$DRY_RUN" -eq 0 ]]; then
                rm -rf "$RUN_DIR"
            else
                echo "[dry-run] rm -rf $RUN_DIR" >&2
            fi
        elif [[ "$REUSE_RUN_DIR" -eq 1 ]]; then
            echo "info: --reuse-run-dir → keeping existing run folder:" >&2
            echo "      $RUN_DIR" >&2
            echo "      Each submitted step overwrites the files it writes;" >&2
            echo "      other files are preserved in place." >&2
        else
            echo "error: run folder already exists: $RUN_DIR" >&2
            echo "       Pass --reuse-run-dir to keep it (resume-style; new writes overwrite files in place)," >&2
            echo "       --force to wipe it and re-run from step 1," >&2
            echo "       or pick a different --run-id." >&2
            exit 3
        fi
    elif [[ "$START_STEP" -ne 1 ]]; then
        # --start-step 4 with both --test-object and --reference-rds is the
        # one legal way to resume into a non-existent run folder: the two
        # artifacts step 4 would normally auto-discover from the run folder
        # are being replaced by the explicit paths, so there is nothing
        # step-1/step-3 needs to have produced in-tree. The folder + logs/
        # get mkdir -p'd below at LOG_DIR creation time.
        if [[ "$START_STEP" == "4" && -n "$TEST_OBJECT" && -n "$REFERENCE_RDS" ]]; then
            echo "info: --start-step 4 + explicit --test-object/--reference-rds → creating fresh run folder:" >&2
            echo "      $RUN_DIR" >&2
        else
            echo "error: --start-step $START_STEP but the run folder does not exist: $RUN_DIR" >&2
            echo "       There is nothing to resume from — did you mean --start-step 1?" >&2
            exit 3
        fi
    fi
fi

# ---------------------------------------------------------------------------
# Prerequisite preflight for --start-step > 1. We check the artifacts each
# resumed step reads from the run folder. Fail-loud with the missing paths
# listed so an operator can eyeball which upstream step didn't complete.
# ---------------------------------------------------------------------------

_check_prereqs() {
    local -a want=("$@")
    local -a missing=()
    local p
    for p in "${want[@]}"; do
        [[ -e "$p" ]] || missing+=("$p")
    done
    if (( ${#missing[@]} > 0 )); then
        echo "error: --start-step $START_STEP requires these prior outputs, but they are missing:" >&2
        for p in "${missing[@]}"; do
            echo "         $p" >&2
        done
        echo "       Re-run the earlier step, or pick a lower --start-step." >&2
        exit 3
    fi
}

if [[ "$START_STEP" == "3" ]]; then
    _check_prereqs \
        "$RUN_DIR/spatial_adata/${SAMPLE_ID}_proseg_raw.h5ad" \
        "$RUN_DIR/spatial_adata/${SAMPLE_ID}_xenium_ranger.h5ad" \
        "$RUN_DIR/rctd/${SAMPLE_ID}_test_object.rds"
elif [[ "$START_STEP" == "4" ]]; then
    # Explicit --test-object/--reference-rds bypass ALL step-4 prereq checks:
    # the two rds files are supplied out-of-band, and the spatial_adata h5ads
    # step 4 augments in-place may also be brought in externally (or produced
    # under the fresh run folder created by the mkdir -p at LOG_DIR time).
    # Downstream failure surfaces immediately in step-4's own log if anything
    # is actually missing at runtime.
    if [[ -n "$TEST_OBJECT" && -n "$REFERENCE_RDS" ]]; then
        :
    else
        _check_prereqs \
            "$RUN_DIR/spatial_adata/${SAMPLE_ID}_proseg_raw.h5ad" \
            "$RUN_DIR/spatial_adata/${SAMPLE_ID}_xenium_ranger.h5ad" \
            "$RUN_DIR/rctd/${SAMPLE_ID}_test_object.rds" \
            "$RUN_DIR/rctd/${SAMPLE_ID}_reference.rds"
    fi
fi

# ---------------------------------------------------------------------------
# Build the --export payload common to all submitted jobs. RUN_ID is
# APPENDED only when we can bind it up front (either --run-id/RUN_ID given,
# or after JOB1 is submitted). Extra step-1 inputs are threaded via env
# vars consumed by submit_step1.sbatch.
# ---------------------------------------------------------------------------

# Common exports (SAMPLE, OUTPUT_ROOT, FLEX_H5AD, CELLTYPE_MARKER_JSON
# threaded through — even to steps that don't consume them, so an operator
# inspecting a downstream job knows the origin).
COMMON_EXPORTS="ALL,SAMPLE=$SAMPLE_ID,OUTPUT_ROOT=$OUTPUT_ROOT,FLEX_H5AD=$FLEX_H5AD,CELLTYPE_MARKER_JSON=$CELLTYPE_MARKER_JSON"
if [[ -n "$PROSEG_DIR" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,PROSEG_DIR=$PROSEG_DIR"
fi
if [[ -n "$XENIUM_CELLS" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,XENIUM_CELLS=$XENIUM_CELLS"
fi
if [[ -n "$XENIUM_RANGER_DIR" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,XENIUM_RANGER_DIR=$XENIUM_RANGER_DIR"
fi
# Donor + fallback donor h5ads: numbered env vars, one per path. Empty
# arrays => COUNT=0 (step 3 skips the corresponding --donor-h5ad /
# --fallback-donor-h5ad thread). We ALWAYS emit COUNT (including 0) so
# submit_step3.sbatch can rely on `${DONOR_H5AD_COUNT:-0}` returning a
# canonical value rather than a leaked-in stale one from the shell env.
COMMON_EXPORTS="$COMMON_EXPORTS,DONOR_H5AD_COUNT=${#DONOR_H5ADS[@]}"
for i in "${!DONOR_H5ADS[@]}"; do
    COMMON_EXPORTS="$COMMON_EXPORTS,DONOR_H5AD_$((i+1))=${DONOR_H5ADS[$i]}"
done
COMMON_EXPORTS="$COMMON_EXPORTS,FALLBACK_H5AD_COUNT=${#FALLBACK_H5ADS[@]}"
for i in "${!FALLBACK_H5ADS[@]}"; do
    COMMON_EXPORTS="$COMMON_EXPORTS,FALLBACK_H5AD_$((i+1))=${FALLBACK_H5ADS[$i]}"
done
# Step-3 celltype-column override — thread only when set, so submit_step3.sbatch
# falls through to ref-build's own default when the caller omits the flag.
if [[ -n "$CELLTYPE_COL_FOR_REF_BUILD" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,CELLTYPE_COL_FOR_REF_BUILD=$CELLTYPE_COL_FOR_REF_BUILD"
fi
# Step-4 max-cores override — thread only when set, so submit_step4.sbatch
# falls through to its own MAX_CORES default (12) when the caller omits both
# --max-cores and $MAX_CORES env.
if [[ -n "$MAX_CORES_OVERRIDE" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,MAX_CORES=$MAX_CORES_OVERRIDE"
fi
# Step-4 explicit inputs (--test-object / --reference-rds). Validated
# both-or-neither above, so either both are set or neither is; the sbatch
# script gates on `[[ -n "$STEP4_TEST_OBJECT" ]]` and appends the two flags
# together when they arrive. Unset ⇒ rctd-split's own layout auto-discovery.
if [[ -n "$TEST_OBJECT" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP4_TEST_OBJECT=$TEST_OBJECT"
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP4_REFERENCE_RDS=$REFERENCE_RDS"
fi

# ---------------------------------------------------------------------------
# Per-step named-flag exports (STEPN_<PARAM>=<value>) — thread only when set,
# so each sbatch script's `if [[ -n "$STEPN_..." ]]` gate falls through to the
# step CLI's own default when the caller omits the flag. Values are always
# atomic (int / bool-string / short identifier) — never contain commas — so
# they thread safely through slurm's `--export=ALL,K=V,K=V` payload.
# ---------------------------------------------------------------------------
# Step 1
if [[ -n "$STEP1_X_SOURCE" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP1_X_SOURCE=$STEP1_X_SOURCE"
fi
if [[ -n "$STEP1_QC_MIN_COUNTS_CELL" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP1_QC_MIN_COUNTS_CELL=$STEP1_QC_MIN_COUNTS_CELL"
fi
if [[ -n "$STEP1_GEX_ONLY" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP1_GEX_ONLY=$STEP1_GEX_ONLY"
fi
if [[ "$STEP1_FORCE_RERUN" -eq 1 ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP1_FORCE_RERUN=1"
fi
# Step 3
if [[ -n "$STEP3_DONOR_BORROW_CAP" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP3_DONOR_BORROW_CAP=$STEP3_DONOR_BORROW_CAP"
fi
if [[ -n "$STEP3_CELL_MIN_INSTANCE" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP3_CELL_MIN_INSTANCE=$STEP3_CELL_MIN_INSTANCE"
fi
if [[ -n "$STEP3_MIN_UMI" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP3_MIN_UMI=$STEP3_MIN_UMI"
fi
if [[ -n "$STEP3_RANDOM_SEED" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP3_RANDOM_SEED=$STEP3_RANDOM_SEED"
fi
if [[ -n "$STEP3_CELLTYPE_TARGET_LIST" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP3_CELLTYPE_TARGET_LIST=$STEP3_CELLTYPE_TARGET_LIST"
fi
# Step 4
if [[ -n "$STEP4_UMI_MIN" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP4_UMI_MIN=$STEP4_UMI_MIN"
fi
if [[ -n "$STEP4_COUNTS_MIN" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP4_COUNTS_MIN=$STEP4_COUNTS_MIN"
fi
if [[ -n "$STEP4_CELL_MIN_INSTANCE" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP4_CELL_MIN_INSTANCE=$STEP4_CELL_MIN_INSTANCE"
fi
if [[ -n "$STEP4_DOUBLET_MODE" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP4_DOUBLET_MODE=$STEP4_DOUBLET_MODE"
fi
if [[ -n "$STEP4_POSTPROCESS_MIN_COUNTS" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP4_POSTPROCESS_MIN_COUNTS=$STEP4_POSTPROCESS_MIN_COUNTS"
fi
if [[ "$STEP4_KEEP_INTERMEDIATE" -eq 1 ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP4_KEEP_INTERMEDIATE=1"
fi

# ---------------------------------------------------------------------------
# Per-step config-file passthrough (--stepN-config) and --override
# base64-encoded yaml payload. The yaml is built by pasting the config file
# (if any) on top of the empty dict, then walking each `dotted.key=yaml-val`
# override into a nested dict, then base64-encoding the yaml.safe_dump for
# safe transport through slurm's comma-separated `--export` payload (yaml
# lists, string values with commas, and nested dicts round-trip cleanly).
#
# Both --stepN-config and --override are OPTIONAL — when neither is set for
# a step, STEPN_OVERRIDES_B64 is unset and the sbatch script skips the
# `--config <tmp.yaml>` addendum entirely, so a caller that touches no
# nested config sees the same CLI invocation as before this feature landed.
# ---------------------------------------------------------------------------

_encode_step_overrides() {
    # Print a base64-encoded yaml blob for one step, built from:
    #   arg 1: --stepN-config path (may be empty)
    #   arg 2+: --override entries (dotted-key=yaml-val), zero or more
    # Prints an empty string when there are no overrides at all — the
    # driver uses that to gate whether to emit STEPN_OVERRIDES_B64.
    local config_path="$1"; shift
    if [[ -z "$config_path" && $# -eq 0 ]]; then
        return 0
    fi
    STEP_CONFIG_PATH="$config_path" \
    STEP_OVERRIDE_COUNT="$#" \
    STEP_OVERRIDES_JOINED="$(printf '%s\n' "$@")" \
    python3 - <<'PY'
import base64, os, sys, yaml

cfg = {}
config_path = os.environ.get("STEP_CONFIG_PATH", "")
if config_path:
    with open(config_path) as f:
        cfg = yaml.safe_load(f) or {}

# Newline-separated is safe: yaml-parsed override values are printed via
# `printf '%s\n'`, which never internally emits newlines (values are
# short scalars / flow-style lists on a single line by driver contract).
count = int(os.environ.get("STEP_OVERRIDE_COUNT", "0"))
if count:
    lines = os.environ.get("STEP_OVERRIDES_JOINED", "").splitlines()
    assert len(lines) == count, (
        f"expected {count} override lines, got {len(lines)}: {lines!r}"
    )
    for kv in lines:
        key, _, val = kv.partition("=")
        parts = key.split(".")
        d = cfg
        for p in parts[:-1]:
            d = d.setdefault(p, {})
        d[parts[-1]] = yaml.safe_load(val)

blob = yaml.safe_dump(cfg, default_flow_style=False).encode()
sys.stdout.write(base64.b64encode(blob).decode())
PY
}

STEP1_OVERRIDES_B64=$(_encode_step_overrides "$STEP1_CONFIG" "${STEP1_OVERRIDES[@]}")
if [[ -n "$STEP1_OVERRIDES_B64" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP1_OVERRIDES_B64=$STEP1_OVERRIDES_B64"
fi
STEP3_OVERRIDES_B64=$(_encode_step_overrides "$STEP3_CONFIG" "${STEP3_OVERRIDES[@]}")
if [[ -n "$STEP3_OVERRIDES_B64" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP3_OVERRIDES_B64=$STEP3_OVERRIDES_B64"
fi
STEP4_OVERRIDES_B64=$(_encode_step_overrides "$STEP4_CONFIG" "${STEP4_OVERRIDES[@]}")
if [[ -n "$STEP4_OVERRIDES_B64" ]]; then
    COMMON_EXPORTS="$COMMON_EXPORTS,STEP4_OVERRIDES_B64=$STEP4_OVERRIDES_B64"
fi

_sbatch() {
    # Wrapper for sbatch that honours --dry-run: prints the command
    # (whitespace-separated) and emits a fake --parsable id when in
    # dry-run mode. Real sbatch call otherwise.
    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "[dry-run] sbatch $*" >&2
        # Fake, monotonically-increasing job id so downstream deps stay
        # coherent when printed.
        _DRY_JOBID=$((${_DRY_JOBID:-999999} + 1))
        echo "$_DRY_JOBID"
    else
        sbatch "$@"
    fi
}

# ---------------------------------------------------------------------------
# Slurm log routing (settylab/TracyY123-nexus#26 comments 5259180881 +
# 5274187257). Route step-N stdout/stderr into the run-scoped
# <output_root>/<sample_id>/<sample_id>_<run_id>/logs/slurm-<jobid>-<stage>.log
# where <stage> is the pipeline's package name — xenium-preprocess (step 1),
# ref-build (step 3), rctd-split (step 4) — not the internal "stepN" label.
# The step-numbering variables in this script (JOB1/JOB3/JOB4, --start-step,
# step 1/3/4 summary lines, submit_stepN.sbatch filenames, config.yaml
# stepN: keys) are unchanged — only the on-disk slurm log suffix.
# Path template depends on when RUN_ID is bound:
#   * Override case (--run-id / $RUN_ID): full path known up front —
#     mkdir the logs dir and pass an explicit --output= for every job.
#   * Auto case (RUN_ID = JOB1's SLURM_JOB_ID): JOB1 uses sbatch's %j
#     substitution in BOTH the run folder and filename slots so the
#     runtime path resolves correctly on the compute node. Immediately
#     after --parsable returns JOB1's id, mkdir the full logs/ folder
#     (races the job's start; slurm queue latency gives us plenty of
#     headroom). JOB3/JOB4 then get the explicit path.
# ---------------------------------------------------------------------------

if [[ -n "$RUN_ID_OVERRIDE" ]]; then
    LOG_DIR="$OUTPUT_ROOT/$SAMPLE_ID/${SAMPLE_ID}_${RUN_ID_OVERRIDE}/logs"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        mkdir -p "$LOG_DIR"
    fi
    JOB1_OUTPUT="$LOG_DIR/slurm-%j-xenium-preprocess.log"
else
    # LOG_DIR is only knowable AFTER JOB1's --parsable id comes back.
    LOG_DIR=""
    JOB1_OUTPUT="$OUTPUT_ROOT/$SAMPLE_ID/${SAMPLE_ID}_%j/logs/slurm-%j-xenium-preprocess.log"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        mkdir -p "$OUTPUT_ROOT/$SAMPLE_ID"
    fi
fi

# ---------------------------------------------------------------------------
# Submit step 1 (unless --start-step > 1). If the caller pinned a RUN_ID,
# thread it. Otherwise omit it: submit_step1.sbatch's own internal default
# (RUN_ID=${RUN_ID:-$SLURM_JOB_ID}) picks up JOB1's SLURM_JOB_ID.
# ---------------------------------------------------------------------------

JOB1=""
if [[ "$START_STEP" -le 1 ]]; then
    if [[ -n "$RUN_ID_OVERRIDE" ]]; then
        JOB1_EXPORTS="$COMMON_EXPORTS,RUN_ID=$RUN_ID_OVERRIDE"
    else
        JOB1_EXPORTS="$COMMON_EXPORTS"
    fi
    JOB1=$(_sbatch --parsable \
        --output="$JOB1_OUTPUT" \
        --export="$JOB1_EXPORTS" \
        "$SCRIPT_DIR/submit_step1.sbatch")
fi

# Bind RUN_ID for downstream jobs.
if [[ -n "$RUN_ID_OVERRIDE" ]]; then
    RUN_ID="$RUN_ID_OVERRIDE"
else
    # No override implies we submitted step 1 (--start-step > 1 blocks the
    # no-override path above), so JOB1 is set and its id names the run.
    RUN_ID="$JOB1"
    LOG_DIR="$OUTPUT_ROOT/$SAMPLE_ID/${SAMPLE_ID}_${RUN_ID}/logs"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        mkdir -p "$LOG_DIR"
    fi
fi
DOWNSTREAM_EXPORTS="$COMMON_EXPORTS,RUN_ID=$RUN_ID"

# ---------------------------------------------------------------------------
# Step 3 — afterok:JOB1, unless --start-step >= 4 (skip entirely) or
# --start-step 3 (submit with no dependency, since step 1 was skipped).
# ---------------------------------------------------------------------------

JOB3=""
if [[ "$START_STEP" -le 3 ]]; then
    _dep_args=()
    if [[ -n "$JOB1" ]]; then
        _dep_args=(--dependency="afterok:$JOB1")
    fi
    JOB3=$(_sbatch --parsable \
        "${_dep_args[@]}" \
        --output="$LOG_DIR/slurm-%j-ref-build.log" \
        --export="$DOWNSTREAM_EXPORTS" \
        "$SCRIPT_DIR/submit_step3.sbatch")
fi

# ---------------------------------------------------------------------------
# Step 4 — afterok:JOB3, unless --start-step 4 (submit with no dependency,
# since step 3 was skipped).
# ---------------------------------------------------------------------------

_dep_args=()
if [[ -n "$JOB3" ]]; then
    _dep_args=(--dependency="afterok:$JOB3")
fi
JOB4=$(_sbatch --parsable \
    "${_dep_args[@]}" \
    --output="$LOG_DIR/slurm-%j-rctd-split.log" \
    --export="$DOWNSTREAM_EXPORTS" \
    "$SCRIPT_DIR/submit_step4.sbatch")

# ---------------------------------------------------------------------------
# Report — echoed to the caller AND mirrored to <run>/logs/workflow-submit.log
# so an operator inspecting the run folder later has an authoritative record
# of which slurm jobs made up the run (jobids, dependency chain, submit time).
# Lines for skipped steps say "skipped" instead of a jobid so the summary
# still records which steps this invocation covered.
# ---------------------------------------------------------------------------

_fmt_step() {
    local jobid="$1" dep="$2"
    if [[ -z "$jobid" ]]; then
        echo "skipped (--start-step $START_STEP)"
    elif [[ -n "$dep" ]]; then
        echo "$jobid   ($dep)"
    else
        echo "$jobid   (no dependency)"
    fi
}

_step1_line=$(_fmt_step "$JOB1" "")
_step3_line=$(_fmt_step "$JOB3" "${JOB1:+afterok:$JOB1}")
_step4_line=$(_fmt_step "$JOB4" "${JOB3:+afterok:$JOB3}")

_summary=$(cat <<EOF
Workflow chain submitted:
  sample_id  = $SAMPLE_ID
  run_id     = $RUN_ID
  output_dir = $OUTPUT_ROOT/$SAMPLE_ID/${SAMPLE_ID}_${RUN_ID}/
  start_step = $START_STEP
  step 1     = $_step1_line
  step 3     = $_step3_line
  step 4     = $_step4_line
  submitted  = $(date -Iseconds 2>/dev/null || date)
EOF
)
echo "$_summary"
if [[ "$DRY_RUN" -eq 0 ]]; then
    printf '%s\n' "$_summary" >> "$LOG_DIR/workflow-submit.log"
fi
