#!/usr/bin/env bash
# submit_workflow.sh — top-level driver for the step-1 → step-3 → step-4
# spatial-genomics workflow (internal issue review).
#
# Submits up to three sbatch jobs and chains them with --dependency=afterok:
#     JOB1 (step 1, xenium-preprocess)
#       └── JOB3 (step 3, ref-build)          afterok:JOB1
#             └── JOB4 (step 4, rctd-split)   afterok:JOB3
#
# --start-step lets an operator resume the chain when earlier steps have
# already produced their outputs in the run folder (internal issue review)
# comment (internal)). --start-step 3 skips JOB1 (step 3 submits with no
# dependency); --start-step 4 skips JOB1 + JOB3.
#
# All submitted jobs share a single RUN_ID, propagated via --export=ALL,RUN_ID=…
# so outputs land in the run-scoped folder:
#     <output_root>/<sample_id>/<sample_id>_<run_id>/
#         ├── spatial_adata/   ← step 1 + step 4 augmentations
#         ├── rctd/            ← step 1 + step 3 + step 4
#         ├── resolved_config.yaml  ← merged across steps
#         └── logs/            ← step{1,3,4}.log
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
#     ./submit_workflow.sh --sample-id SAMPLE1 --flex-h5ad /path/flex.h5ad \
#                          --celltype-marker-json /path/markers.json
#     # Resume from step 3 (step 1 already ran, keep its outputs):
#     ./submit_workflow.sh --sample-id SAMPLE1 --flex-h5ad /path/flex.h5ad \
#                          --celltype-marker-json /path/markers.json \
#                          --run-id my_experiment_v2 \
#                          --start-step 3 --reuse-run-dir
#     # From-scratch re-run under an existing run-id (destructive):
#     ./submit_workflow.sh --sample-id SAMPLE1 --flex-h5ad /path/flex.h5ad \
#                          --celltype-marker-json /path/markers.json \
#                          --run-id my_experiment_v2 --force
#     RUN_ID=my_id ./submit_workflow.sh --sample-id SAMPLE1 \
#                          --flex-h5ad /path/flex.h5ad \
#                          --celltype-marker-json /path/markers.json
#     ./submit_workflow.sh --sample-id SAMPLE1 --flex-h5ad /path/flex.h5ad \
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

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults + arg parsing
# ---------------------------------------------------------------------------

# Script dir — location of submit_stepN.sbatch stubs.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

DEFAULT_OUTPUT_ROOT=${OUTPUT_ROOT:-}

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
# omitted, but the standard flow needs them named).
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

usage() {
    cat <<'EOF'
Usage: submit_workflow.sh --sample-id <S> --flex-h5ad <path>
                          --celltype-marker-json <path> [OPTIONS]

Required:
  --sample-id <S>              Sample identifier (e.g. SAMPLE1).
  --flex-h5ad <path>           Flex scRNA h5ad (recorded verbatim under
                               step3.flex_h5ad_path in resolved_config.yaml;
                               no copy / no symlink).
  --celltype-marker-json <p>   Marker-gene JSON declaring the expected
                               celltype set (step-3 required input;
                               threaded to ref-build run --celltype-marker-json).

Optional:
  --output-root <dir>          Root output directory. Required unless
                               the OUTPUT_ROOT env var is set.
  --run-id <id>                Explicit run identifier. Precedence:
                               --run-id > $RUN_ID env > JOB1 SLURM_JOB_ID.
                               Required when --start-step > 1.
  --start-step <1|3|4>         Skip earlier steps and start submission at
                               step N. Default 1 (full chain). --start-step 3
                               submits step 3 (no dep) then step 4; assumes
                               step-1 outputs already exist in the run
                               folder. --start-step 4 submits only step 4;
                               assumes step-1 + step-3 outputs exist. Both
                               imply --reuse-run-dir. Fails loud if the
                               required prior outputs are missing.
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
if [[ -z "$OUTPUT_ROOT" ]]; then
    echo "error: --output-root is required (or set OUTPUT_ROOT env)" >&2
    usage >&2
    exit 2
fi

# ---------------------------------------------------------------------------
# --start-step validation. Legal values are 1, 3, 4 — matches the pipeline's
# step numbering (there is no step 2 in this workflow). > 1 needs a
# pre-bound RUN_ID (nothing to resume without one) and implies
# --reuse-run-dir (the whole point is to keep the earlier steps' outputs).
# ---------------------------------------------------------------------------

case "$START_STEP" in
    1|3|4) ;;
    *)
        echo "error: --start-step must be one of 1, 3, 4 (got: $START_STEP)" >&2
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
        echo "error: --start-step $START_STEP but the run folder does not exist: $RUN_DIR" >&2
        echo "       There is nothing to resume from — did you mean --start-step 1?" >&2
        exit 3
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
    _check_prereqs \
        "$RUN_DIR/spatial_adata/${SAMPLE_ID}_proseg_raw.h5ad" \
        "$RUN_DIR/spatial_adata/${SAMPLE_ID}_xenium_ranger.h5ad" \
        "$RUN_DIR/rctd/${SAMPLE_ID}_test_object.rds" \
        "$RUN_DIR/rctd/${SAMPLE_ID}_reference.rds"
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
# Slurm log routing (internal issue review)
# (internal issue review)). Route step-N stdout/stderr into the run-scoped
# <output_root>/<sample_id>/<sample_id>_<run_id>/logs/slurm-<jobid>-<stage>.out
# where <stage> is the pipeline's package name — xenium-preprocess (step 1),
# ref-build (step 3), rctd-split (step 4) — not the internal "stepN" label.
# The step-numbering variables in this script (JOB1/JOB3/JOB4, --start-step,
# step 1/3/4 summary lines, submit_stepN.sbatch filenames, resolved_config.yaml
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
    JOB1_OUTPUT="$LOG_DIR/slurm-%j-xenium-preprocess.out"
else
    # LOG_DIR is only knowable AFTER JOB1's --parsable id comes back.
    LOG_DIR=""
    JOB1_OUTPUT="$OUTPUT_ROOT/$SAMPLE_ID/${SAMPLE_ID}_%j/logs/slurm-%j-xenium-preprocess.out"
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
        --output="$LOG_DIR/slurm-%j-ref-build.out" \
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
    --output="$LOG_DIR/slurm-%j-rctd-split.out" \
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
