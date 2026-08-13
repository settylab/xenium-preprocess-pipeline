#!/bin/bash -l
# ---------------------------------------------------------------------
# Slurm submission wrapper for the `rctd-split` pipeline (step 4 of the
# Xenium spatial-data preprocessing pipeline).
#
# `bash -l` makes this a LOGIN shell so ~/.bash_profile (and indirectly
# ~/.bashrc on most setups) gets sourced — that's what initialises
# micromamba in interactive sessions but is otherwise skipped in
# non-interactive Slurm batch jobs.
#
# Usage:
#   sbatch scripts/submit.slurm.sh SAMPLE_ID TEST_OBJECT_RDS REFERENCE_RDS [--key value ...]
#
# Example:
#   sbatch scripts/submit.slurm.sh SAMPLE1 \
#       /data/SAMPLE1/xenium_preprocess/rctd_prep/test_object.rds \
#       /data/SAMPLE1/ref_build/rctd_reference/SAMPLE1_scRNA_ref.rds \
#       --output-root /data/rctd_split_runs
#
# Logs land at: <output_root>/<sample_id>/logs/<job_name>_<jobid>.log
#
# No #SBATCH --output / --error directives — Slurm parses those at
# submit time using your current cwd as the base, and if the cwd isn't
# writable the log file can't be created and your job fails silently.
# Instead, we redirect everything to a path under YOUR output root once
# we know it (which the pipeline owns and `mkdir -p`s on demand).
# ---------------------------------------------------------------------
#SBATCH --job-name=rctd-split
#SBATCH --partition=YOUR_PARTITION
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=6-00:00:00

set -euo pipefail

# ----- Positional args (parsed FIRST so we can route the log) ---------
SAMPLE_ID="${1:?usage: sbatch scripts/submit.slurm.sh SAMPLE_ID TEST_OBJECT_RDS REFERENCE_RDS [extra --flags]}"
TEST_OBJECT="${2:?missing TEST_OBJECT_RDS}"
REFERENCE_RDS="${3:?missing REFERENCE_RDS}"
shift $(( $# < 3 ? $# : 3 ))
EXTRA_ARGS=("$@")

# ----- Resolve OUTPUT_ROOT (mirrors rctd-split precedence) ------------
# Precedence: --output-root in EXTRA_ARGS > $OUTPUT_ROOT env > default.
DEFAULT_OUTPUT_ROOT="/data/rctd_split_runs"
OUTPUT_ROOT="${OUTPUT_ROOT:-$DEFAULT_OUTPUT_ROOT}"
for ((i=0; i<${#EXTRA_ARGS[@]}; i++)); do
    if [[ "${EXTRA_ARGS[i]}" == "--output-root" && $((i+1)) -lt ${#EXTRA_ARGS[@]} ]]; then
        OUTPUT_ROOT="${EXTRA_ARGS[i+1]}"
        break
    fi
done

# ----- Create the user-facing log location + REDIRECT (early) --------
SAMPLE_OUT="$OUTPUT_ROOT/$SAMPLE_ID"
LOG_DIR="$SAMPLE_OUT/logs"
mkdir -p "$LOG_DIR"

JOB_TAG="${SLURM_JOB_NAME:-rctd-split}_${SLURM_JOB_ID:-local-$(date +%Y%m%d-%H%M%S)}"
LOG_FILE="$LOG_DIR/${JOB_TAG}.log"

# ----- Validate paths BEFORE the tee redirect ------------------------
fail=0
if [[ ! -e "$TEST_OBJECT" ]]; then
    echo "[submit] ERROR: TEST_OBJECT does not exist: $TEST_OBJECT" >&2
    fail=1
fi
if [[ ! -e "$REFERENCE_RDS" ]]; then
    echo "[submit] ERROR: REFERENCE_RDS does not exist: $REFERENCE_RDS" >&2
    fail=1
fi
if [[ $fail -ne 0 ]]; then
    echo "[submit] aborting before redirect; fix the path(s) above and resubmit." >&2
    exit 2
fi

# Redirect EVERYTHING from this point on to $LOG_FILE.
#
# TWO PATHS, chosen by execution context (ported verbatim from step 1's
# xenium-preprocess submit.slurm.sh — see comment there for the deadlock
# root cause that motivates the split):
#
#   sbatch batch script  (SLURM_JOB_ID set AND SLURM_STEP_ID unset)
#       -> plain file redirect: `exec >> "$LOG_FILE" 2>&1`.
#       No tee, no process substitution. python -u keeps output
#       unbuffered, so we don't need tee's line-buffered flush.
#
#   interactive `bash submit.slurm.sh` or `srun bash submit.slurm.sh`
#       -> keep the tee-in-process-substitution for live console echo.
#       stdbuf -oL forces per-line flush.
#
# WHY split? Under sbatch (no controlling terminal + shared-filesystem
# latency), a chatty subprocess can burst output faster than tee can
# flush; tee's stdin pipe buffer (64 KB kernel default) fills, tee's
# write blocks on FS latency, python's next write blocks on the full
# pipe, and with `set -euo pipefail` there is no SIGPIPE escape hatch --
# the whole job deadlocks silently holding its full allocation.
if [[ -n "${SLURM_JOB_ID:-}" && -z "${SLURM_STEP_ID:-}" ]]; then
    exec >> "$LOG_FILE" 2>&1
else
    exec > >(stdbuf -oL -eL tee -a "$LOG_FILE") 2> >(stdbuf -oL -eL tee -a "$LOG_FILE" >&2)
fi

echo "[submit] log file: $LOG_FILE"

# ----- Conda / micromamba / mamba env activation ---------------------
# Slurm batch jobs run non-interactive shells. Interactive setups put
# micromamba's init in ~/.bashrc (which defines a SHELL FUNCTION
# `micromamba`, not a binary) — non-interactive shells don't source
# .bashrc, so the function is missing and `micromamba activate` fails.
ENV_NAME="${ENV_NAME:-xenium}"
activated=0

if [[ -f "$HOME/.bashrc" ]]; then
    echo "[submit] sourcing $HOME/.bashrc"
    set +e; set +u
    # shellcheck disable=SC1091
    source "$HOME/.bashrc"
    set -e; set -u
fi

: "${MAMBA_ROOT_PREFIX:=$HOME/micromamba}"
export MAMBA_ROOT_PREFIX
if ! command -v micromamba >/dev/null 2>&1; then
    for bindir in \
        "$MAMBA_ROOT_PREFIX/bin" \
        "$HOME/.local/bin" \
        "$HOME/micromamba/bin" \
        "/app/software/micromamba/bin"; do
        if [[ -x "$bindir/micromamba" ]]; then
            echo "[submit] adding micromamba bin dir to PATH: $bindir"
            export PATH="$bindir:$PATH"
            break
        fi
    done
fi

if [[ "$(type -t micromamba 2>/dev/null)" != "function" ]]; then
    for hook in \
        "$MAMBA_ROOT_PREFIX/etc/profile.d/micromamba.sh" \
        "$HOME/micromamba/etc/profile.d/micromamba.sh" \
        "$HOME/.local/share/mamba/etc/profile.d/micromamba.sh" \
        "/app/software/micromamba/etc/profile.d/micromamba.sh"; do
        if [[ -f "$hook" ]]; then
            echo "[submit] sourcing micromamba hook: $hook"
            # shellcheck disable=SC1090
            source "$hook"
            break
        fi
    done
fi

if command -v micromamba >/dev/null 2>&1; then
    if micromamba activate "$ENV_NAME" 2>/dev/null; then
        activated=1
        echo "[submit] activated env via micromamba: $ENV_NAME"
    else
        echo "[submit] WARN: micromamba activate '$ENV_NAME' failed; trying mamba/conda"
    fi
fi

if [[ $activated -eq 0 ]] && command -v mamba >/dev/null 2>&1; then
    eval "$(mamba shell hook -s bash 2>/dev/null)" || true
    if mamba activate "$ENV_NAME" 2>/dev/null; then
        activated=1
        echo "[submit] activated env via mamba: $ENV_NAME"
    fi
fi

if [[ $activated -eq 0 ]] && command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" 2>/dev/null || true
    if conda activate "$ENV_NAME" 2>/dev/null; then
        activated=1
        echo "[submit] activated env via conda: $ENV_NAME"
    fi
fi

if [[ $activated -eq 0 ]]; then
    echo "[submit] ERROR: could not activate env '$ENV_NAME'." >&2
    echo "[submit]   Override env name: ENV_NAME=your_env sbatch scripts/submit.slurm.sh ..." >&2
    exit 1
fi

echo "[submit] which python:      $(command -v python)"
echo "[submit] which rctd-split:  $(command -v rctd-split 2>/dev/null || echo 'NOT FOUND — did you `pip install -e .` inside the env?')"
echo "[submit] CONDA_PREFIX:      ${CONDA_PREFIX:-unset}"

# ----- LMOD R for the RCTD + SPLIT + export_mtx stages ---------------
# Three of the four stages shell out to Rscript. The reference
# workflow uses the Lmod R module `fhR/4.4.1-foss-2023b` which
# carries Seurat + Matrix + readr preinstalled. spacexr and SPLIT need
# a one-time user-local install at ~/.claude/r_libs/4.4.1 — see
# docs/install.md.
#
# Run `module load` on its OWN LINE — NEVER pipe it. `module` is a
# shell function; piping forks a subshell and the eval that sets env
# vars is silently discarded (see nexus common-gotchas).
R_MODULE="${R_MODULE:-fhR/4.4.1-foss-2023b}"
if command -v ml >/dev/null 2>&1; then
    echo "[submit] loading R module: $R_MODULE"
    ml "$R_MODULE" || echo "[submit] WARN: 'ml $R_MODULE' failed; rctd_run/split_purify/export_mtx may fail unless Rscript is otherwise on PATH."
elif command -v module >/dev/null 2>&1; then
    module load "$R_MODULE" || echo "[submit] WARN: 'module load $R_MODULE' failed; R stages may fail."
else
    echo "[submit] WARN: no 'ml' / 'module' on PATH. R stages depend on Rscript; make sure it's reachable."
fi
echo "[submit] which Rscript: $(command -v Rscript 2>/dev/null || echo 'NOT FOUND — the R stages will fail')"

# ----- Slurm sanity ---------------------------------------------------
echo "===== SLURM INFO ====="
echo "JOBID:         ${SLURM_JOB_ID:-NA}"
echo "NODELIST:      ${SLURM_NODELIST:-NA}"
echo "HOST:          $(hostname)"
echo "PWD:           $(pwd)"
echo "SUBMIT_DIR:    ${SLURM_SUBMIT_DIR:-NA}"
echo "DATE:          $(date)"
echo "PYTHON:        $(which python)"
echo "OUTPUT_ROOT:   $OUTPUT_ROOT"
echo "SAMPLE_OUT:    $SAMPLE_OUT"
echo "LOG_FILE:      $LOG_FILE"
echo "TEST_OBJECT:   $TEST_OBJECT"
echo "REFERENCE_RDS: $REFERENCE_RDS"
python --version || true
echo "======================"

# ----- Run ------------------------------------------------------------
ARGS=(
    run
    --sample-id "$SAMPLE_ID"
    --test-object "$TEST_OBJECT"
    --reference-rds "$REFERENCE_RDS"
)
# Only auto-add --output-root when the user didn't already pass one
# (prevents duplicate --output-root on the CLI).
if ! printf '%s\n' "${EXTRA_ARGS[@]}" | grep -qxF -- "--output-root"; then
    ARGS+=(--output-root "$OUTPUT_ROOT")
fi
ARGS+=( "${EXTRA_ARGS[@]}" )

# Reset PYTHONPATH to ONLY the pipeline source. `ml fhR/…` prepends
# LMOD SciPy-bundle/numpy paths compiled against Python 3.11, which
# corrupt the env's Python 3.13 numpy import (C-extension mismatch).
export PYTHONPATH="<path-to-xenium-preprocess-pipeline>/packages/rctd-split/src"

echo "[submit] running:"
echo "    python -u -m rctd_split.cli ${ARGS[*]}"

# PYTHONUNBUFFERED=1 + python -u together force unbuffered stdout/stderr.
# Without this, Python's prints sit in OS pipe buffers and are LOST when
# Slurm SIGKILLs at time-limit.
export PYTHONUNBUFFERED=1
"${CONDA_PREFIX:-~/micromamba/envs/xenium}/bin/python" -u -m rctd_split.cli "${ARGS[@]}"

echo "[submit] DONE at $(date)"
