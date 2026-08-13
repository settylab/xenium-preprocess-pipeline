#!/usr/bin/env bash
#SBATCH --job-name=rs_step4
#SBATCH --partition=YOUR_PARTITION
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=24:00:00
# NOTE: --output is set by submit_workflow.sh at submission time so the log
#       lands under <output_root>/<sample>/<sample>_<run_id>/logs/, not the
#       driver's CWD. A direct `sbatch submit_step4.sbatch` invocation will
#       fall through to slurm's default (slurm-<jobid>.out in the caller's
#       CWD) — invoke via submit_workflow.sh instead.
# submit_step4.sbatch — RCTD run + split + purify + writebacks (rctd-split).
#
# Reads env vars set by submit_workflow.sh --export:
#     SAMPLE               (required)  Sample id.
#     OUTPUT_ROOT          (required)  Root output directory.
#     RUN_ID               (required)  Shared run id.
#
# Uses the run-scoped layout to locate step-1 test_object.rds and step-3
# reference.rds automatically:
#     <OUTPUT_ROOT>/<SAMPLE>/<SAMPLE>_<RUN_ID>/rctd/<SAMPLE>_test_object.rds
#     <OUTPUT_ROOT>/<SAMPLE>/<SAMPLE>_<RUN_ID>/rctd/<SAMPLE>_reference.rds
#
# Writes / augments (all in-place under the run folder):
#     rctd/<SAMPLE>_rctd_results.rds                    (new)
#     spatial_adata/<SAMPLE>_proseg_purified.h5ad       (new)
#     spatial_adata/<SAMPLE>_proseg_raw.h5ad            (obs cols augmented
#                                                       by writeback_to_step1_raw)
#     spatial_adata/<SAMPLE>_xenium_ranger.h5ad         (celltype + distance
#                                                       cols augmented by
#                                                       celltype_writeback)
#     resolved_config.yaml (step4: key)
#     logs/step4.log

set -euo pipefail

# --- Python + R environments -------------------------------------------------
# The default env name is `xenium` (internal issue review). Override
# RCTD_SPLIT_ENV in your calling shell (or here) if you use a different env
# for rctd-split. Tries micromamba first, falls back to conda, so the same
# block works for either. `ml fhR/4.4.1-foss-2023b` provides Seurat + spacexr
# + SPLIT for the rctd_run and split_purify R stages (see
# docs/installation.md).
export PATH="$HOME/.local/bin:$PATH"        # so `micromamba`/`conda` resolves
_env_name="${RCTD_SPLIT_ENV:-xenium}"
if command -v micromamba >/dev/null 2>&1; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate "$_env_name"
elif command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate "$_env_name"
else
    echo "error: neither micromamba nor conda found on PATH." >&2
    echo "       See docs/installation.md." >&2
    exit 4
fi
ml fhR/4.4.1-foss-2023b
# fhR pulls in Python/3.11 + SciPy-bundle/2023.11 as deps and sets
# PYTHONPATH to their python3.11/site-packages. The micromamba env's
# python (3.13) then imports numpy from that path and crashes on the
# C-extension mismatch (see (internal issue review) comment
# (internal issue review)). R has no need for PYTHONPATH — drop it.
unset PYTHONPATH

# Diagnostic breadcrumb — makes any downstream failure log self-explaining.
echo "[step4] env=$_env_name  python=$(command -v python || echo '<none>')"
echo "[step4] rctd-split -> $(command -v rctd-split || echo '<none>')"

# Preflight: the step's CLI must actually be installed in the activated env.
# Without this check a missing install surfaces only as an opaque
# `rctd-split: command not found` with no hint at the fix.
_cli="rctd-split"
_pkg_src="${XENIUM_PREPROCESS_PIPELINE_ROOT:-<path-to-xenium-preprocess-pipeline>}/packages/rctd-split"
if ! command -v "$_cli" >/dev/null 2>&1; then
    cat >&2 <<EOF
error: '$_cli' is not installed in env '$_env_name'.
       Install it once, from a login node:
           micromamba activate $_env_name
           uv pip install -e $_pkg_src
       Then re-submit the workflow.
EOF
    exit 5
fi
# -----------------------------------------------------------------------------

: "${SAMPLE:?SAMPLE is required (set by submit_workflow.sh --export)}"
: "${OUTPUT_ROOT:?OUTPUT_ROOT is required (set by submit_workflow.sh --export)}"

RUN_ID="${RUN_ID:-$SLURM_JOB_ID}"

# MAX_CORES: RCTD parallelism. sbatch alloc is 16 CPUs (--cpus-per-task
# above); default here is 12 (headroom for R workers + parent process).
# Editing rctd_run.max_cores in <run>/resolved_config.yaml does NOT flow
# back into the pipeline — that file is a write-only snapshot the pipeline
# overwrites at every stage. Override per-run by exporting MAX_CORES=<N>
# before invoking submit_workflow.sh.
MAX_CORES="${MAX_CORES:-12}"

rctd-split run \
    --sample-id   "$SAMPLE" \
    --output-root "$OUTPUT_ROOT" \
    --run-id      "$RUN_ID" \
    --max-cores   "$MAX_CORES"
