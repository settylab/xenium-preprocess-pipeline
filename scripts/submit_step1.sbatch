#!/usr/bin/env bash
#SBATCH --job-name=xp_step1
#SBATCH --partition=YOUR_PARTITION
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
# NOTE: --output is set by submit_workflow.sh at submission time so the log
#       lands under <output_root>/<sample>/<sample>_<run_id>/logs/, not the
#       driver's CWD. A direct `sbatch submit_step1.sbatch` invocation will
#       fall through to slurm's default (slurm-<jobid>.out in the caller's
#       CWD) — invoke via submit_workflow.sh instead.
# submit_step1.sbatch — spatial preprocessing (xenium-preprocess).
#
# Reads env vars set by submit_workflow.sh --export:
#     SAMPLE               (required)  Sample id, e.g. SAMPLE1.
#     OUTPUT_ROOT          (required)  Root output directory.
#     RUN_ID               (optional)  Shared run id; defaults to this
#                                      job's SLURM_JOB_ID when unset.
#     PROSEG_DIR           (optional)  Overrides step-1 proseg-dir.
#     XENIUM_CELLS         (optional)  Overrides step-1 xenium-cells.
#     XENIUM_RANGER_DIR    (optional)  Overrides step-1 xenium-ranger-dir.
#     FLEX_H5AD            (unused here) Passed through by the driver so
#                                      operators inspecting env can trace
#                                      the origin; step 3 consumes it.
#
# Writes:
#     <OUTPUT_ROOT>/<SAMPLE>/<SAMPLE>_<RUN_ID>/spatial_adata/*.h5ad
#     <OUTPUT_ROOT>/<SAMPLE>/<SAMPLE>_<RUN_ID>/rctd/<SAMPLE>_test_object.rds
#     <OUTPUT_ROOT>/<SAMPLE>/<SAMPLE>_<RUN_ID>/resolved_config.yaml (step1: key)
#     <OUTPUT_ROOT>/<SAMPLE>/<SAMPLE>_<RUN_ID>/logs/step1.log

set -euo pipefail

# --- Python + R environments -------------------------------------------------
# The default env name is `xenium` (internal issue review). Override
# XENIUM_PREPROCESS_ENV in your calling shell (or here) if you named it
# differently. Tries micromamba first, falls back to conda, so the same block
# works for either. `ml fhR/4.4.1-foss-2023b` provides Seurat + spacexr for
# the stage-4 R shell-out (see an internal source pathREADME.md).
export PATH="$HOME/.local/bin:$PATH"        # so `micromamba`/`conda` resolves
_env_name="${XENIUM_PREPROCESS_ENV:-xenium}"
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
# `command -v python` reports Lmod's Python (which the shebang overrides),
# so also record the shebang-python that xenium-preprocess will actually use.
echo "[step1] env=$_env_name  python=$(command -v python || echo '<none>')"
echo "[step1] xenium-preprocess -> $(command -v xenium-preprocess || echo '<none>')"

# Preflight: the step's CLI must actually be installed in the activated env.
# Without this check a missing install surfaces only as an opaque
# `xenium-preprocess: command not found` (see (internal issue review)
# comment (internal)) with no hint at the fix.
_cli="xenium-preprocess"
_pkg_src="${XENIUM_PREPROCESS_PIPELINE_ROOT:-<path-to-xenium-preprocess-pipeline>}/packages/xenium-preprocess"
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

# Assemble optional flags without word-splitting surprises.
extra=()
if [[ -n "${PROSEG_DIR:-}" ]]; then
    extra+=(--proseg-dir "$PROSEG_DIR")
fi
if [[ -n "${XENIUM_CELLS:-}" ]]; then
    extra+=(--xenium-cells "$XENIUM_CELLS")
fi
if [[ -n "${XENIUM_RANGER_DIR:-}" ]]; then
    extra+=(--xenium-ranger-dir "$XENIUM_RANGER_DIR")
fi

xenium-preprocess run \
    --sample-id   "$SAMPLE" \
    --output-root "$OUTPUT_ROOT" \
    --run-id      "$RUN_ID" \
    "${extra[@]}"
