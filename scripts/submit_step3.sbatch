#!/usr/bin/env bash
#SBATCH --job-name=rb_step3
#SBATCH --partition=YOUR_PARTITION
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=6:00:00
# NOTE: --output is set by submit_workflow.sh at submission time so the log
#       lands under <output_root>/<sample>/<sample>_<run_id>/logs/, not the
#       driver's CWD. A direct `sbatch submit_step3.sbatch` invocation will
#       fall through to slurm's default (slurm-<jobid>.out in the caller's
#       CWD) — invoke via submit_workflow.sh instead.
# submit_step3.sbatch — RCTD reference build (ref-build).
#
# Reads env vars set by submit_workflow.sh --export:
#     SAMPLE               (required)  Sample id.
#     OUTPUT_ROOT          (required)  Root output directory.
#     FLEX_H5AD            (required)  Absolute path to the flex scRNA h5ad;
#                                      recorded verbatim under
#                                      step3.flex_h5ad_path in
#                                      resolved_config.yaml. NO copy, NO
#                                      symlink (reference-in-place per
#                                      comment (internal) §1).
#     CELLTYPE_MARKER_JSON (required)  Marker-gene JSON declaring the
#                                      expected celltype set; threaded to
#                                      ref-build run --celltype-marker-json.
#     RUN_ID               (required)  Shared run id (bound by JOB1 or by
#                                      the caller). Defaults to this job's
#                                      SLURM_JOB_ID only if the driver was
#                                      bypassed and left it unset.
#     DONOR_H5AD_COUNT     (optional)  Number of donor h5ads (default 0).
#     DONOR_H5AD_N         (optional)  Nth donor h5ad path (N in 1..COUNT).
#                                      Each is passed as one --donor-h5ad
#                                      arg to `ref-build run`.
#     FALLBACK_H5AD_COUNT  (optional)  Number of fallback donor h5ads (0).
#     FALLBACK_H5AD_N      (optional)  Nth fallback donor h5ad path;
#                                      passed as --fallback-donor-h5ad.
#     CELLTYPE_COL_FOR_REF_BUILD (optional)  .obs column ref-build reads
#                                      for the per-cell celltype label.
#                                      Passed as --celltype-col to
#                                      `ref-build run`. Unset ⇒ ref-build's
#                                      config default
#                                      (Final_level1_celltype_annotation).
#
# Writes:
#     <OUTPUT_ROOT>/<SAMPLE>/<SAMPLE>_<RUN_ID>/rctd/<SAMPLE>_reference_post_rules.h5ad
#     <OUTPUT_ROOT>/<SAMPLE>/<SAMPLE>_<RUN_ID>/rctd/<SAMPLE>_reference.rds
#     <OUTPUT_ROOT>/<SAMPLE>/<SAMPLE>_<RUN_ID>/resolved_config.yaml (step3: key)
#     <OUTPUT_ROOT>/<SAMPLE>/<SAMPLE>_<RUN_ID>/logs/step3.log

set -euo pipefail

# --- Python + R environments -------------------------------------------------
# The default env name is `xenium` (internal issue review). Override
# REF_BUILD_ENV in your calling shell (or here) if you use a different env
# for ref-build. Tries micromamba first, falls back to conda, so the same
# block works for either. `ml fhR/4.4.1-foss-2023b` provides Seurat + spacexr
# for the rctd_reference_build R stage (see docs/installation.md).
export PATH="$HOME/.local/bin:$PATH"        # so `micromamba`/`conda` resolves
_env_name="${REF_BUILD_ENV:-xenium}"
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
echo "[step3] env=$_env_name  python=$(command -v python || echo '<none>')"
echo "[step3] ref-build -> $(command -v ref-build || echo '<none>')"

# Preflight: the step's CLI must actually be installed in the activated env.
# Without this check a missing install surfaces only as an opaque
# `ref-build: command not found` with no hint at the fix.
_cli="ref-build"
_pkg_src="${XENIUM_PREPROCESS_PIPELINE_ROOT:-<path-to-xenium-preprocess-pipeline>}/packages/ref-build"
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
: "${FLEX_H5AD:?FLEX_H5AD is required (set by submit_workflow.sh --export)}"
: "${CELLTYPE_MARKER_JSON:?CELLTYPE_MARKER_JSON is required (set by submit_workflow.sh --export)}"

RUN_ID="${RUN_ID:-$SLURM_JOB_ID}"

# Expand the numbered DONOR_H5AD_* / FALLBACK_H5AD_* env vars set by
# submit_workflow.sh into repeated --donor-h5ad / --fallback-donor-h5ad
# args. COUNT=0 (or unset) yields zero extra args — safe default.
DONOR_ARGS=()
for i in $(seq 1 "${DONOR_H5AD_COUNT:-0}"); do
    var="DONOR_H5AD_$i"
    DONOR_ARGS+=(--donor-h5ad "${!var}")
done
FALLBACK_ARGS=()
for i in $(seq 1 "${FALLBACK_H5AD_COUNT:-0}"); do
    var="FALLBACK_H5AD_$i"
    FALLBACK_ARGS+=(--fallback-donor-h5ad "${!var}")
done

# Optional celltype-col override — only appended when the driver was
# invoked with --celltype-col-for-ref-build; otherwise ref-build falls
# through to its own config default (Final_level1_celltype_annotation).
CELLTYPE_COL_ARGS=()
if [[ -n "${CELLTYPE_COL_FOR_REF_BUILD:-}" ]]; then
    CELLTYPE_COL_ARGS=(--celltype-col "$CELLTYPE_COL_FOR_REF_BUILD")
fi

ref-build run \
    --sample-id            "$SAMPLE" \
    --output-root          "$OUTPUT_ROOT" \
    --run-id               "$RUN_ID" \
    --flex-h5ad            "$FLEX_H5AD" \
    --celltype-marker-json "$CELLTYPE_MARKER_JSON" \
    "${CELLTYPE_COL_ARGS[@]}" \
    "${DONOR_ARGS[@]}" \
    "${FALLBACK_ARGS[@]}"
