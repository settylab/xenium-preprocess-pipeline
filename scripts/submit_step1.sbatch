#!/usr/bin/env bash
#SBATCH --job-name=xp_step1
#SBATCH --partition=campus-new
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
# NOTE: --output is set by submit_workflow.sh at submission time so the log
#       lands under <output_root>/<sample>/<sample>_<run_id>/logs/, not the
#       driver's CWD. A direct `sbatch submit_step1.sbatch` invocation will
#       fall through to slurm's default (slurm-<jobid>.log in the caller's
#       CWD) — invoke via submit_workflow.sh instead.
# submit_step1.sbatch — spatial preprocessing (xenium-preprocess).
#
# Reads env vars set by submit_workflow.sh --export:
#     SAMPLE               (required)  Sample id, e.g. MH10.
#     OUTPUT_ROOT          (required)  Root output directory.
#     RUN_ID               (optional)  Shared run id; defaults to this
#                                      job's SLURM_JOB_ID when unset.
#     PROSEG_DIR           (optional)  Overrides xenium-preprocess proseg-dir.
#     XENIUM_CELLS         (optional)  Overrides xenium-preprocess xenium-cells.
#     XENIUM_RANGER_DIR    (optional)  Overrides xenium-preprocess xenium-ranger-dir.
#     FLEX_H5AD            (unused here) Passed through by the driver so
#                                      operators inspecting env can trace
#                                      the origin; ref-build consumes it.
#     STEP1_X_SOURCE       (optional)  Passed as --x-source.
#     STEP1_QC_MIN_COUNTS_CELL (opt.)  Passed as --qc-min-counts-cell.
#     STEP1_GEX_ONLY       (optional)  Passed as --gex-only.
#     STEP1_FORCE_RERUN    (optional)  If '1', passes --force-rerun.
#     STEP1_STAGE_COUNT    (optional)  Number of --stages entries
#                                      (default 0 ⇒ xenium-preprocess's
#                                      DEFAULT_STAGES).
#     STEP1_STAGE_N        (optional)  Nth stage name (N in 1..COUNT).
#                                      Each is passed as one positional
#                                      value to `xenium-preprocess run --stages …`.
#     STEP1_OVERRIDES_B64  (optional)  base64-encoded YAML built from the
#                                      driver's --xenium-preprocess-config +
#                                      --override xenium_preprocess.<key>=<val>. When
#                                      set, decoded to a tmpfile and
#                                      passed to `xenium-preprocess run
#                                      --config <tmp>` before the named
#                                      flags (so named flags win).
#
# Writes:
#     <OUTPUT_ROOT>/<SAMPLE>/<SAMPLE>_<RUN_ID>/spatial_adata/*.h5ad
#     <OUTPUT_ROOT>/<SAMPLE>/<SAMPLE>_<RUN_ID>/rctd/<SAMPLE>_test_object.rds
#     <OUTPUT_ROOT>/<SAMPLE>/<SAMPLE>_<RUN_ID>/config.yaml (xenium_preprocess: key)
#     <OUTPUT_ROOT>/<SAMPLE>/<SAMPLE>_<RUN_ID>/logs/xenium-preprocess.log

set -euo pipefail

# --- Python + R environments -------------------------------------------------
# Tracy's env is named `xenium` (issue #26 comment 5258774668). Override
# XENIUM_PREPROCESS_ENV in your calling shell (or here) if you named it
# differently. Tries micromamba first, falls back to conda, so the same block
# works for either. `ml fhR/4.4.1-foss-2023b` provides Seurat + spacexr for
# the stage-4 R shell-out (see scratch/xenium-preprocess/README.md).
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
    echo "       See scratch/xenium-preprocess/README.md § Installation." >&2
    exit 4
fi
ml fhR/4.4.1-foss-2023b
# fhR pulls in Python/3.11 + SciPy-bundle/2023.11 as deps and sets
# PYTHONPATH to their python3.11/site-packages. The micromamba env's
# python (3.13) then imports numpy from that path and crashes on the
# C-extension mismatch (see settylab/TracyY123-nexus#26 comment
# 5258986672). R has no need for PYTHONPATH — drop it.
unset PYTHONPATH

# Diagnostic breadcrumb — makes any downstream failure log self-explaining.
# `command -v python` reports Lmod's Python (which the shebang overrides),
# so also record the shebang-python that xenium-preprocess will actually use.
echo "[step1] env=$_env_name  python=$(command -v python || echo '<none>')"
echo "[step1] xenium-preprocess -> $(command -v xenium-preprocess || echo '<none>')"

# Preflight: the step's CLI must actually be installed in the activated env.
# Without this check a missing install surfaces only as an opaque
# `xenium-preprocess: command not found` (see settylab/TracyY123-nexus#26
# comment 5258859987) with no hint at the fix.
_cli="xenium-preprocess"
_pkg_src="/fh/fast/setty_m/user/ryang/nexus/nexus-code/scratch/xenium-preprocess"
if ! command -v "$_cli" >/dev/null 2>&1; then
    cat >&2 <<EOF
error: '$_cli' is not installed in env '$_env_name'.
       Install it once, from a login node:
           micromamba activate $_env_name
           pip install -e $_pkg_src
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

# --config passthrough from driver's --xenium-preprocess-config + --override
# xenium_preprocess.<key>=<val>. When neither was set the driver leaves
# STEP1_OVERRIDES_B64 unset, and this block is a no-op.
if [[ -n "${STEP1_OVERRIDES_B64:-}" ]]; then
    STEP1_TMP_CFG=$(mktemp -t step1_cfg.XXXXXX.yaml)
    trap 'rm -f "$STEP1_TMP_CFG"' EXIT
    printf '%s' "$STEP1_OVERRIDES_B64" | base64 -d > "$STEP1_TMP_CFG"
    extra+=(--config "$STEP1_TMP_CFG")
fi

# Named per-step params — each maps 1:1 to an existing xenium-preprocess CLI
# flag. Materialized AFTER --config so xenium-preprocess's own precedence
# (CLI flag > user YAML > default.yaml) lets these override the yaml
# transparently.
if [[ -n "${STEP1_X_SOURCE:-}" ]]; then
    extra+=(--x-source "$STEP1_X_SOURCE")
fi
if [[ -n "${STEP1_QC_MIN_COUNTS_CELL:-}" ]]; then
    extra+=(--qc-min-counts-cell "$STEP1_QC_MIN_COUNTS_CELL")
fi
if [[ -n "${STEP1_GEX_ONLY:-}" ]]; then
    extra+=(--gex-only "$STEP1_GEX_ONLY")
fi
if [[ "${STEP1_FORCE_RERUN:-0}" == "1" ]]; then
    extra+=(--force-rerun)
fi

# Expand the numbered STEP1_STAGE_* env vars set by submit_workflow.sh
# --xenium-preprocess-stages into a single `--stages s1 s2 ...` block.
# COUNT=0 (or unset) yields zero extra args — xenium-preprocess falls
# through to its DEFAULT_STAGES.
if (( ${STEP1_STAGE_COUNT:-0} > 0 )); then
    extra+=(--stages)
    for i in $(seq 1 "${STEP1_STAGE_COUNT}"); do
        var="STEP1_STAGE_$i"
        extra+=("${!var}")
    done
fi

xenium-preprocess run \
    --sample-id   "$SAMPLE" \
    --output-root "$OUTPUT_ROOT" \
    --run-id      "$RUN_ID" \
    "${extra[@]}"
