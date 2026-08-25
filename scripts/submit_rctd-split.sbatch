#!/usr/bin/env bash
#SBATCH --job-name=rctd-split
#SBATCH --partition=campus-new
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=6-00:00:00
# NOTE: --output is set by submit_workflow.sh at submission time so the log
#       lands under <output_root>/<sample>/<sample>_<run_id>/logs/, not the
#       driver's CWD. A direct `sbatch submit_rctd-split.sbatch` invocation will
#       fall through to slurm's default (slurm-<jobid>.log in the caller's
#       CWD) — invoke via submit_workflow.sh instead.
# submit_rctd-split.sbatch — RCTD run + split + purify + writebacks (rctd-split).
#
# Reads env vars set by submit_workflow.sh --export:
#     SAMPLE               (required)  Sample id.
#     OUTPUT_ROOT          (required)  Root output directory.
#     RUN_ID               (required)  Shared run id.
#     MAX_CORES            (optional)  Passed as --max-cores. Default:
#                                      $SLURM_CPUS_PER_TASK (the actual
#                                      --cpus-per-task allocation above),
#                                      falling back to 12 outside Slurm.
#     STEP4_UMI_MIN        (optional)  Passed as --umi-min.
#     STEP4_COUNTS_MIN     (optional)  Passed as --counts-min.
#     STEP4_CELL_MIN_INSTANCE (optional) Passed as --cell-min-instance.
#     STEP4_DOUBLET_MODE   (optional)  Passed as --doublet-mode.
#     STEP4_POSTPROCESS_MIN_COUNTS (opt.) Passed as --postprocess-min-counts.
#     STEP4_KEEP_INTERMEDIATE (optional) If '1', passes --keep-intermediate.
#     STEP4_FORCE_RERUN    (optional)  If '1', passes --force-rerun.
#     STEP4_STAGE_COUNT    (optional)  Number of --stages entries
#                                      (default 0 ⇒ rctd-split's DEFAULT_STAGES).
#     STEP4_STAGE_N        (optional)  Nth stage name (N in 1..COUNT).
#                                      Each is passed as one positional
#                                      value to `rctd-split run --stages …`.
#     STEP4_TEST_OBJECT    (optional)  Explicit test_object.rds path. When
#                                      set, passed as `--test-object <path>`;
#                                      rctd-split then skips the run-folder
#                                      layout auto-discovery for the test
#                                      object. MUTUAL with STEP4_REFERENCE_RDS
#                                      (both-or-neither, enforced upstream).
#     STEP4_REFERENCE_RDS  (optional)  Explicit reference.rds path. When set,
#                                      passed as `--reference-rds <path>`;
#                                      same semantics as STEP4_TEST_OBJECT.
#     STEP4_RCTD_RESULTS_RDS (optional) Explicit rctd_results.rds path. When
#                                      set, passed as `--rctd-results-rds
#                                      <path>`; rctd-split auto-drops the
#                                      `rctd_run` stage (external RCTD
#                                      result supplied) and feeds the file
#                                      straight into split_purify.
#                                      Standalone (not mutual with the
#                                      test_object/reference pair).
#     STEP4_OVERRIDES_B64  (optional)  base64-encoded YAML from the driver's
#                                      --rctd-split-config + --override
#                                      rctd_split.<key>=<val>. When set, decoded
#                                      to a tmpfile and passed to `rctd-split
#                                      run --config <tmp>` before the named
#                                      flags (so named flags win).
#
# When STEP4_TEST_OBJECT / STEP4_REFERENCE_RDS are unset, rctd-split uses the
# run-scoped layout to locate xenium-preprocess test_object.rds and ref-build reference.rds
# automatically:
#     <OUTPUT_ROOT>/<SAMPLE>/<SAMPLE>_<RUN_ID>/rctd/<SAMPLE>_test_object.rds
#     <OUTPUT_ROOT>/<SAMPLE>/<SAMPLE>_<RUN_ID>/rctd/<SAMPLE>_reference.rds
#
# Writes / augments (all in-place under the run folder):
#     rctd/<SAMPLE>_rctd_results.rds                    (new)
#     spatial_adata/<SAMPLE>_proseg_purified.h5ad       (new)
#     spatial_adata/<SAMPLE>_proseg_raw.h5ad            (obs cols augmented
#                                                       by writeback_to_raw)
#     spatial_adata/<SAMPLE>_xenium_ranger.h5ad         (celltype + distance
#                                                       cols augmented by
#                                                       celltype_writeback)
#     config.yaml (rctd_split: key)
#     logs/slurm-<jobid>-rctd-split.log
#       (slurm-captured stdout+stderr; the stub does not write a separate app log)

set -euo pipefail

# --- Python + R environments -------------------------------------------------
# Resolved from scripts/env.local.conf (scripts/write-env-config.sh),
# NEVER re-derived from `$HOME` here: the two failures this pipeline
# has actually hit (missing spacexr's R library; micromamba's own
# root-prefix resolving to the wrong place) were both cases of a
# `$HOME`-relative path resolving differently INSIDE the job than in
# the submitting shell. A file on shared storage, holding already-
# resolved absolute values, is immune to that; `--export=ALL` is not
# (it forwards env VAR values faithfully, but does nothing to stop
# `$HOME` itself from resolving to a different instance once the job's
# own shell starts). See scripts/lib/env_config.sh for the full
# rationale. Activation is by absolute PREFIX, never by name, for the
# same reason.
# PIPELINE_SCRIPT_DIR is set by submit_workflow.sh (the only place
# BASH_SOURCE resolution is reliable — this job's own BASH_SOURCE
# array can be empty under some sbatch launchers, which would make
# `${BASH_SOURCE[0]}` an unbound array element under `set -u`).
# Fall back to self-resolution for a manual/direct `sbatch` invocation.
SCRIPT_DIR="${PIPELINE_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)}"
# shellcheck source=lib/env_config.sh
source "$SCRIPT_DIR/lib/env_config.sh"

eval "$("$MICROMAMBA_BIN" shell hook --shell bash)"
micromamba activate "$XENIUM_ENV_PREFIX"

ml "$R_MODULE"
# `fhR` provides Seurat + Matrix, but NOT spacexr or SPLIT (verified
# against a live `.libPaths()`). Both are supplied via
# `--r-lib-paths "$R_LIB_DIR"` below (config key r_lib_paths — this
# package already had the mechanism; see rctd_split/stages/*.py).
# fhR pulls in Python/3.11 + SciPy-bundle/2023.11 as deps and sets
# PYTHONPATH to their python3.11/site-packages. The micromamba env's
# python (3.13) then imports numpy from that path and crashes on the
# C-extension mismatch (see settylab/TracyY123-nexus#26 comment
# 5258986672). R has no need for PYTHONPATH — drop it.
unset PYTHONPATH

# Diagnostic breadcrumb — makes any downstream failure log self-explaining.
echo "[rctd-split] env=$XENIUM_ENV_PREFIX  python=$(command -v python || echo '<none>')"
echo "[rctd-split] rctd-split -> $(command -v rctd-split || echo '<none>')"

# Preflight: the step's CLI must actually be installed in the activated env.
# Without this check a missing install surfaces only as an opaque
# `rctd-split: command not found` with no hint at the fix.
_cli="rctd-split"
_pkg_src="/fh/fast/setty_m/user/ryang/nexus/nexus-code/scratch/rctd-split"
if ! command -v "$_cli" >/dev/null 2>&1; then
    cat >&2 <<EOF
error: '$_cli' is not installed in env '$XENIUM_ENV_PREFIX'.
       Install it once, from a login node:
           micromamba activate $XENIUM_ENV_PREFIX
           pip install -e $_pkg_src
       Then re-submit the workflow.
EOF
    exit 5
fi
# -----------------------------------------------------------------------------

: "${SAMPLE:?SAMPLE is required (set by submit_workflow.sh --export)}"
: "${OUTPUT_ROOT:?OUTPUT_ROOT is required (set by submit_workflow.sh --export)}"

RUN_ID="${RUN_ID:-$SLURM_JOB_ID}"

# MAX_CORES: RCTD parallelism. Defaults to the ACTUAL allocation
# ($SLURM_CPUS_PER_TASK, set by --cpus-per-task above) rather than a
# constant, so raising/lowering --cpus-per-task changes RCTD's core
# count automatically instead of silently doing nothing (operator
# request, 2026-08-21: "make sure RCTD uses the number of cores as
# requested for job submission"). Trade-off: this reserves NO headroom
# for the parent R process alongside its `max_cores` workers — if you
# want headroom, either bump --cpus-per-task above what you want RCTD
# to use, or export MAX_CORES=<N> explicitly (still wins over both).
# 12 remains the fallback for a direct `sbatch submit_rctd-split.sbatch`
# invocation outside Slurm's allocation context, where
# $SLURM_CPUS_PER_TASK is unset.
#
# Editing rctd_run.max_cores in <run>/config.yaml does NOT flow back
# into the pipeline — that file is a write-only snapshot the pipeline
# overwrites at every stage; THIS env var (MAX_CORES, from
# $SLURM_CPUS_PER_TASK or an explicit export) is what actually wins.
# Override per-run by exporting MAX_CORES=<N> before invoking
# submit_workflow.sh.
MAX_CORES="${MAX_CORES:-${SLURM_CPUS_PER_TASK:-12}}"

# --config passthrough from driver's --rctd-split-config + --override
# rctd_split.<key>=<val>. When neither was set the driver leaves
# STEP4_OVERRIDES_B64 unset, and CONFIG_ARGS stays empty.
CONFIG_ARGS=()
if [[ -n "${STEP4_OVERRIDES_B64:-}" ]]; then
    STEP4_TMP_CFG=$(mktemp -t rctd-split_cfg.XXXXXX.yaml)
    trap 'rm -f "$STEP4_TMP_CFG"' EXIT
    printf '%s' "$STEP4_OVERRIDES_B64" | base64 -d > "$STEP4_TMP_CFG"
    CONFIG_ARGS=(--config "$STEP4_TMP_CFG")
fi

# Named per-step params — each maps 1:1 to an existing rctd-split CLI flag.
# Materialized AFTER --config so rctd-split's own precedence (CLI flag >
# user YAML > default.yaml) lets these override the yaml transparently.
NAMED_ARGS=()
if [[ -n "${STEP4_UMI_MIN:-}" ]]; then
    NAMED_ARGS+=(--umi-min "$STEP4_UMI_MIN")
fi
if [[ -n "${STEP4_COUNTS_MIN:-}" ]]; then
    NAMED_ARGS+=(--counts-min "$STEP4_COUNTS_MIN")
fi
if [[ -n "${STEP4_CELL_MIN_INSTANCE:-}" ]]; then
    NAMED_ARGS+=(--cell-min-instance "$STEP4_CELL_MIN_INSTANCE")
fi
if [[ -n "${STEP4_DOUBLET_MODE:-}" ]]; then
    NAMED_ARGS+=(--doublet-mode "$STEP4_DOUBLET_MODE")
fi
if [[ -n "${STEP4_POSTPROCESS_MIN_COUNTS:-}" ]]; then
    NAMED_ARGS+=(--postprocess-min-counts "$STEP4_POSTPROCESS_MIN_COUNTS")
fi
if [[ "${STEP4_KEEP_INTERMEDIATE:-0}" == "1" ]]; then
    NAMED_ARGS+=(--keep-intermediate)
fi
if [[ "${STEP4_FORCE_RERUN:-0}" == "1" ]]; then
    NAMED_ARGS+=(--force-rerun)
fi
# Explicit test_object / reference_rds — MUTUAL both-or-neither (enforced
# by submit_workflow.sh upstream, so we trust the invariant here and add
# both flags together when the first is present).
if [[ -n "${STEP4_TEST_OBJECT:-}" ]]; then
    NAMED_ARGS+=(--test-object   "$STEP4_TEST_OBJECT"
                 --reference-rds "$STEP4_REFERENCE_RDS")
fi
# Explicit rctd_results.rds — standalone; rctd-split's CLI auto-drops
# the rctd_run stage on its side when this is set.
if [[ -n "${STEP4_RCTD_RESULTS_RDS:-}" ]]; then
    NAMED_ARGS+=(--rctd-results-rds "$STEP4_RCTD_RESULTS_RDS")
fi

# Expand the numbered STEP4_STAGE_* env vars set by submit_workflow.sh
# --rctd-split-stages into a single `--stages s1 s2 ...` block. COUNT=0
# (or unset) yields zero extra args — rctd-split falls through to its
# DEFAULT_STAGES.
STAGES_ARGS=()
if (( ${STEP4_STAGE_COUNT:-0} > 0 )); then
    STAGES_ARGS=(--stages)
    for i in $(seq 1 "${STEP4_STAGE_COUNT}"); do
        var="STEP4_STAGE_$i"
        STAGES_ARGS+=("${!var}")
    done
fi

rctd-split run \
    --sample-id   "$SAMPLE" \
    --output-root "$OUTPUT_ROOT" \
    --run-id      "$RUN_ID" \
    --max-cores   "$MAX_CORES" \
    --r-lib-paths "$R_LIB_DIR" \
    "${CONFIG_ARGS[@]}" \
    "${NAMED_ARGS[@]}" \
    "${STAGES_ARGS[@]}"
