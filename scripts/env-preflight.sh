#!/usr/bin/env bash
# scripts/env-preflight.sh — submit-time preflight for the resolved
# environment. Run this ONCE, before dispatching any sbatch job (the
# driver, submit_workflow.sh, calls it automatically unless
# --skip-preflight is passed) so a broken/stale environment fails in
# seconds instead of after a queue wait.
#
# Checks:
#   1. scripts/env.local.conf exists and every key resolves to a real
#      path (delegated to scripts/lib/env_config.sh — the same check
#      every sbatch stub does at job start, just run here first).
#   2. The R library it names ACTUALLY resolves spacexr + SPLIT — the
#      expensive check (`ml` + an Rscript subprocess), which is why
#      this is a separate, explicitly-invoked step rather than
#      something every sbatch job repeats.
#
# Usage: scripts/env-preflight.sh [--skip-r-check]
#   --skip-r-check   Only run the cheap config-file checks (useful for
#                     xenium-preprocess-only runs, whose R side never
#                     touches spacexr/SPLIT — see rctd_prep.R).

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

SKIP_R_CHECK=0
if [[ "${1:-}" == "--skip-r-check" ]]; then
    SKIP_R_CHECK=1
fi

# shellcheck source=lib/env_config.sh
source "$SCRIPT_DIR/lib/env_config.sh"
echo "[env-preflight] scripts/env.local.conf OK:"
echo "[env-preflight]   MICROMAMBA_BIN=$MICROMAMBA_BIN"
echo "[env-preflight]   XENIUM_ENV_PREFIX=$XENIUM_ENV_PREFIX"
echo "[env-preflight]   R_LIB_DIR=$R_LIB_DIR"
echo "[env-preflight]   R_MODULE=$R_MODULE"

if [[ "$SKIP_R_CHECK" == "1" ]]; then
    echo "[env-preflight] --skip-r-check: not checking spacexr/SPLIT resolution."
    exit 0
fi

if ! command -v ml >/dev/null 2>&1; then
    echo "error: [env-preflight] 'ml' (Lmod) not found on PATH — cannot load $R_MODULE to check R packages." >&2
    echo "       Pass --skip-r-check if this host has no Lmod (off-cluster dev)." >&2
    exit 5
fi

# `ml` is a shell function, not piped — run on its own line (see
# bash-footgun-guard: piping ml silently discards its env-changing eval).
ml "$R_MODULE"

if ! command -v Rscript >/dev/null 2>&1; then
    echo "error: [env-preflight] Rscript not found on PATH after 'ml $R_MODULE'." >&2
    exit 5
fi

R_CHECK_OUT=$(Rscript --vanilla -e '
  args <- commandArgs(trailingOnly = TRUE)
  lib <- args[[1]]
  .libPaths(c(lib, .libPaths()))
  missing <- Filter(function(pkg) !requireNamespace(pkg, quietly = TRUE),
                     c("spacexr", "SPLIT"))
  if (length(missing) > 0) {
    cat("MISSING:", paste(missing, collapse = ","), "\n")
    quit(status = 1)
  }
  cat("OK\n")
' "$R_LIB_DIR" 2>&1) || {
    echo "error: [env-preflight] R library at $R_LIB_DIR does not resolve spacexr + SPLIT:" >&2
    echo "$R_CHECK_OUT" | sed 's/^/       /' >&2
    echo "       Install them there (docs/installation.md § R side), or point" >&2
    echo "       --r-lib-dir at the right location and re-run write-env-config.sh." >&2
    exit 6
}
echo "[env-preflight] R library resolves spacexr + SPLIT: $R_CHECK_OUT"
echo "[env-preflight] all checks passed."
