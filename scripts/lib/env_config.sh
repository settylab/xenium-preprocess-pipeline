#!/usr/bin/env bash
# scripts/lib/env_config.sh — load + validate scripts/env.local.conf.
#
# Meant to be SOURCED, not executed: `source ".../lib/env_config.sh"`.
# Both the sbatch stubs (submit_{xenium-preprocess,ref-build,rctd-split}.sbatch)
# and the driver (submit_workflow.sh, for its submit-time preflight) source
# this so there is exactly one place that knows the config's shape.
#
# Rationale for why this file exists at all: the two failures that hit
# this pipeline for real (missing spacexr's R library; micromamba's
# root-prefix resolving to the wrong place) were BOTH env-propagation
# failures — resolving a `$HOME`-relative path *inside the job*, where
# `$HOME` can be a different instance than the submitting shell's. A
# file on shared storage, read for its literal (already-absolute)
# values, is immune to that class of bug; `sbatch --export=ALL` forwards
# env VAR values faithfully but does nothing to stop `$HOME` itself
# from resolving differently once the job's own shell starts — so
# anything that re-derives a path from `$HOME` at job time is still
# exposed even with `--export=ALL`. Resolving once, at install time,
# and reading the resolved absolute values back is what actually fixes
# it.
#
# On success, exports:
#   MICROMAMBA_BIN      absolute path to the micromamba binary
#   XENIUM_ENV_PREFIX   absolute path to the micromamba env PREFIX (not name)
#   R_LIB_DIR           absolute path to the R library holding spacexr/SPLIT
#   R_MODULE            Lmod module name (e.g. fhR/4.4.1-foss-2023b)
#
# Fails loud (exit/return 1) with a "run the install step" message if
# scripts/env.local.conf is missing, unreadable, or missing any key.

_env_config_fail() {
    echo "error: $1" >&2
    echo "       Run scripts/write-env-config.sh to (re)generate "\
"scripts/env.local.conf — see docs/installation.md." >&2
    return 1 2>/dev/null || exit 1
}

_ENV_CONFIG_LIB_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# XENIUM_ENV_LOCAL_CONF overrides the default location — used by the
# test suite to point at a fixture conf instead of the real
# scripts/env.local.conf; production callers should never need to set it.
ENV_LOCAL_CONF="${XENIUM_ENV_LOCAL_CONF:-$_ENV_CONFIG_LIB_DIR/../env.local.conf}"

if [[ ! -f "$ENV_LOCAL_CONF" ]]; then
    _env_config_fail "scripts/env.local.conf not found (looked at $ENV_LOCAL_CONF)."
fi

# shellcheck source=/dev/null
source "$ENV_LOCAL_CONF"

for _var in MICROMAMBA_BIN XENIUM_ENV_PREFIX R_LIB_DIR R_MODULE; do
    if [[ -z "${!_var:-}" ]]; then
        _env_config_fail "scripts/env.local.conf is missing or has an empty '$_var'."
    fi
done
unset _var

if [[ ! -x "$MICROMAMBA_BIN" ]]; then
    _env_config_fail "MICROMAMBA_BIN in scripts/env.local.conf is not an executable file: $MICROMAMBA_BIN"
fi

if [[ ! -d "$XENIUM_ENV_PREFIX" ]]; then
    _env_config_fail "XENIUM_ENV_PREFIX in scripts/env.local.conf does not exist: $XENIUM_ENV_PREFIX"
fi

if [[ ! -d "$R_LIB_DIR" ]]; then
    _env_config_fail "R_LIB_DIR in scripts/env.local.conf does not exist: $R_LIB_DIR"
fi

export MICROMAMBA_BIN XENIUM_ENV_PREFIX R_LIB_DIR R_MODULE
