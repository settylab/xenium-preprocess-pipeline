#!/usr/bin/env bash
# scripts/uv-pip-install.sh — run `uv pip install` with an isolated
# UV_CACHE_DIR, self-verifying that ~/.cache/uv stayed untouched.
#
# Same shape as scripts/create-env.sh's `CONDA_PKGS_DIRS` override:
# `uv pip install` always writes its content-addressed cache to
# `~/.cache/uv` unless `UV_CACHE_DIR` is set — MAMBA_ROOT_PREFIX being
# set (docs/installation.md § 2's isolation guidance) does nothing to
# stop that on its own, since uv doesn't read micromamba's config at
# all. Nothing errors when this happens; a "fully isolated" install
# can still quietly leave files behind in the user's home directory.
# This wraps `uv pip install` with an explicit `UV_CACHE_DIR` override
# (when MAMBA_ROOT_PREFIX is set) and asserts, after the fact, that
# `~/.cache/uv` was not touched.
#
# If MAMBA_ROOT_PREFIX is NOT set, this is a thin passthrough to `uv
# pip install` (no isolation was requested, so there's nothing to
# isolate or assert).
#
# Usage:
#   scripts/uv-pip-install.sh -e packages/xenium-preprocess
#   scripts/uv-pip-install.sh -e "packages/xenium-preprocess[test]"
#   MAMBA_ROOT_PREFIX=/abs/isolated/root scripts/uv-pip-install.sh -e packages/ref-build
#
# All arguments are passed through to `uv pip install` verbatim.

set -euo pipefail

if ! command -v uv >/dev/null 2>&1; then
    echo "error: [uv-pip-install] uv not found on PATH." >&2
    exit 4
fi

HOME_UV_CACHE="$HOME/.cache/uv"

if [[ -z "${MAMBA_ROOT_PREFIX:-}" ]]; then
    echo "[uv-pip-install] MAMBA_ROOT_PREFIX not set — using uv's default cache; nothing to isolate."
    uv pip install "$@"
    exit 0
fi

if [[ ! -d "$MAMBA_ROOT_PREFIX" ]]; then
    mkdir -p "$MAMBA_ROOT_PREFIX"
fi
MAMBA_ROOT_PREFIX=$(cd "$MAMBA_ROOT_PREFIX" && pwd)
export UV_CACHE_DIR="$MAMBA_ROOT_PREFIX/uv-cache"
mkdir -p "$UV_CACHE_DIR"

echo "[uv-pip-install] UV_CACHE_DIR=$UV_CACHE_DIR"

MARKER=$(mktemp)
trap 'rm -f "$MARKER"' EXIT

uv pip install "$@"

if [[ -d "$HOME_UV_CACHE" ]]; then
    TOUCHED=$(find "$HOME_UV_CACHE" -newer "$MARKER" 2>/dev/null || true)
    if [[ -n "$TOUCHED" ]]; then
        echo "error: [uv-pip-install] $HOME_UV_CACHE was modified during uv pip install" >&2
        echo "       despite UV_CACHE_DIR=$UV_CACHE_DIR — isolation did NOT hold." >&2
        echo "       Files touched:" >&2
        echo "$TOUCHED" | sed 's/^/       /' >&2
        exit 1
    fi
fi
echo "[uv-pip-install] verified: $HOME_UV_CACHE was not modified. Isolation held."
