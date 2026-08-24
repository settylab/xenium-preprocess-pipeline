#!/usr/bin/env bash
# scripts/check-python-env.sh — asserts that `python`/`pytest` on PATH
# actually resolve INSIDE the target micromamba env, and that no stray
# PYTHONPATH entries are shadowing it.
#
# The one symptom that has bitten this pipeline for real: `ml
# fhR/...` (docs/installation.md § R side) leaking PATH *and*
# PYTHONPATH into the calling shell. `micromamba activate` fixes PATH
# but never touches PYTHONPATH, so even the venv's own `pytest`
# binary can end up importing a Lmod-managed `_pytest` package instead
# of the venv's — and the resulting failure (`ModuleNotFoundError: No
# module named 'yaml'`, or an `ImportError` importing
# `_pytest.config`) names something that is not actually broken,
# sending you to debug the wrong thing. Run this immediately before
# `pytest` (step 5) so a shadowed env fails loudly and accurately
# instead of silently.
#
# Usage: scripts/check-python-env.sh [--env-name xenium]

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ENV_NAME="xenium"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --env-name) ENV_NAME="$2"; shift 2 ;;
        -h|--help) sed -n '2,18p' "$0"; exit 0 ;;
        *) echo "error: unknown argument: $1" >&2; exit 2 ;;
    esac
done

RECEIPT="$SCRIPT_DIR/.env-prefix-$ENV_NAME"
if [[ ! -f "$RECEIPT" ]]; then
    echo "error: [check-python-env] no receipt at $RECEIPT." >&2
    echo "       Run scripts/create-env.sh -n $ENV_NAME -f environments/xenium.yml first (step 2)." >&2
    exit 3
fi
EXPECTED_PREFIX=$(<"$RECEIPT")

if ! command -v python >/dev/null 2>&1; then
    echo "error: [check-python-env] no 'python' on PATH. Did you run 'micromamba activate $ENV_NAME'?" >&2
    exit 4
fi

CHECK_ENV_NAME="$ENV_NAME" EXPECTED_ENV_PREFIX="$EXPECTED_PREFIX" python - <<'PYEOF'
import os
import sys

expected = os.environ["EXPECTED_ENV_PREFIX"].rstrip("/")
env_name = os.environ["CHECK_ENV_NAME"]
actual = sys.prefix.rstrip("/")

problems = []
if actual != expected:
    problems.append(f"sys.prefix is {actual!r}, expected {expected!r} (python: {sys.executable})")

try:
    import pytest
    pytest_file = getattr(pytest, "__file__", "") or ""
    if not pytest_file.startswith(expected):
        problems.append(f"pytest module resolved from {pytest_file!r}, outside the env")
except ImportError as exc:
    problems.append(f"cannot import pytest at all: {exc}")

if problems:
    print("error: [check-python-env] wrong Python/pytest resolved for this shell.", file=sys.stderr)
    for p in problems:
        print(f"       - {p}", file=sys.stderr)
    pythonpath = os.environ.get("PYTHONPATH", "")
    if pythonpath:
        print(f"       PYTHONPATH is currently set to:\n         {pythonpath}", file=sys.stderr)
        print("       This is almost always 'ml <module>' (Lmod) leaking PATH/PYTHONPATH into", file=sys.stderr)
        print("       this shell — micromamba activate fixes PATH but never touches", file=sys.stderr)
        print("       PYTHONPATH. Fix:", file=sys.stderr)
        print(f"         micromamba activate {env_name} && unset PYTHONPATH", file=sys.stderr)
        print("       Then re-run this check. See docs/installation.md § R side for why", file=sys.stderr)
        print("       running 'ml' inside a subshell avoids this in the first place.", file=sys.stderr)
    else:
        print("       PYTHONPATH is not set, so this looks like a different problem (env not", file=sys.stderr)
        print("       activated, or a stale PATH entry from something else). Re-run:", file=sys.stderr)
        print(f"         micromamba activate {env_name}", file=sys.stderr)
    sys.exit(1)

print(f"[check-python-env] OK: python + pytest resolve inside {expected}")
PYEOF
