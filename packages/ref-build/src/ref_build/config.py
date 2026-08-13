"""YAML config loading, merging, and validation for ref-build.

The canonical default config lives at `<repo>/config/default.yaml`
(top-level of the source tree). `load_default()` returns it as a dict;
`load_yaml(path)` returns any YAML file as a dict; `deep_update(base, override)`
does a recursive merge. `validate(cfg)` raises with a readable list of
missing required top-level keys.
"""
from __future__ import annotations

from pathlib import Path

import yaml

# Package layout:
#   src/ref_build/config.py       <- this file
#   src/ref_build/
#   config/default.yaml           <- default config (repo top-level)
#
# From this file, walk up: parents[0]=ref_build, [1]=src, [2]=repo root.
_REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH = _REPO_ROOT / "config" / "default.yaml"

REQUIRED_KEYS = (
    "sample_id",
    "primary_h5ad",
    "output_root",
    "celltype_marker_json",
)

VALID_STAGES = (
    "load_primary_and_donors",
    "census",
    "assemble",
    "export_mtx",
    "rctd_reference_build",
)


def deep_update(base: dict, override: dict) -> dict:
    """Recursive dict merge — override wins."""
    out = dict(base)
    for k, v in override.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = deep_update(out[k], v)
        else:
            out[k] = v
    return out


def load_yaml(path: Path) -> dict:
    """Read a YAML file, return `{}` if the file is empty."""
    with open(path) as f:
        return yaml.safe_load(f) or {}


def load_default() -> dict:
    """Load the package-shipped default config."""
    return load_yaml(DEFAULT_CONFIG_PATH)


def validate(cfg: dict) -> None:
    """Raise SystemExit if any required top-level key is missing/null."""
    missing = [k for k in REQUIRED_KEYS if not cfg.get(k)]
    if missing:
        raise SystemExit(f"missing required config keys: {missing}")
