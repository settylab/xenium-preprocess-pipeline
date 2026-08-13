"""YAML config loading, merging, and validation for rctd-split.

The canonical default config lives at `<package>/config/default.yaml`
(top-level of the source tree). `load_default()` returns it as a dict;
`load_yaml(path)` returns any YAML file as a dict; `deep_update(base, override)`
does a recursive merge. `validate(cfg)` raises with a readable list of
missing required top-level keys.
"""
from __future__ import annotations

from pathlib import Path

import yaml

_REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH = _REPO_ROOT / "config" / "default.yaml"

REQUIRED_KEYS = ("sample_id", "output_root")
# `test_object` and `reference_rds` are OPTIONAL at the config layer:
# when unset, `rctd_split.pipeline.run` auto-derives them from the
# run-folder layout convention
#     <output_root>/<sample_id>/<sample_id>_<run_id>/rctd/<sample_id>_test_object.rds
#     <output_root>/<sample_id>/<sample_id>_<run_id>/rctd/<sample_id>_reference.rds
# and only fails if the derived file does not exist on disk.

VALID_STAGES = (
    "rctd_run",
    "split_purify",
    "export_mtx",
    "mtx_to_h5ad",
    "filter_status",
    "postprocess",
    "writeback_to_step1_raw",
    "celltype_writeback",
    "qc_report",
)

# Default stages: everything except the two writeback stages, which
# depend on step-1 outputs and are opt-in.
DEFAULT_STAGES = (
    "rctd_run",
    "split_purify",
    "export_mtx",
    "mtx_to_h5ad",
    "filter_status",
    "postprocess",
    "writeback_to_step1_raw",
    "celltype_writeback",
    "qc_report",
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
