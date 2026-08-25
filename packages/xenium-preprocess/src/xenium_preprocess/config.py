"""YAML config loading, merging, and validation for xenium-preprocess.

The canonical default config lives at `<package>/config/default.yaml`
(top-level of the source tree). `load_default()` returns it as a dict;
`load_yaml(path)` returns any YAML file as a dict; `deep_update(base, override)`
does a recursive merge. `validate(cfg)` raises with a readable list of
missing required top-level keys.
"""
from __future__ import annotations

from pathlib import Path

import yaml

# Package layout:
#   src/xenium_preprocess/config.py       <- this file
#   src/xenium_preprocess/
#   config/default.yaml                   <- default config (repo top-level)
#
# From this file, walk up: parents[0]=xenium_preprocess, [1]=src, [2]=repo root.
_REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH = _REPO_ROOT / "config" / "default.yaml"

REQUIRED_KEYS = ("sample_id", "proseg_dir", "output_root")

# Every stage this pipeline knows how to dispatch. `preprocess` is
# retained here for reversibility (2026-08-11) but is NOT in the
# default list — re-enable by explicit `--stages ... preprocess ...`.
# `celltype_writeback` is on the roadmap but NOT part of the xenium-preprocess
# implementation worker's scope; it lands in a follow-up worker so
# `VALID_STAGES` here does NOT list it yet — adding it would falsely
# advertise a stage that has no dispatch block in pipeline.run().
VALID_STAGES = (
    "proseg_to_anndata",
    "enrich_xenium_id",
    "qc_filter",
    "xenium_ranger_to_anndata",
    "preprocess",
    "split_prep",
    "rctd_prep",
)

# Stages that run when `--stages` is not passed on the CLI (the new,
# post-refactor default). No `preprocess` — no PCA/UMAP/Leiden/celltype
# in xenium-preprocess per the roadmap. Order matters: it's the on-disk dispatch
# order in `pipeline.run()`.
#
# History of the split_prep / rctd_prep defaults:
#   - 2026-08-11 (internal issue review): BOTH stages
#     dropped as a bundle. Interpreted the user's "remove all 10x bundle"
#     ask as targeting the RDS as well.
#   - 2026-08-12 (internal issue review): rctd_prep
#     PARTIALLY restored. the user clarified she still needs the R
#     `test_object.rds` — only the mtx/features/barcodes side-artifact
#     under `split_prep/` was meant to disappear. `rctd_prep` is back
#     in defaults; `split_prep` stays out as an opt-in intermediate.
#     The pipeline auto-produces the mtx bundle as a transient
#     prerequisite for rctd_prep and deletes it after the RDS lands.
#
# `enrich_xenium_id` runs AFTER `xenium_ranger_to_anndata`: the reversed
# direction (proseg cell_id → xenium-ranger h5ad, 2026-08-11) needs
# both source h5ads on disk before it can NN-map.
DEFAULT_STAGES = (
    "proseg_to_anndata",
    "qc_filter",
    "xenium_ranger_to_anndata",
    "enrich_xenium_id",
    "rctd_prep",
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
