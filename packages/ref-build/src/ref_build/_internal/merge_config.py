"""Shared merged `config.yaml` helper.

Merge invariant (from (internal issue review):
each pipeline reads the current `config.yaml`, writes back
ONLY its own top-level key (`xenium_preprocess:` / `ref_build:` /
`rctd_split:` / `driver:`), and PRESERVES all sibling top-level keys.
This lets the three pipelines run in sequence under a shared `--run-id`
folder without any pipeline clobbering another's section.

Duplicated per repo by design — there is no shared parent package the
three pipelines depend on. This file is a byte-identical copy of
`xenium_preprocess._internal.merge_config`; the test in
`tests/test_merge_config.py` covers the invariant here just as it does
in the sister repo so a regression in ANY one pipeline is caught
locally.
"""
from __future__ import annotations

import os
from pathlib import Path

import yaml


VALID_TOP_LEVEL_KEYS = (
    "driver",
    "xenium_preprocess",
    "ref_build",
    "rctd_split",
)


def merge_config(
    resolved_yaml_path: Path,
    step_key: str,
    step_cfg: dict,
) -> dict:
    """Read the merged `config.yaml`, replace ONLY `step_key`,
    preserve every other top-level key, atomically write back.

    Returns the full merged dict (post-write) so the caller can log
    the sibling keys that survived.

    Refuses `step_key` values outside `VALID_TOP_LEVEL_KEYS` so a
    typo can't silently land a phantom section.
    """
    if step_key not in VALID_TOP_LEVEL_KEYS:
        raise SystemExit(
            f"[merge_config] step_key={step_key!r} is not one of "
            f"{VALID_TOP_LEVEL_KEYS}"
        )
    resolved_yaml_path = Path(resolved_yaml_path)
    if resolved_yaml_path.exists():
        with open(resolved_yaml_path) as f:
            merged = yaml.safe_load(f) or {}
        if not isinstance(merged, dict):
            raise SystemExit(
                f"[merge_config] existing {resolved_yaml_path} is not a "
                f"mapping (got {type(merged).__name__})"
            )
    else:
        merged = {}

    merged[step_key] = step_cfg

    resolved_yaml_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = resolved_yaml_path.with_suffix(resolved_yaml_path.suffix + ".tmp")
    with open(tmp, "w") as f:
        yaml.safe_dump(merged, f, sort_keys=False)
    os.replace(tmp, resolved_yaml_path)
    return merged
