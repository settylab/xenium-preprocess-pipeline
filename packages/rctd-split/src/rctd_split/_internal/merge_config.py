"""Shared merged `config.yaml` helper.

Merge invariant: each pipeline reads the current `config.yaml`, writes
back ONLY its own top-level key (`xenium_preprocess:` / `ref_build:` /
`rctd_split:` / `driver:`), and PRESERVES sibling top-level keys that
are themselves in that semantic set. This lets the three pipelines
run in sequence under a shared `--run-id` folder without clobbering
each other. Unknown top-level keys (e.g. numeric `step1:` / `step3:` /
`step4:` left behind by a pre-semantic-migration run folder) are
dropped so a re-run against an older run dir self-heals rather than
carrying the legacy junk forward.

Duplicated per repo by design — there is no shared parent package the
three pipelines depend on. `ref-build` and `rctd-split` each ship an
identical copy of this file under their own `_internal/`. The test in
`tests/test_merge_config.py` covers the invariant; the sister repos
MUST carry a parallel test so a regression in any one pipeline is
caught locally.
"""
from __future__ import annotations

import os
import sys
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
    preserve every valid sibling top-level key, drop any legacy /
    unknown top-level key, atomically write back.

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
            existing = yaml.safe_load(f) or {}
        if not isinstance(existing, dict):
            raise SystemExit(
                f"[merge_config] existing {resolved_yaml_path} is not a "
                f"mapping (got {type(existing).__name__})"
            )
    else:
        existing = {}

    merged = {k: v for k, v in existing.items() if k in VALID_TOP_LEVEL_KEYS}
    dropped = sorted(set(existing) - set(merged))
    if dropped:
        print(
            f"[merge_config] dropped stale/unknown top-level keys from "
            f"{resolved_yaml_path}: {dropped}",
            file=sys.stderr,
        )

    merged[step_key] = step_cfg

    resolved_yaml_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = resolved_yaml_path.with_suffix(resolved_yaml_path.suffix + ".tmp")
    with open(tmp, "w") as f:
        yaml.safe_dump(merged, f, sort_keys=False)
    os.replace(tmp, resolved_yaml_path)
    return merged
