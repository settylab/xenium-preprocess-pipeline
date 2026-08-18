"""rctd-split — rctd-split of the Xenium spatial-data preprocessing pipeline.

Public surface:
- __version__ — package version.
- rctd_split.cli.main — CLI entry point.
- rctd_split.config — YAML load + merge + validate.
- rctd_split.pipeline.run — stage-loop orchestrator.
- rctd_split.stages.* — individual stages (each exposes a `run_<stage>` fn).
"""
from __future__ import annotations

__version__ = "0.1.0"
