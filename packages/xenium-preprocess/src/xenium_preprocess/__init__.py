"""xenium-preprocess — step 1 of the Xenium spatial-data preprocessing pipeline.

Default sub-stages executed in order:
    proseg_to_anndata → qc_filter → xenium_ranger_to_anndata → enrich_xenium_id

Legacy 10x-bundle stages (`split_prep`, `rctd_prep`) and the deprecated
`preprocess` stage remain in `VALID_STAGES` for opt-in via `--stages`
but are no longer part of the default flow.
"""
from __future__ import annotations

__version__ = "0.1.0"

__all__ = ["__version__"]
