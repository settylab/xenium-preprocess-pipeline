"""exploration — downstream / optional exploratory analyses for Xenium / Proseg
preprocessed data.

Conceptual boundary (vs `xenium-preprocess`):

- `xenium-preprocess`: stable, general-purpose preprocessing that should
  apply broadly across Xenium / Proseg datasets. Produces the canonical
  processed `.h5ad`.
- `exploration`: downstream, question-specific or optional analyses that
  are NOT required for producing the canonical preprocessed dataset.

Flow:

    raw data -> xenium-preprocess -> canonical processed AnnData
             -> exploration -> downstream analyses

See `README.md` for the full invariants (additive, non-mutating, loud,
zero cross-package deps) and the extension recipe.

Public surface:
- __version__ — package version.
- annotate_regions_and_tumor_distance — region + signed tumor-distance
  annotation (Tracy #31).
"""
from __future__ import annotations

__version__ = "0.1.0"

from exploration.annotate import annotate_regions_and_tumor_distance

__all__ = ["__version__", "annotate_regions_and_tumor_distance"]
