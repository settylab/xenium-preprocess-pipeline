# exploration

Downstream, question-specific or optional analyses that run **after** the
canonical preprocessed AnnData has been produced by `xenium-preprocess`.
This package is explicitly opt-in: no stage of `xenium-preprocess` depends
on anything here, and the canonical processed `.h5ad` is produced without
any `exploration` function ever being invoked.

## Conceptual boundary

- **`xenium-preprocess`** — stable, general-purpose preprocessing that
  should apply broadly across Xenium / Proseg datasets. Produces the
  canonical processed `.h5ad`.
- **`exploration`** — downstream, question-specific or optional analyses
  that are **not required** for producing the canonical preprocessed
  dataset. May include reusable functions for exploratory spatial
  analyses, QC, scoring, visualization, cross-sample comparisons, and
  hypothesis-driven analyses.

Flow:

```
raw data  →  xenium-preprocess  →  canonical processed AnnData  →  exploration  →  downstream analyses
```

## Current public surface

| Function | Purpose |
|---|---|
| `exploration.annotate.annotate_regions_and_tumor_distance(adata, region_csv, ...)` | Add three strictly-additive `.obs` columns — `selected_region` (all cells), `tumor_distance_bin` + `tumor_boundary_distance_um` (hepatocytes only) — from a user-authored region-polygon CSV. Returns a new AnnData; the input is never mutated. Tracy `#31`. |

## Design invariants (shared across all `exploration` functions)

1. **Strictly additive** — functions in this package augment `.obs` /
   `.uns` / `.obsm` but never mutate pre-existing columns, never touch
   `.X` / `.layers` / `.var`, never drop cells.
2. **Input is never mutated** — functions take an AnnData, return a new
   AnnData (`adata.copy()`). The caller writes the result wherever they
   want; the canonical preprocessed `.h5ad` on disk is never overwritten
   by this package.
3. **Loud failure on malformed inputs** — downstream analyses should
   fail fast rather than silently produce all-NaN annotations. CSV
   paths that don't exist, polygon files with zero polygons, obs
   columns that are missing but required by a function all raise
   `RuntimeError` with an actionable message.
4. **Zero cross-package import dependency** — `exploration` is
   installable and usable without `xenium-preprocess` on the Python
   path. (The ONLY shared helper, the `log()` function, is inlined
   under `exploration._internal.logging`.)

## Install

```
pip install -e packages/exploration/
```

or add `exploration` to your conda env from `packages/exploration/environment.yml`.

## Usage

```python
import anndata
from exploration import annotate_regions_and_tumor_distance

adata = anndata.read_h5ad("MH6_BR_2328_proseg_purified.h5ad")
annotated = annotate_regions_and_tumor_distance(adata, "MH6_Annotations.csv")
annotated.write_h5ad("MH6_BR_2328_proseg_purified_with_regions.h5ad")  # NOT the canonical path
```

## Extending

Adding a new exploratory function:

1. Create `src/exploration/<topic>.py` (or a subpackage if the topic
   warrants multiple files).
2. Honour the four invariants above (additive, non-mutating, loud, no
   cross-package deps).
3. Add a corresponding `tests/test_<topic>.py` with unit tests + a
   preservation test that verifies `.X` / `.layers` / `.var` / `.obsm`
   / pre-existing `.obs` columns are unchanged.
4. Re-export the public symbol from `src/exploration/__init__.py` and
   add a row to the "Current public surface" table above.
