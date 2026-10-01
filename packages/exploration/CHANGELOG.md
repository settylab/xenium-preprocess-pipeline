# Changelog

All notable changes to the `exploration` package will be documented here.
Follows [Keep a Changelog](https://keepachangelog.com/) shape.

## [0.1.0] — 2026-10-01

### Added (Tracy `#31`)

- Initial `exploration` package scaffold — downstream / optional analyses
  that run **after** the canonical `xenium-preprocess` output. See
  `README.md` for the full conceptual boundary.
- `exploration.annotate.annotate_regions_and_tumor_distance` — moved
  verbatim from `xenium_preprocess.annotate` (originally added on
  v0.2.7 of `xenium-preprocess` in commit `9a8530c`). Takes a
  preprocessed AnnData + a region CSV and returns a NEW AnnData with
  three additive `.obs` columns (`selected_region`,
  `tumor_distance_bin`, `tumor_boundary_distance_um`). Hepatocyte-only
  for the two distance columns (`first_type == "liver"`). F1-patched
  tumor-selection classifier (ported verbatim from the audited
  `scratch/tracy-i31-mh-all-metx-workflow/run_sample.py`). Refuses to
  overwrite pre-existing columns; input adata never mutated.
- `exploration._internal.logging.log` — timestamped-flushed stdout helper
  inlined from `xenium_preprocess._internal.logging` so this package has
  zero import-side dependency on `xenium-preprocess`.
- `tests/test_annotate.py` — 10 unit tests covering overlap rule, bin
  edges, hep-only distances, F1 classifier, preservation invariants,
  and parity against the AUDITED reference algorithm.
- `shapely>=2` carried over from the original `xenium-preprocess`
  dependency set — now lives here instead.
