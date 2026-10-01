"""Reusable region/distance annotation for an AnnData — Tracy #31.

The pipeline accepts a user-authored CSV of manually traced region
polygons (`Selection`, `X`, `Y`, optional `Class`) and augments a
preprocessed AnnData with three strictly additive ``.obs`` columns:

- ``selected_region`` — ALL cells. The ``Selection`` name of the
  CSV polygon containing the cell centroid, or ``not_selected``
  when no polygon encloses the cell. Overlap rule below.
- ``tumor_distance_bin`` — HEPATOCYTES only. One of
  ``intratumoral`` (signed<0), ``0-100um`` (0 ≤ signed < 100),
  ``100-300um`` (100 ≤ signed < 300), ``>300um`` (signed ≥ 300).
  Non-hepatocyte cells are ``NA``.
- ``tumor_boundary_distance_um`` — HEPATOCYTES only. Signed
  Euclidean distance in µm from the cell centroid to the nearest
  tumor-boundary point (negative inside the tumor union, positive
  outside). Non-hepatocyte cells are ``NaN``.

Design invariants (Tracy #31 — "annotation-only"):

- ADDITIVE ONLY: never mutate pre-existing ``.obs`` columns, never
  touch ``.X``/``.layers``/``.var``/``.obsm``, never drop cells.
  If any of the three target columns already exist on the input,
  raise ``RuntimeError`` immediately — the caller must decide.
- Centroids come from ``adata.obsm['spatial']`` (the pipeline's
  stable interface). CSV ``X``/``Y`` and ``obsm['spatial']`` are
  both in µm (no unit conversion).
- Tumor-selection classification uses the F1-patched gate from
  ``scratch/tracy-i31-mh-all-metx-workflow/run_sample.py:411-430``:

    1. Any Selection whose name starts with ``tumor`` (case-insensitive)
       is tumor.
    2. Any Selection whose name starts with ``non-tumor`` (case-insensitive)
       is non-tumor.
    3. Everything else is ambiguous: majority-vote the ``Class`` column;
       the top class is tumor iff ``"tumor" in top and not
       top.startswith("non-tumor") and top != "non-tumor"``. This correctly
       classifies ``non-fibrotic tumor`` as tumor while rejecting
       ``non-tumor``.

- Overlap rule for ``selected_region``: if two or more polygons
  enclose the same cell centroid, the alphabetically-first
  ``Selection`` name wins. Overlap counts are logged (one summary
  line per annotation call).

Reusable interface:

    annotate_regions_and_tumor_distance(adata, region_csv, ...) -> AnnData

Returns a NEW AnnData (``adata.copy()`` with the three columns added
on the copy); the caller writes the result wherever they want. The
original adata is never mutated — guaranteeing Tracy's "original file
UNCHANGED" invariant at the function boundary.
"""
from __future__ import annotations

from io import StringIO
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from xenium_preprocess._internal.logging import log


# Signed-distance bin edges — kept in sync with the AUDITED algorithm in
# scratch/tracy-i31-mh-all-metx-workflow/run_sample.py:442-446.
_DEFAULT_BINS: tuple[tuple[str, float, float], ...] = (
    ("intratumoral", -float("inf"), 0.0),
    ("0-100um",       0.0,          100.0),
    ("100-300um",   100.0,          300.0),
    (">300um",      300.0,          float("inf")),
)

# Column names (configurable by caller but match Tracy's spec by default).
_DEFAULT_SELECTED_REGION_COL = "selected_region"
_DEFAULT_DISTANCE_BIN_COL = "tumor_distance_bin"
_DEFAULT_DISTANCE_UM_COL = "tumor_boundary_distance_um"
_DEFAULT_NOT_SELECTED_LABEL = "not_selected"
_DEFAULT_HEPATOCYTE_COL = "first_type"
_DEFAULT_HEPATOCYTE_VALUE = "liver"
_DEFAULT_OBSM_KEY = "spatial"
_DEFAULT_CLASS_COL = "Class"


def _read_region_csv(path: Path) -> pd.DataFrame:
    """Parse the region CSV: strip ``#``-prefixed metadata lines, keep the rest."""
    with open(path) as fh:
        raw = fh.read().splitlines()
    data_lines = [ln for ln in raw if not ln.startswith("#")]
    if not data_lines:
        raise RuntimeError(f"[annotate] region CSV empty after dropping # lines: {path}")
    df = pd.read_csv(StringIO("\n".join(data_lines)))
    required = {"Selection", "X", "Y"}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(
            f"[annotate] region CSV {path} missing required columns {sorted(missing)}; "
            f"present: {list(df.columns)}"
        )
    return df


def _build_polygons(csv_df: pd.DataFrame) -> dict:
    """One ``shapely.Polygon`` per ``Selection``. Groups with <3 vertices are
    dropped with a WARN log line."""
    from shapely.geometry import Polygon

    polygons: dict[str, "Polygon"] = {}
    dropped: list[tuple[str, int, str]] = []
    for sel, grp in csv_df.groupby("Selection"):
        pts = grp[["X", "Y"]].to_numpy(dtype=float)
        if len(pts) < 3:
            dropped.append((str(sel), len(pts), "too few vertices"))
            continue
        try:
            polygons[str(sel)] = Polygon(pts)
        except Exception as e:  # noqa: BLE001 — surface the shapely error contextually
            dropped.append((str(sel), len(pts), f"Polygon() raised: {e}"))
    if dropped:
        for name, n, why in dropped:
            log(f"[annotate] WARN: dropping Selection {name!r} (n={n}): {why}")
    return polygons


def _classify_tumor_selections(
    polygons: dict,
    csv_df: pd.DataFrame,
    class_column: str,
) -> list[str]:
    """F1-patched tumor-selection classifier.

    Returns a sorted list of ``Selection`` names deemed tumor. Verbatim
    port of ``scratch/tracy-i31-mh-all-metx-workflow/run_sample.py:411-430``
    (the AUDITED version). Logs a summary of the three classes.
    """
    name_tumor = [s for s in polygons if s.lower().startswith("tumor")]
    name_non_tumor = [s for s in polygons if s.lower().startswith("non-tumor")]
    ambiguous = [s for s in polygons if s not in name_tumor and s not in name_non_tumor]

    data_tumor: list[str] = []
    data_non_tumor: list[str] = []
    data_undetermined: list[str] = []
    for sel in ambiguous:
        if class_column not in csv_df.columns:
            data_undetermined.append(sel)
            continue
        cls = csv_df.loc[csv_df["Selection"] == sel, class_column].astype(str)
        if len(cls) == 0:
            data_undetermined.append(sel)
            continue
        top = cls.str.lower().value_counts().index[0].strip()
        # The F1 gate: primary term is 'tumor' AND label does NOT negate it
        # ('non-tumor' / 'non-tumor foo' → non-tumor; 'non-fibrotic tumor' →
        # tumor; 'tumor-adjacent stroma' → still tumor by this gate, matching
        # the audited run_sample.py:428 behaviour).
        if "tumor" in top and not top.startswith("non-tumor") and top != "non-tumor":
            data_tumor.append(sel)
        else:
            data_non_tumor.append(sel)

    tumor_sels = sorted(set(name_tumor) | set(data_tumor))
    log(
        f"[annotate] tumor classification: name-based tumor={sorted(name_tumor)} "
        f"name-based non-tumor={sorted(name_non_tumor)} "
        f"data-driven tumor={sorted(data_tumor)} "
        f"data-driven non-tumor={sorted(data_non_tumor)} "
        f"undetermined (no Class)={sorted(data_undetermined)}"
    )
    return tumor_sels


def _assign_selected_region(
    coords: np.ndarray,
    polygons: dict,
    not_selected_label: str,
) -> np.ndarray:
    """Vectorised ``selected_region`` assignment with alphabetical
    tie-break for overlap. Returns an object-dtype array of length n_cells."""
    import shapely as shp

    n = coords.shape[0]
    labels = np.full(n, not_selected_label, dtype=object)
    counts = np.zeros(n, dtype=np.int32)

    pts = shp.points(coords[:, 0], coords[:, 1])

    # Iterate selections in REVERSE alphabetical order so the earlier
    # name wins ties (later writes overwrite). Deterministic + explicit.
    for sel in sorted(polygons, reverse=True):
        poly = polygons[sel]
        inside = np.asarray(shp.contains(poly, pts), dtype=bool)
        counts[inside] += 1
        labels[inside] = sel

    overlap_cells = int((counts > 1).sum())
    if overlap_cells:
        log(
            f"[annotate] overlap: {overlap_cells:,} cells fell inside >1 polygon; "
            f"assigned the alphabetically-first Selection among the overlapping polygons."
        )
    return labels


def _signed_distance_to_union(
    coords: np.ndarray,
    tumor_union,
) -> np.ndarray:
    """Signed Euclidean distance from each point to the boundary of
    ``tumor_union`` (negative inside, positive outside)."""
    import shapely as shp

    pts = shp.points(coords[:, 0], coords[:, 1])
    boundary = tumor_union.boundary
    inside = np.asarray(shp.contains(tumor_union, pts), dtype=bool)
    d = np.asarray(shp.distance(pts, boundary), dtype=float)
    return np.where(inside, -d, d)


def _bin_signed_distance(
    signed: np.ndarray,
    bins: Sequence[tuple[str, float, float]],
) -> np.ndarray:
    """Map signed distances to bin labels. Returns object-dtype array; any
    value not falling in a bin (e.g. NaN) stays ``None``.

    Convention: a point is in bin ``(label, lo, hi)`` iff ``lo <= s < hi``.
    """
    out = np.full(len(signed), None, dtype=object)
    for label, lo, hi in bins:
        mask = (signed >= lo) & (signed < hi) & np.isfinite(signed)
        out[mask] = label
    return out


def _extract_centroids(adata, obsm_key: str) -> np.ndarray:
    if obsm_key not in adata.obsm:
        raise RuntimeError(
            f"[annotate] adata.obsm[{obsm_key!r}] not found — "
            f"the pipeline requires cell centroids under this key."
        )
    arr = np.asarray(adata.obsm[obsm_key])
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise RuntimeError(
            f"[annotate] adata.obsm[{obsm_key!r}] has unexpected shape {arr.shape} — "
            f"expected (n_cells, >=2)."
        )
    return arr[:, :2].astype(float, copy=False)


def annotate_regions_and_tumor_distance(
    adata,
    region_csv: str | Path,
    *,
    hepatocyte_column: str = _DEFAULT_HEPATOCYTE_COL,
    hepatocyte_value: str = _DEFAULT_HEPATOCYTE_VALUE,
    obsm_key: str = _DEFAULT_OBSM_KEY,
    selected_region_column: str = _DEFAULT_SELECTED_REGION_COL,
    tumor_distance_bin_column: str = _DEFAULT_DISTANCE_BIN_COL,
    tumor_boundary_distance_column: str = _DEFAULT_DISTANCE_UM_COL,
    not_selected_label: str = _DEFAULT_NOT_SELECTED_LABEL,
    class_column: str = _DEFAULT_CLASS_COL,
    distance_bins: Sequence[tuple[str, float, float]] = _DEFAULT_BINS,
):
    """Add Tracy's three region/distance annotations to an AnnData.

    See the module docstring for the semantics of each column. The input
    adata is NEVER mutated — a shallow copy is returned with the three
    new columns written into its ``.obs``. If any of the three target
    columns already exists on the input, a ``RuntimeError`` is raised.
    """
    import anndata as ad  # noqa: F401 — imported for isinstance clarity if callers subclass
    from shapely.ops import unary_union

    region_csv = Path(region_csv)
    log(f"[annotate] region_csv={region_csv}")
    log(
        f"[annotate] adata: n_obs={adata.n_obs:,} n_vars={adata.n_vars:,} "
        f"obsm_key={obsm_key!r} hep_column={hepatocyte_column!r} "
        f"hep_value={hepatocyte_value!r}"
    )

    # ---- 1. Pre-flight: refuse to overwrite existing columns ---------
    collisions = [
        c for c in (
            selected_region_column,
            tumor_distance_bin_column,
            tumor_boundary_distance_column,
        )
        if c in adata.obs.columns
    ]
    if collisions:
        raise RuntimeError(
            f"[annotate] refusing to overwrite pre-existing obs columns {collisions}. "
            f"This stage is strictly additive — delete or rename them first "
            f"if a re-annotation is intended."
        )

    if hepatocyte_column not in adata.obs.columns:
        raise RuntimeError(
            f"[annotate] hepatocyte_column={hepatocyte_column!r} not found on "
            f"adata.obs (present: {list(adata.obs.columns)[:12]}"
            f"{'...' if len(adata.obs.columns) > 12 else ''}). "
            f"The hepatocyte definition is `{hepatocyte_column} == {hepatocyte_value!r}` "
            f"per Tracy #31 — fix the column or pass hepatocyte_column=..."
        )

    # ---- 2. Load centroids + CSV, build polygons ---------------------
    coords = _extract_centroids(adata, obsm_key)
    csv_df = _read_region_csv(region_csv)
    polygons = _build_polygons(csv_df)
    if not polygons:
        raise RuntimeError(
            f"[annotate] region CSV {region_csv} produced zero polygons "
            f"(all Selections had <3 vertices?)."
        )
    log(f"[annotate] built {len(polygons)} polygons: {sorted(polygons)}")

    # ---- 3. selected_region (ALL cells) ------------------------------
    selected = _assign_selected_region(coords, polygons, not_selected_label)
    _log_counts(f"[annotate] selected_region counts: ", pd.Series(selected))

    # ---- 4. Tumor union + signed distance (hepatocytes only) ---------
    tumor_sels = _classify_tumor_selections(polygons, csv_df, class_column)
    if not tumor_sels:
        raise RuntimeError(
            f"[annotate] no tumor Selections resolved from {region_csv}. "
            f"Classifier found no tumor-named polygons AND no ambiguous "
            f"polygons whose Class majority gates as tumor. "
            f"tumor_distance_bin / tumor_boundary_distance_um cannot be computed."
        )
    log(f"[annotate] tumor union over: {tumor_sels}")
    tumor_union = unary_union([polygons[s] for s in tumor_sels])

    hep_mask = adata.obs[hepatocyte_column].astype(str).to_numpy() == str(hepatocyte_value)
    n_hep = int(hep_mask.sum())
    log(f"[annotate] hepatocytes ({hepatocyte_column}={hepatocyte_value!r}): "
        f"{n_hep:,} / {adata.n_obs:,}")

    # Build full-length arrays; non-hep stays NaN / None.
    signed_full = np.full(adata.n_obs, np.nan, dtype=float)
    bin_full = np.full(adata.n_obs, None, dtype=object)

    if n_hep > 0:
        signed_hep = _signed_distance_to_union(coords[hep_mask], tumor_union)
        signed_full[hep_mask] = signed_hep
        bin_full[hep_mask] = _bin_signed_distance(signed_hep, distance_bins)
        log(
            f"[annotate] hep signed-distance (µm): "
            f"min={np.nanmin(signed_hep):.2f} "
            f"median={float(np.nanmedian(signed_hep)):.2f} "
            f"max={np.nanmax(signed_hep):.2f}"
        )
        _log_counts(
            f"[annotate] tumor_distance_bin counts (hepatocytes): ",
            pd.Series(bin_full[hep_mask]),
        )
    else:
        log(
            f"[annotate] WARN: zero hepatocytes under {hepatocyte_column}="
            f"{hepatocyte_value!r} — distance columns will be all-NA."
        )

    # ---- 5. Return a NEW adata with the three columns added ----------
    out = adata.copy()
    out.obs[selected_region_column] = pd.Categorical(selected)
    out.obs[tumor_boundary_distance_column] = signed_full
    out.obs[tumor_distance_bin_column] = pd.Categorical(
        bin_full,
        categories=[b[0] for b in distance_bins],
        ordered=True,
    )
    return out


def _log_counts(prefix: str, s: pd.Series) -> None:
    """Compact counts summary (object / categorical / string tolerated)."""
    try:
        vc = s.value_counts(dropna=False).to_dict()
    except Exception:  # noqa: BLE001
        vc = {}
    log(prefix + ", ".join(f"{k}={v}" for k, v in vc.items()))
