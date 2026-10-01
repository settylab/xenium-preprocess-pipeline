"""Tests for ``exploration.annotate.annotate_regions_and_tumor_distance``.

Covers:
  - All cells get a ``selected_region`` label; ``not_selected`` fallback.
  - Overlap polygons: alphabetically-first Selection wins; count reported.
  - F1-patched tumor classifier: name-based + Class-majority gate.
  - Hepatocyte-only distance columns; non-hep cells NaN / NA.
  - Signed-distance sign (negative inside, positive outside) + bin edges.
  - Preservation: raises if any target column already exists; refuses to
    mutate input adata (``.X``, ``.obsm``, ``.var``, existing ``.obs``).
  - Parity: our implementation matches the AUDITED
    ``run_sample.py:load_distance_bins`` signed-distance algorithm on
    a synthetic CSV + adata (same primitive, same output).
"""
from __future__ import annotations

import pytest


def _skip_if_missing(*mods):
    import importlib
    for m in mods:
        try:
            importlib.import_module(m)
        except ImportError:
            pytest.skip(f"missing dependency: {m}")


# ---------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------

def _make_adata(coords, first_type):
    """Minimal AnnData with ``obsm['spatial']`` + ``obs['first_type']``."""
    import anndata
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    coords = np.asarray(coords, dtype=float)
    n = coords.shape[0]
    obs = pd.DataFrame(
        {"first_type": np.asarray(first_type, dtype=object)},
        index=pd.Index([f"cell_{i}" for i in range(n)], name="cell_id"),
    )
    var = pd.DataFrame(index=[f"gene_{k}" for k in range(2)])
    X = csr_matrix(np.ones((n, 2), dtype=np.float32))
    ad = anndata.AnnData(X=X, obs=obs, var=var)
    ad.obsm["spatial"] = coords
    return ad


def _write_csv(path, rows, header_lines=None):
    """Write a region CSV with optional ``#`` prefixed header lines.

    `rows` is a list of (Selection, X, Y, Class) tuples.
    """
    header = [
        "#Selection names: synthetic",
        "#Areas (µm^2): irrelevant",
    ] if header_lines is None else list(header_lines)
    with open(path, "w") as fh:
        for ln in header:
            fh.write(ln + "\n")
        fh.write("Selection,X,Y,Class,Color\n")
        for sel, x, y, cls in rows:
            fh.write(f"{sel},{x},{y},{cls},#000000\n")


def _square(sel, cls, x0, y0, side):
    """Return CSV rows for a square polygon with corners (x0,y0)..(x0+side,y0+side)."""
    pts = [(x0, y0), (x0 + side, y0), (x0 + side, y0 + side), (x0, y0 + side)]
    return [(sel, x, y, cls) for x, y in pts]


# ---------------------------------------------------------------------
# Basic contract: selected_region + hep-only distance columns
# ---------------------------------------------------------------------

def test_selected_region_assigns_all_cells(tmp_path):
    _skip_if_missing("anndata", "pandas", "numpy", "scipy", "shapely")
    from exploration.annotate import annotate_regions_and_tumor_distance

    # Two non-overlapping squares (tumor + non-tumor) plus three cells:
    #   cell_0 inside Tumor_1
    #   cell_1 inside Non-tumor_1
    #   cell_2 outside both → not_selected
    rows = _square("Tumor_1", "Unclassified", 0, 0, 100) \
         + _square("Non-tumor_1", "Unclassified", 200, 200, 100)
    csv = tmp_path / "regions.csv"
    _write_csv(csv, rows)

    adata = _make_adata(
        coords=[(50, 50), (250, 250), (500, 500)],
        first_type=["liver", "liver", "liver"],
    )
    out = annotate_regions_and_tumor_distance(adata, csv)
    sr = out.obs["selected_region"].astype(str).tolist()
    assert sr == ["Tumor_1", "Non-tumor_1", "not_selected"]


def test_non_hep_cells_have_na_distance_columns(tmp_path):
    _skip_if_missing("anndata", "pandas", "numpy", "scipy", "shapely")
    import numpy as np
    from exploration.annotate import annotate_regions_and_tumor_distance

    rows = _square("Tumor_1", "Unclassified", 0, 0, 100)
    csv = tmp_path / "regions.csv"
    _write_csv(csv, rows)

    adata = _make_adata(
        coords=[(50, 50), (50, 50), (50, 50), (50, 50)],
        first_type=["liver", "stroma", "immune", "liver"],
    )
    out = annotate_regions_and_tumor_distance(adata, csv)
    dist = out.obs["tumor_boundary_distance_um"].to_numpy()
    assert np.isfinite(dist[0]) and np.isfinite(dist[3])
    assert np.isnan(dist[1]) and np.isnan(dist[2])
    bin_col = out.obs["tumor_distance_bin"]
    # non-hep cells → NA on a Categorical column
    assert bin_col.isna().iloc[1] and bin_col.isna().iloc[2]
    assert not bin_col.isna().iloc[0]


def test_zero_non_hep_cells_have_non_na_distance(tmp_path):
    """Invariant from Tracy's validation checklist — must be exactly 0."""
    _skip_if_missing("anndata", "pandas", "numpy", "scipy", "shapely")
    import numpy as np
    from exploration.annotate import annotate_regions_and_tumor_distance

    rows = _square("Tumor_1", "Unclassified", 0, 0, 100)
    csv = tmp_path / "regions.csv"
    _write_csv(csv, rows)

    adata = _make_adata(
        coords=[(50, 50)] * 5,
        first_type=["liver", "stroma", "immune", "endothelial", "liver"],
    )
    out = annotate_regions_and_tumor_distance(adata, csv)
    non_hep = out.obs["first_type"] != "liver"
    non_hep_dist = out.obs.loc[non_hep, "tumor_boundary_distance_um"]
    n_non_hep_with_dist = int(np.isfinite(non_hep_dist.to_numpy()).sum())
    assert n_non_hep_with_dist == 0


# ---------------------------------------------------------------------
# Signed-distance + bin edges
# ---------------------------------------------------------------------

def test_signed_distance_sign_and_bins(tmp_path):
    _skip_if_missing("anndata", "pandas", "numpy", "scipy", "shapely")
    import numpy as np
    from exploration.annotate import annotate_regions_and_tumor_distance

    # Tumor square (0,0)..(1000,1000). Place cells at known distances:
    #   cell_0 at (500,500)   — inside, distance to boundary = 500 (signed -500 → intratumoral)
    #   cell_1 at (1050,500)  — outside, 50 from boundary          → 0-100um
    #   cell_2 at (1200,500)  — outside, 200 from boundary         → 100-300um
    #   cell_3 at (1500,500)  — outside, 500 from boundary         → >300um
    #   cell_4 at (1100,500)  — outside, 100 from boundary — boundary case ← 100-300um (lo-inclusive)
    rows = _square("Tumor_1", "Unclassified", 0, 0, 1000)
    csv = tmp_path / "regions.csv"
    _write_csv(csv, rows)

    adata = _make_adata(
        coords=[
            (500, 500),
            (1050, 500),
            (1200, 500),
            (1500, 500),
            (1100, 500),
        ],
        first_type=["liver"] * 5,
    )
    out = annotate_regions_and_tumor_distance(adata, csv)
    signed = out.obs["tumor_boundary_distance_um"].to_numpy()
    bins = out.obs["tumor_distance_bin"].astype(str).tolist()

    np.testing.assert_allclose(signed[0], -500.0, atol=1e-6)
    np.testing.assert_allclose(signed[1],  50.0, atol=1e-6)
    np.testing.assert_allclose(signed[2], 200.0, atol=1e-6)
    np.testing.assert_allclose(signed[3], 500.0, atol=1e-6)
    np.testing.assert_allclose(signed[4], 100.0, atol=1e-6)

    assert bins[0] == "intratumoral"
    assert bins[1] == "0-100um"
    assert bins[2] == "100-300um"
    assert bins[3] == ">300um"
    # Lo-inclusive, hi-exclusive: s=100 lands in 100-300um (NOT 0-100um).
    assert bins[4] == "100-300um"


# ---------------------------------------------------------------------
# F1-patched tumor classifier
# ---------------------------------------------------------------------

def test_f1_classifier_name_based(tmp_path):
    _skip_if_missing("anndata", "pandas", "numpy", "scipy", "shapely")
    from exploration.annotate import _build_polygons, _classify_tumor_selections, _read_region_csv

    rows = (_square("Tumor_1",      "Unclassified", 0,   0,   100)
          + _square("Non-tumor_1",  "Unclassified", 200, 0,   100)
          + _square("tumor_2",      "Unclassified", 400, 0,   100))
    csv = tmp_path / "regions.csv"
    _write_csv(csv, rows)

    df = _read_region_csv(csv)
    polys = _build_polygons(df)
    tumor = _classify_tumor_selections(polys, df, "Class")
    assert tumor == ["Tumor_1", "tumor_2"]


def test_f1_classifier_class_majority_includes_non_fibrotic_tumor(tmp_path):
    """The pre-F1 ``"non" not in top`` gate wrongly dropped
    ``non-fibrotic tumor``. The F1 gate keeps it."""
    _skip_if_missing("anndata", "pandas", "numpy", "scipy", "shapely")
    from exploration.annotate import _build_polygons, _classify_tumor_selections, _read_region_csv

    # Ambiguous name — classified via Class majority.
    rows = (_square("Region_A", "non-fibrotic tumor", 0,   0,   100)
          + _square("Region_B", "non-tumor",          200, 0,   100)
          + _square("Region_C", "tumor",              400, 0,   100))
    csv = tmp_path / "regions.csv"
    _write_csv(csv, rows)

    df = _read_region_csv(csv)
    polys = _build_polygons(df)
    tumor = _classify_tumor_selections(polys, df, "Class")
    assert tumor == ["Region_A", "Region_C"]  # Region_B excluded


# ---------------------------------------------------------------------
# Overlap rule for selected_region
# ---------------------------------------------------------------------

def test_overlap_assigns_alphabetically_first_selection(tmp_path):
    _skip_if_missing("anndata", "pandas", "numpy", "scipy", "shapely")
    from exploration.annotate import annotate_regions_and_tumor_distance

    # Three squares overlapping at (50,50):
    #   A, B, Tumor_1 all contain cell_0 → "A" wins (alphabetically first)
    rows = (_square("A",       "Unclassified", 0,   0,   100)
          + _square("B",       "Unclassified", 0,   0,   100)
          + _square("Tumor_1", "Unclassified", 0,   0,   100))
    csv = tmp_path / "regions.csv"
    _write_csv(csv, rows)

    adata = _make_adata(
        coords=[(50, 50), (500, 500)],
        first_type=["liver", "liver"],
    )
    out = annotate_regions_and_tumor_distance(adata, csv)
    sr = out.obs["selected_region"].astype(str).tolist()
    assert sr[0] == "A"
    assert sr[1] == "not_selected"


# ---------------------------------------------------------------------
# Preservation requirements
# ---------------------------------------------------------------------

def test_refuses_to_overwrite_existing_columns(tmp_path):
    _skip_if_missing("anndata", "pandas", "numpy", "scipy", "shapely")
    from exploration.annotate import annotate_regions_and_tumor_distance

    rows = _square("Tumor_1", "Unclassified", 0, 0, 100)
    csv = tmp_path / "regions.csv"
    _write_csv(csv, rows)

    adata = _make_adata([(50, 50)], ["liver"])
    adata.obs["selected_region"] = "pre-existing"
    with pytest.raises(RuntimeError, match="refusing to overwrite"):
        annotate_regions_and_tumor_distance(adata, csv)


def test_does_not_mutate_input_adata(tmp_path):
    _skip_if_missing("anndata", "pandas", "numpy", "scipy", "shapely")
    import numpy as np
    from exploration.annotate import annotate_regions_and_tumor_distance

    rows = _square("Tumor_1", "Unclassified", 0, 0, 100)
    csv = tmp_path / "regions.csv"
    _write_csv(csv, rows)

    adata = _make_adata(
        coords=[(50, 50), (1050, 500)],
        first_type=["liver", "stroma"],
    )
    pre_cols = list(adata.obs.columns)
    pre_x = adata.X.toarray().copy()
    pre_spatial = np.asarray(adata.obsm["spatial"]).copy()
    pre_var = list(adata.var.columns)
    pre_n_obs = adata.n_obs

    out = annotate_regions_and_tumor_distance(adata, csv)

    assert list(adata.obs.columns) == pre_cols
    np.testing.assert_array_equal(adata.X.toarray(), pre_x)
    np.testing.assert_array_equal(np.asarray(adata.obsm["spatial"]), pre_spatial)
    assert list(adata.var.columns) == pre_var
    assert adata.n_obs == pre_n_obs

    # And the output has the three new columns + same cell count.
    assert out.n_obs == pre_n_obs
    assert "selected_region" in out.obs.columns
    assert "tumor_distance_bin" in out.obs.columns
    assert "tumor_boundary_distance_um" in out.obs.columns


# ---------------------------------------------------------------------
# Parity with the AUDITED run_sample.py algorithm (reference impl)
# ---------------------------------------------------------------------

def _reference_signed_distance(coords, csv_df):
    """Verbatim subset of ``run_sample.py:load_distance_bins`` lines 383-446.

    Not imported from scratch/ to keep the test self-contained; the
    algorithm here is the exact snippet the pipeline is specified to
    reproduce.
    """
    import numpy as np
    import shapely as shp
    from shapely.geometry import Polygon
    from shapely.ops import unary_union

    polygons = {}
    for sel, grp in csv_df.groupby("Selection"):
        pts = grp[["X", "Y"]].to_numpy(dtype=float)
        if len(pts) < 3:
            continue
        polygons[sel] = Polygon(pts)

    tumor_sels_name = [s for s in polygons if s.lower().startswith("tumor")]
    ambig = [s for s in polygons
             if not (s.lower().startswith("tumor")
                     or s.lower().startswith("non-tumor"))]
    tumor_sels_data = []
    for sel in ambig:
        cls = csv_df.loc[csv_df["Selection"] == sel, "Class"].astype(str)
        if len(cls) == 0:
            continue
        top = cls.str.lower().value_counts().index[0].strip()
        if "tumor" in top and not top.startswith("non-tumor") and top != "non-tumor":
            tumor_sels_data.append(sel)
    tumor_sels = sorted(set(tumor_sels_name + tumor_sels_data))
    union = unary_union([polygons[s] for s in tumor_sels])
    boundary = union.boundary
    pts = shp.points(coords[:, 0], coords[:, 1])
    inside = np.asarray(shp.contains(union, pts), dtype=bool)
    d = np.asarray(shp.distance(pts, boundary), dtype=float)
    signed = np.where(inside, -d, d)

    bins = np.full(len(signed), None, dtype=object)
    bins[signed < 0] = "intratumoral"
    bins[(signed >= 0) & (signed < 100)] = "0-100um"
    bins[(signed >= 100) & (signed < 300)] = "100-300um"
    bins[signed >= 300] = ">300um"
    return signed, bins


def test_parity_with_reference_algorithm_on_synthetic_multiselection(tmp_path):
    """Synthetic test — build a CSV with Tumor/Non-tumor/ambiguous polygons,
    scatter cells, and verify our implementation returns the exact same
    signed distance + bin as the audited reference algorithm."""
    _skip_if_missing("anndata", "pandas", "numpy", "scipy", "shapely")
    import numpy as np
    from exploration.annotate import annotate_regions_and_tumor_distance, _read_region_csv

    rows = (
        _square("Tumor_1",              "Unclassified", 0,    0,    1000)
      + _square("Non-tumor_1",          "Unclassified", 2000, 0,    1000)
      + _square("tumor_2_lowercase",    "Unclassified", 0,    2000, 500)
      + _square("ambig_fibrotic_tumor", "non-fibrotic tumor", 3000, 3000, 300)
      + _square("ambig_lymph_node",     "lymphocyte",   4000, 4000, 300)
    )
    csv = tmp_path / "regions.csv"
    _write_csv(csv, rows)

    rng = np.random.default_rng(seed=0)
    coords = rng.uniform(-500, 5000, size=(500, 2))

    adata = _make_adata(coords=coords, first_type=["liver"] * len(coords))

    out = annotate_regions_and_tumor_distance(adata, csv)
    our_signed = out.obs["tumor_boundary_distance_um"].to_numpy()
    our_bin = out.obs["tumor_distance_bin"].astype(str).tolist()

    ref_signed, ref_bin = _reference_signed_distance(coords, _read_region_csv(csv))

    np.testing.assert_allclose(our_signed, ref_signed, atol=1e-6)
    assert our_bin == [str(x) if x is not None else "nan" for x in ref_bin]
