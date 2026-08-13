"""Run-scoped layout tests (TracyY123-nexus#26 comment 5251080220).

Locked layout:
    <output_root>/<S>/<S>_<run_id>/{spatial_adata,rctd}/
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def _synth_proseg_dir(tmp_path: Path, n_cells: int = 8, n_genes: int = 5) -> Path:
    rng = np.random.default_rng(7)
    proseg_dir = tmp_path / "proseg"
    proseg_dir.mkdir()
    gene_names = [f"GENE{i:03d}" for i in range(n_genes)]
    expected = rng.uniform(0.0, 3.0, size=(n_cells, n_genes)).astype(np.float32)
    pd.DataFrame(expected, columns=gene_names).to_csv(
        proseg_dir / "expected-counts.csv", index=False
    )
    pd.DataFrame({
        "cell": [f"cell_{i:04d}" for i in range(n_cells)],
        "centroid_x": rng.uniform(0.0, 1000.0, size=n_cells),
        "centroid_y": rng.uniform(0.0, 1000.0, size=n_cells),
    }).to_csv(proseg_dir / "cell-metadata.csv", index=False)
    return proseg_dir


def test_run_dir_shape():
    """`<output_root>/<S>/<S>_<run_id>/` — anchor for every other path."""
    from xenium_preprocess._internal.layout import (
        rctd_dir, run_dir, spatial_adata_dir,
    )
    root = Path("/tmp/out")
    assert run_dir(root, "MH10", "42") == Path("/tmp/out/MH10/MH10_42")
    assert spatial_adata_dir(root, "MH10", "42") == \
        Path("/tmp/out/MH10/MH10_42/spatial_adata")
    assert rctd_dir(root, "MH10", "42") == Path("/tmp/out/MH10/MH10_42/rctd")


def test_spatial_adata_path_basenames():
    from xenium_preprocess._internal.layout import spatial_adata_path
    root = Path("/tmp/out")
    assert spatial_adata_path(root, "MH10", "42", "proseg_raw") == \
        Path("/tmp/out/MH10/MH10_42/spatial_adata/MH10_proseg_raw.h5ad")
    assert spatial_adata_path(root, "MH10", "42", "xenium_ranger") == \
        Path("/tmp/out/MH10/MH10_42/spatial_adata/MH10_xenium_ranger.h5ad")
    assert spatial_adata_path(root, "MH10", "42", "proseg_purified") == \
        Path("/tmp/out/MH10/MH10_42/spatial_adata/MH10_proseg_purified.h5ad")


def test_rctd_path_basenames():
    from xenium_preprocess._internal.layout import rctd_path
    root = Path("/tmp/out")
    assert rctd_path(root, "MH10", "42", "test_object") == \
        Path("/tmp/out/MH10/MH10_42/rctd/MH10_test_object.rds")
    assert rctd_path(root, "MH10", "42", "reference") == \
        Path("/tmp/out/MH10/MH10_42/rctd/MH10_reference.rds")


def test_proseg_to_anndata_writes_new_layout(tmp_path: Path):
    try:
        import anndata  # noqa: F401
        import scipy  # noqa: F401
    except ImportError as e:
        pytest.skip(f"missing dep: {e.name}")

    from xenium_preprocess.stages.proseg_to_anndata import run_proseg_to_anndata

    proseg_dir = _synth_proseg_dir(tmp_path)
    out_root = tmp_path / "out"
    sample = "SYNTH"
    run_id = "test_run_1"

    new_path = run_proseg_to_anndata(
        sample_id=sample,
        run_id=run_id,
        proseg_dir=proseg_dir,
        output_root=out_root,
        count_matrix_path=None,
        cell_metadata_path=None,
        maxpost_matrix_path=None,
        count_matrix_glob="expected-counts*",
        cell_metadata_glob="cell-metadata*",
        maxpost_matrix_glob=None,
        centroid_x_col="centroid_x",
        centroid_y_col="centroid_y",
        cell_id_col="cell",
        force_rerun=False,
    )
    expected = (
        out_root / sample / f"{sample}_{run_id}" / "spatial_adata"
        / f"{sample}_proseg_raw.h5ad"
    )
    assert new_path == expected
    assert new_path.exists()
    assert not new_path.is_symlink()  # No more legacy symlinks in the new layout.


# ---------------------------------------------------------------------------
# atomic_write_h5ad sanitizer coverage (settylab/TracyY123-nexus#26 comments
# 5259204580 and 5259518318). anndata's h5ad writer has no registered handler
# for pandas.arrays.ArrowStringArray on the index code path OR on
# `write_categorical`'s `.categories._values` path. The sanitizer must cover
# both, without permanently mutating the caller's AnnData.
# ---------------------------------------------------------------------------


def _arrow_string_env_or_skip():
    """Skip if pandas doesn't produce ArrowStringArray by default — the
    sanitizer's failure mode only exists in that env. On pandas <3.x with
    default settings the frames stay object-dtype and the write path never
    hits the arrow branch."""
    pd = pytest.importorskip("pandas")
    if not pd.get_option("future.infer_string"):
        pytest.skip("future.infer_string=False; ArrowStringArray not produced")


def test_sanitize_string_column_and_index():
    _arrow_string_env_or_skip()
    import pandas as pd

    from xenium_preprocess._internal.layout import _sanitize_frame_for_h5ad

    df = pd.DataFrame(
        {"seg_method": pd.array(["proseg", "proseg", "manual"], dtype="string")},
        index=pd.Index(
            pd.array(["cell_A", "cell_B", "cell_C"], dtype="string"), name="cell_d"
        ),
    )
    assert isinstance(df["seg_method"].dtype, pd.StringDtype)
    assert isinstance(df.index.dtype, pd.StringDtype)

    fixed = _sanitize_frame_for_h5ad(df)
    assert fixed is not df  # copy, not in-place
    assert fixed["seg_method"].dtype == object
    assert fixed.index.dtype == object
    assert fixed.index.name == "cell_d"
    assert fixed["seg_method"].tolist() == ["proseg", "proseg", "manual"]
    assert fixed.index.tolist() == ["cell_A", "cell_B", "cell_C"]


def test_sanitize_categorical_with_arrow_string_categories():
    """Regression for comment 5259518318 — /obs/fov categorical whose
    categories are backed by ArrowStringArray."""
    _arrow_string_env_or_skip()
    import pandas as pd

    from xenium_preprocess._internal.layout import _sanitize_frame_for_h5ad

    fov = pd.array(["A", "B", "A", "C", "B"], dtype="string")
    df = pd.DataFrame({"fov": pd.Series(fov).astype("category")})
    assert isinstance(df["fov"].dtype, pd.CategoricalDtype)
    assert isinstance(df["fov"].cat.categories.dtype, pd.StringDtype)

    fixed = _sanitize_frame_for_h5ad(df)
    assert fixed is not df
    assert isinstance(fixed["fov"].dtype, pd.CategoricalDtype)  # still categorical
    assert fixed["fov"].cat.categories.dtype == object
    # Values + codes preserved.
    assert fixed["fov"].tolist() == ["A", "B", "A", "C", "B"]
    assert sorted(fixed["fov"].cat.categories.tolist()) == ["A", "B", "C"]


def test_atomic_write_h5ad_roundtrips_arrow_and_categorical_obs(tmp_path):
    """End-to-end: atomic_write_h5ad succeeds on obs holding both a plain
    ArrowString column and a Categorical-with-arrow-string-categories
    column, values round-trip, and caller's AnnData is left untouched."""
    _arrow_string_env_or_skip()
    anndata = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    import pandas as pd

    from xenium_preprocess._internal.layout import atomic_write_h5ad

    n = 6
    obs = pd.DataFrame({
        "fov": pd.Series(
            pd.array(["A", "B", "A", "C", "B", "A"], dtype="string"),
        ).astype("category"),
        "seg_method": pd.array(
            ["proseg"] * n, dtype="string",
        ),
    })
    ad = anndata.AnnData(
        X=np.zeros((n, 1), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["g0"]),
    )
    # Snapshot pre-write dtypes so we can assert restoration.
    pre_fov_cat_type = type(ad.obs["fov"].cat.categories._values)
    pre_seg_dtype = ad.obs["seg_method"].dtype

    out = tmp_path / "cat_arrow_roundtrip.h5ad"
    atomic_write_h5ad(ad, out)
    assert out.exists()

    back = anndata.read_h5ad(out)
    assert isinstance(back.obs["fov"].dtype, pd.CategoricalDtype)
    assert back.obs["fov"].tolist() == ["A", "B", "A", "C", "B", "A"]
    assert back.obs["seg_method"].tolist() == ["proseg"] * n

    # Caller's AnnData restored to its original arrow-backed state.
    assert type(ad.obs["fov"].cat.categories._values) is pre_fov_cat_type
    assert ad.obs["seg_method"].dtype == pre_seg_dtype


def test_atomic_write_h5ad_handles_strings_with_missing_values(tmp_path):
    """Regression for comment 5259730464 — Tracy's third-round failure.

    `atomic_write_h5ad` must succeed on obs with a plain ArrowString column
    that contains missing values AND is NOT pre-categorized. The failure
    mode: anndata.write_h5ad runs adata.strings_to_categoricals() on entry
    (`convert_strings_to_categoricals=True` default), which — under pandas
    3.x with future.infer_string=True — categorizes the ArrowString column
    into a Categorical whose categories are ArrowString-backed, defeating
    our sanitize that ran before write. The fix runs strings_to_categoricals
    ourselves BEFORE sanitize so sanitize catches the arrow-backed
    categories.
    """
    _arrow_string_env_or_skip()
    anndata = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    import pandas as pd

    from xenium_preprocess._internal.layout import atomic_write_h5ad

    n = 8
    fov_values = ["G7", "L7", None, "K8", None, "O4", "P7", "G7"]
    obs = pd.DataFrame({
        "fov": pd.array(fov_values, dtype="string"),
    })
    ad = anndata.AnnData(
        X=np.zeros((n, 1), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["g0"]),
    )
    # Pre-condition: fov is a plain (non-categorical) ArrowString column
    # with real missing values — the exact shape Tracy's cell-metadata CSV
    # produces from pd.read_csv under future.infer_string=True.
    assert isinstance(ad.obs["fov"].dtype, pd.StringDtype)
    assert ad.obs["fov"].isna().sum() == 2

    out = tmp_path / "strings_with_na.h5ad"
    atomic_write_h5ad(ad, out)
    assert out.exists()

    back = anndata.read_h5ad(out)
    # Missing preserved as NaN; string values preserved.
    assert back.obs["fov"].isna().sum() == 2
    got = back.obs["fov"].astype(object).where(back.obs["fov"].notna(), None).tolist()
    assert got == fov_values

    # Caller's AnnData: fov dtype restored (not permanently categorized).
    assert isinstance(ad.obs["fov"].dtype, pd.StringDtype)


def test_sanitize_is_noop_when_no_string_ext():
    """Frames with no arrow-string columns/index and no arrow-backed
    categoricals should short-circuit (return-input identity) so we don't
    pay a copy on the common path."""
    import numpy as np
    import pandas as pd

    from xenium_preprocess._internal.layout import _sanitize_frame_for_h5ad

    df = pd.DataFrame(
        {
            "x": np.arange(3, dtype=np.int64),
            "y": np.arange(3.0, dtype=np.float64),
        },
        index=pd.Index([0, 1, 2], name="idx"),
    )
    fixed = _sanitize_frame_for_h5ad(df)
    assert fixed is df  # no copy on the fast path
