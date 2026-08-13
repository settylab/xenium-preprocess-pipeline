"""Run-scoped layout tests for rctd-split (post-2026-08-11 refactor).

Layout: ``<output_root>/<sample>/<sample>_<run_id>/{spatial_adata,rctd,intermediate}``.
"""
from __future__ import annotations

from pathlib import Path

import pytest


def test_run_dir_shape(tmp_path: Path):
    from rctd_split._internal.layout import run_dir, spatial_adata_dir, rctd_dir

    root = tmp_path / "runs"
    assert run_dir(root, "MHTEST", "42") == root / "MHTEST" / "MHTEST_42"
    assert spatial_adata_dir(root, "MHTEST", "42") == root / "MHTEST" / "MHTEST_42" / "spatial_adata"
    assert rctd_dir(root, "MHTEST", "42") == root / "MHTEST" / "MHTEST_42" / "rctd"


def test_spatial_adata_paths(tmp_path: Path):
    from rctd_split._internal.layout import spatial_adata_path

    root = tmp_path / "runs"
    assert spatial_adata_path(root, "MHTEST", "42", "proseg_purified") == (
        root / "MHTEST" / "MHTEST_42" / "spatial_adata" / "MHTEST_proseg_purified.h5ad"
    )
    assert spatial_adata_path(root, "MHTEST", "42", "proseg_raw") == (
        root / "MHTEST" / "MHTEST_42" / "spatial_adata" / "MHTEST_proseg_raw.h5ad"
    )
    assert spatial_adata_path(root, "MHTEST", "42", "xenium_ranger") == (
        root / "MHTEST" / "MHTEST_42" / "spatial_adata" / "MHTEST_xenium_ranger.h5ad"
    )


def test_rctd_paths(tmp_path: Path):
    from rctd_split._internal.layout import rctd_path

    root = tmp_path / "runs"
    # Renamed: rctd_results.rds now carries the sample_id prefix.
    assert rctd_path(root, "MHTEST", "42", "rctd_results") == (
        root / "MHTEST" / "MHTEST_42" / "rctd" / "MHTEST_rctd_results.rds"
    )


def test_intermediate_paths(tmp_path: Path):
    from rctd_split._internal.layout import intermediate_path

    root = tmp_path / "runs"
    prefix = root / "MHTEST" / "MHTEST_42" / "intermediate"
    assert intermediate_path(root, "MHTEST", "42", "unpurified_rds") == (
        prefix / "split" / "MHTEST_unpurified.rds"
    )
    assert intermediate_path(root, "MHTEST", "42", "purified_rds") == (
        prefix / "split" / "MHTEST_purified.rds"
    )
    assert intermediate_path(root, "MHTEST", "42", "unpurified_h5ad") == (
        prefix / "adata" / "MHTEST_step4_unpurified.h5ad"
    )
    assert intermediate_path(root, "MHTEST", "42", "filter_status_csv") == (
        prefix / "adata" / "MHTEST_filter_status.csv"
    )


# ---------------------------------------------------------------------------
# atomic_write_h5ad sanitizer coverage (settylab/TracyY123-nexus#26 comment
# 5260723528). Step 4 reads h5ads produced by step 1 / step 3 whose obs/var
# are backed by pandas.arrays.ArrowStringArray under pandas 3.x
# (`future.infer_string=True`); anndata 0.11's writer registry has no handler
# for that class on either the index code path or the categorical-categories
# path. `atomic_write_h5ad` must sanitize every write, without permanently
# mutating the caller's AnnData.
# ---------------------------------------------------------------------------


def _arrow_string_env_or_skip():
    """Skip if pandas doesn't produce ArrowStringArray by default — the
    sanitizer's failure mode only exists in that env."""
    pd = pytest.importorskip("pandas")
    if not pd.get_option("future.infer_string"):
        pytest.skip("future.infer_string=False; ArrowStringArray not produced")


def test_sanitize_string_column_and_index():
    _arrow_string_env_or_skip()
    import pandas as pd

    from rctd_split._internal.layout import _sanitize_frame_for_h5ad

    df = pd.DataFrame(
        {"seg_method": pd.array(["proseg", "proseg", "manual"], dtype="string")},
        index=pd.Index(
            pd.array(["cell_A", "cell_B", "cell_C"], dtype="string"), name="cell_d"
        ),
    )
    assert isinstance(df["seg_method"].dtype, pd.StringDtype)
    assert isinstance(df.index.dtype, pd.StringDtype)

    fixed = _sanitize_frame_for_h5ad(df)
    assert fixed is not df
    assert fixed["seg_method"].dtype == object
    assert fixed.index.dtype == object
    assert fixed.index.name == "cell_d"
    assert fixed["seg_method"].tolist() == ["proseg", "proseg", "manual"]
    assert fixed.index.tolist() == ["cell_A", "cell_B", "cell_C"]


def test_sanitize_categorical_with_arrow_string_categories():
    _arrow_string_env_or_skip()
    import pandas as pd

    from rctd_split._internal.layout import _sanitize_frame_for_h5ad

    fov = pd.array(["A", "B", "A", "C", "B"], dtype="string")
    df = pd.DataFrame({"fov": pd.Series(fov).astype("category")})
    assert isinstance(df["fov"].dtype, pd.CategoricalDtype)
    assert isinstance(df["fov"].cat.categories.dtype, pd.StringDtype)

    fixed = _sanitize_frame_for_h5ad(df)
    assert fixed is not df
    assert isinstance(fixed["fov"].dtype, pd.CategoricalDtype)
    assert fixed["fov"].cat.categories.dtype == object
    assert fixed["fov"].tolist() == ["A", "B", "A", "C", "B"]
    assert sorted(fixed["fov"].cat.categories.tolist()) == ["A", "B", "C"]


def test_atomic_write_h5ad_handles_strings_with_missing_values(tmp_path: Path):
    """Regression for comment 5260723528 — step 4's ``atomic_write_h5ad``
    must succeed on obs with a plain ArrowString column that contains
    missing values AND is NOT pre-categorized, mirroring what step 4 sees
    after ``ad.read_h5ad(...)`` on step 1 / step 3 outputs."""
    _arrow_string_env_or_skip()
    anndata = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    import pandas as pd

    from rctd_split._internal.layout import atomic_write_h5ad

    n = 8
    fov_values = ["G7", "L7", None, "K8", None, "O4", "P7", "G7"]
    obs = pd.DataFrame({"fov": pd.array(fov_values, dtype="string")})
    ad = anndata.AnnData(
        X=np.zeros((n, 1), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["g0"]),
    )
    assert isinstance(ad.obs["fov"].dtype, pd.StringDtype)
    assert ad.obs["fov"].isna().sum() == 2

    out = tmp_path / "strings_with_na.h5ad"
    atomic_write_h5ad(ad, out)
    assert out.exists()

    back = anndata.read_h5ad(out)
    assert back.obs["fov"].isna().sum() == 2
    got = back.obs["fov"].astype(object).where(back.obs["fov"].notna(), None).tolist()
    assert got == fov_values

    # Caller's AnnData: fov dtype restored (not permanently categorized).
    assert isinstance(ad.obs["fov"].dtype, pd.StringDtype)


def test_atomic_write_h5ad_handles_arrow_string_index(tmp_path: Path):
    """Regression for comment 5260631648 (step 3 crash) — arrow-string obs
    index, the exact shape ``ad.read_h5ad(...)`` produces from a step 1 /
    step 3 output under pandas 3.x."""
    _arrow_string_env_or_skip()
    anndata = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    import pandas as pd

    from rctd_split._internal.layout import atomic_write_h5ad

    n = 5
    obs = pd.DataFrame(
        {"x": np.arange(n, dtype=np.int64)},
        index=pd.Index(
            pd.array([f"cell_{i}" for i in range(n)], dtype="string"),
            name="cell_id",
        ),
    )
    assert isinstance(obs.index.dtype, pd.StringDtype)

    ad = anndata.AnnData(
        X=np.zeros((n, 1), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=["g0"]),
    )

    out = tmp_path / "arrow_index.h5ad"
    atomic_write_h5ad(ad, out)
    assert out.exists()

    back = anndata.read_h5ad(out)
    assert back.obs.index.tolist() == [f"cell_{i}" for i in range(n)]
