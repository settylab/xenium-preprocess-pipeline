"""Unit tests for the qc_report stage (Stage 9, terminal).

Fabricates the three source h5ads with the minimum shape the stage
depends on:

  * ``proseg_purified.h5ad`` — carries ``obsm['X_umap_clipped_norm']``,
    ``obs['leiden_0.5']``, ``obs['first_type']``,
    ``obs['purification_status']`` (SPLIT-native categorical) and
    ``obs['nCount_Proseg']`` (Seurat-managed purified count column).
  * ``proseg_raw.h5ad`` — superset of purified cells; carries
    ``obs['total_counts']`` and ``layers['maxpost_counts']`` (the
    argmax-posterior counts matrix Tracy tracks).
  * ``xenium_ranger.h5ad`` — arbitrary counts + ``obs['total_counts']``
    + ``obsm['spatial']``.

Skipped when matplotlib is not importable — the stage renders PNGs.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import pytest


pytest.importorskip("anndata")
pytest.importorskip("scipy")
pytest.importorskip("pandas")
pytest.importorskip("matplotlib")


RUN_ID = "42"
SAMPLE = "MHTEST"

# purified cells: obs_names 0..7. raw is a superset (0..11) so the
# purification funnel is non-trivial. xenium is disjoint (Q-names).
_N_PURIFIED = 8
_N_RAW = 12
_N_XENIUM = 6
_N_GENES = 5


def _write_purified(
    output_root: Path,
    *,
    with_purification_status: bool = True,
    with_umap: bool = True,
    with_leiden: bool = True,
    with_first_type: bool = True,
    with_ncount_proseg: bool = True,
    extra_obs_cols: dict | None = None,
):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    from rctd_split._internal.layout import spatial_adata_path

    rng = np.random.default_rng(0)
    counts = rng.poisson(2.0, size=(_N_PURIFIED, _N_GENES)).astype("float32")
    X = csr_matrix(counts)

    obs = pd.DataFrame(index=pd.Index(
        [f"C{i:03d}" for i in range(_N_PURIFIED)], name="cell_id",
    ))
    if with_first_type:
        obs["first_type"] = (["tumor"] * 4) + (["stroma"] * 4)
    if with_leiden:
        # Deterministic leiden labels, 3 clusters.
        obs["leiden_0.5"] = pd.Categorical(
            (["0"] * 3) + (["1"] * 3) + (["2"] * 2),
        )
    if with_purification_status:
        obs["purification_status"] = pd.Categorical(
            (["purified"] * 5) + (["flagged"] * 2) + (["singlet"] * 1),
        )
    if with_ncount_proseg:
        # Non-zero so log10 histogram has a real distribution.
        obs["nCount_Proseg"] = counts.sum(axis=1).astype(float) + 1.0
    if extra_obs_cols:
        for k, v in extra_obs_cols.items():
            obs[k] = v
    var = pd.DataFrame(index=[f"gene{i}" for i in range(_N_GENES)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    if with_umap:
        # Trivial (deterministic) UMAP layout: 2D grid.
        xs = np.linspace(-1.0, 1.0, _N_PURIFIED, dtype=float)
        ys = np.linspace(-1.0, 1.0, _N_PURIFIED, dtype=float)[::-1]
        adata.obsm["X_umap_clipped_norm"] = np.column_stack([xs, ys])

    out = spatial_adata_path(output_root, SAMPLE, RUN_ID, "proseg_purified")
    out.parent.mkdir(parents=True, exist_ok=True)
    from rctd_split._internal.layout import atomic_write_h5ad
    atomic_write_h5ad(adata, out, compression="gzip")
    return out


def _write_raw(
    output_root: Path,
    *,
    with_maxpost_layer: bool = True,
    with_total_counts_obs: bool = True,
    with_spot_class: bool = True,
    with_first_type: bool = True,
    include_zero_count_cell: bool = False,
):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    from rctd_split._internal.layout import spatial_adata_path

    rng = np.random.default_rng(1)
    x_arr = rng.poisson(2.0, size=(_N_RAW, _N_GENES)).astype("float32")
    X = csr_matrix(x_arr)

    # Give maxpost a distinctly different distribution than X so tests
    # can prove the stage reads the layer, not .X.
    maxpost_arr = rng.poisson(5.0, size=(_N_RAW, _N_GENES)).astype("float32")

    cell_ids = [f"C{i:03d}" for i in range(_N_RAW)]
    obs = pd.DataFrame(index=pd.Index(cell_ids, name="cell_id"))
    total_counts_col = maxpost_arr.sum(axis=1).astype(float) + 1.0
    if include_zero_count_cell:
        # Force the last raw cell to zero so the histogram
        # "did not pass the threshold or generally low count" branch
        # is exercised.
        total_counts_col[-1] = 0.0
    if with_total_counts_obs:
        obs["total_counts"] = total_counts_col
    if with_spot_class:
        # Mix of canonical spacexr categories so both the summary
        # tabulation and the rejected-first_type breakdown are exercised.
        # Cells 0..2 → singlet
        # Cells 3..5 → doublet_certain
        # Cells 6..7 → doublet_uncertain
        # Cells 8..9 → reject
        # Cells 10..11 → "" (step-1 QC dropped them BEFORE RCTD;
        # empty spot_class is the sentinel).
        spot = (
            (["singlet"] * 3)
            + (["doublet_certain"] * 3)
            + (["doublet_uncertain"] * 2)
            + (["reject"] * 2)
            + ([""] * 2)
        )
        assert len(spot) == _N_RAW
        obs["spot_class"] = spot
    if with_first_type:
        # Assign first_type per cell — rejected cells get a mix of
        # tumor/stroma so the rejected-first_type table has ≥1 row per.
        first_type = (
            (["tumor"] * 5)
            + (["stroma"] * 3)
            # rejected cells (indices 8, 9): one tumor, one stroma.
            + ["tumor", "stroma"]
            + ([""] * 2)
        )
        assert len(first_type) == _N_RAW
        obs["first_type"] = first_type
    var = pd.DataFrame(index=[f"gene{i}" for i in range(_N_GENES)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    if with_maxpost_layer:
        adata.layers["maxpost_counts"] = csr_matrix(maxpost_arr)

    out = spatial_adata_path(output_root, SAMPLE, RUN_ID, "proseg_raw")
    out.parent.mkdir(parents=True, exist_ok=True)
    from rctd_split._internal.layout import atomic_write_h5ad
    atomic_write_h5ad(adata, out, compression="gzip")
    return out


def _write_xenium(output_root: Path, *, with_total_counts_obs: bool = True):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    from rctd_split._internal.layout import spatial_adata_path

    rng = np.random.default_rng(2)
    counts = rng.poisson(2.0, size=(_N_XENIUM, _N_GENES)).astype("float32")
    X = csr_matrix(counts)
    obs = pd.DataFrame(index=pd.Index(
        [f"Q{i}" for i in range(_N_XENIUM)], name="cell_id",
    ))
    if with_total_counts_obs:
        obs["total_counts"] = counts.sum(axis=1).astype(float) + 1.0
    var = pd.DataFrame(index=[f"gene{i}" for i in range(_N_GENES)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm["spatial"] = rng.uniform(0, 100, size=(_N_XENIUM, 2))

    out = spatial_adata_path(output_root, SAMPLE, RUN_ID, "xenium_ranger")
    out.parent.mkdir(parents=True, exist_ok=True)
    from rctd_split._internal.layout import atomic_write_h5ad
    atomic_write_h5ad(adata, out, compression="gzip")
    return out


def _write_resolved_config(
    output_root: Path,
    *,
    step1_min_counts_cell=10,
    step4_postprocess_min_counts=50,
):
    """Write a minimal ``config.yaml`` under the run dir so
    ``_read_hist_thresholds`` picks up the two threshold values.

    Pass ``None`` for either threshold to omit that config key (tests
    the "no threshold line" fallback path).
    """
    import yaml

    from rctd_split._internal.layout import resolved_config_path

    p = resolved_config_path(output_root, SAMPLE, RUN_ID)
    p.parent.mkdir(parents=True, exist_ok=True)

    merged: dict = {}
    if step1_min_counts_cell is not None:
        merged["step1"] = {"qc_filter": {
            "min_counts_cell": step1_min_counts_cell,
        }}
    if step4_postprocess_min_counts is not None:
        merged["step4"] = {"postprocess": {"qc": {
            "min_counts": step4_postprocess_min_counts,
        }}}
    with open(p, "w") as f:
        yaml.safe_dump(merged, f, sort_keys=False)
    return p


def _setup_run(tmp_path: Path, *, purified_kwargs=None, raw_kwargs=None,
               xenium_kwargs=None, resolved_config_kwargs=None):
    output_root = tmp_path / "runs"
    purified = _write_purified(output_root, **(purified_kwargs or {}))
    raw = _write_raw(output_root, **(raw_kwargs or {}))
    xen = _write_xenium(output_root, **(xenium_kwargs or {}))
    if resolved_config_kwargs is not None:
        _write_resolved_config(output_root, **resolved_config_kwargs)
    else:
        # Default: write the canonical demo thresholds so the histogram
        # threshold-line assertions have data to grab.
        _write_resolved_config(output_root)
    return output_root, purified, raw, xen


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


_DEFAULT_CALL = dict(
    purification_status_column="purification_status",
    raw_layer="maxpost_counts",
)


def test_qc_report_happy_path(tmp_path: Path):
    """Full happy path: HTML + 3 UMAP PNGs + 3 histogram PNGs + metrics
    CSV + sentinel emitted."""
    from rctd_split._internal.layout import (
        intermediate_path,
        summary_path,
        summary_plots_dir,
    )
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(tmp_path)

    sentinel = run_qc_report(
        sample_id=SAMPLE,
        run_id=RUN_ID,
        output_root=output_root,
        force_rerun=False,
        **_DEFAULT_CALL,
    )
    assert sentinel.exists()
    assert sentinel == intermediate_path(
        output_root, SAMPLE, RUN_ID, "qc_report_sentinel",
    )

    html_out = summary_path(output_root, SAMPLE, RUN_ID, "html_report")
    metrics_out = summary_path(output_root, SAMPLE, RUN_ID, "metrics_csv")
    assert html_out.exists()
    assert metrics_out.exists()

    plots = summary_plots_dir(output_root, SAMPLE, RUN_ID)
    assert (plots / f"{SAMPLE}_umap_purification_status.png").exists()
    assert (plots / f"{SAMPLE}_umap_leiden_res0.5.png").exists()
    assert (plots / f"{SAMPLE}_umap_first_type.png").exists()
    assert (plots / f"{SAMPLE}_hist_total_counts_proseg_raw.png").exists()
    assert (plots / f"{SAMPLE}_hist_total_counts_xenium_ranger.png").exists()
    assert (
        plots / f"{SAMPLE}_hist_nCount_Proseg_proseg_purified.png"
    ).exists()

    body = html_out.read_text()
    assert "Per-h5ad metrics" in body
    assert "Purification funnel" in body
    assert "Count distributions" in body
    assert "purification_status" in body
    assert "nCount_Proseg" in body
    assert "resolution = 0.5" in body
    assert "first_type" in body
    # The user-facing UMAP + rejected-cells section headings use the
    # display label "SPLIT-inferred celltype" (Tracy's ask on
    # #26, comment 5320970944). The underlying data column is still
    # `first_type`, which is why the substring above still appears in
    # the body via the descriptive paragraphs' <code>first_type</code>.
    assert "SPLIT-inferred celltype" in body
    # The raw metrics row records the layer used (maxpost_counts).
    assert "maxpost_counts" in body
    # RCTD summary section + rejected SPLIT-inferred celltype breakdown.
    assert "RCTD summary" in body
    assert "SPLIT-inferred celltype in rejected cells" in body
    assert "singlet" in body
    assert "doublet_certain" in body
    assert "doublet_uncertain" in body
    assert "reject" in body
    # Threshold captions include the config-source annotation.
    assert "step1.qc_filter.min_counts_cell" in body
    assert "step4.postprocess.qc.min_counts" in body

    # RCTD summary CSV is a separate sidecar.
    from rctd_split._internal.layout import summary_path
    rctd_csv = summary_path(output_root, SAMPLE, RUN_ID, "rctd_summary_csv")
    assert rctd_csv.exists()

    # Color map JSON sidecar is written alongside the HTML report;
    # both categorical keys are present so downstream consumers can
    # reproduce the fixed per-level coloring.
    import json
    from rctd_split._internal.layout import summary_dir
    color_map_path = summary_dir(output_root, SAMPLE, RUN_ID) / (
        f"{SAMPLE}_color_map.json"
    )
    assert color_map_path.exists()
    color_map = json.loads(color_map_path.read_text())
    assert set(color_map.keys()) == {"purification_status", "first_type"}
    # The HTML captions link to the JSON sidecar by basename.
    assert color_map_path.name in body


def test_qc_report_html_provenance_section(tmp_path: Path):
    """Provenance section in the HTML captures the invoking command line
    and embeds the full merged config content (settylab/TracyY123-nexus#26
    comment 5320970944, item 5)."""
    from rctd_split._internal.layout import (
        resolved_config_path,
        summary_path,
    )
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(tmp_path)
    argv = [
        "rctd-split", "run",
        "--sample-id", SAMPLE, "--run-id", RUN_ID,
        "--output-root", str(output_root),
    ]

    run_qc_report(
        sample_id=SAMPLE, run_id=RUN_ID,
        output_root=output_root, force_rerun=False,
        invoking_argv=argv,
        **_DEFAULT_CALL,
    )

    body = summary_path(
        output_root, SAMPLE, RUN_ID, "html_report",
    ).read_text()

    # Provenance heading + invocation row.
    assert "<h2>Provenance</h2>" in body
    assert "Invocation" in body
    assert " ".join(argv) in body

    # Resolved config path + full content embedded via <details>/<pre>.
    cfg_path = str(resolved_config_path(output_root, SAMPLE, RUN_ID))
    assert cfg_path in body
    cfg_text = resolved_config_path(
        output_root, SAMPLE, RUN_ID,
    ).read_text()
    # The fixture writes step1.qc_filter.min_counts_cell — pick a
    # substring that would only appear if the config text made it into
    # the report body.
    assert "min_counts_cell" in cfg_text
    assert "min_counts_cell" in body


def test_qc_report_leaves_source_h5ads_untouched(tmp_path: Path):
    """Reproducibility invariant: qc_report is READ-ONLY on the three
    source h5ads. Hash before/after to prove no accidental rewrites."""
    from rctd_split.stages.qc_report import run_qc_report

    output_root, purified, raw, xen = _setup_run(tmp_path)
    hashes_before = {p.name: _sha256(p) for p in (purified, raw, xen)}

    run_qc_report(
        sample_id=SAMPLE,
        run_id=RUN_ID,
        output_root=output_root,
        force_rerun=False,
        **_DEFAULT_CALL,
    )

    hashes_after = {p.name: _sha256(p) for p in (purified, raw, xen)}
    assert hashes_before == hashes_after


def test_qc_report_raw_metrics_use_maxpost_not_X(tmp_path: Path):
    """The `proseg_raw` metrics row must be computed off
    ``layers['maxpost_counts']`` (not ``.X``). Since the fixtures give
    the two matrices distinctly different distributions, the median-
    counts number lands on the layer's value, not X's."""
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    from anndata import read_h5ad

    from rctd_split._internal.layout import summary_path
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, raw_p, _ = _setup_run(tmp_path)

    run_qc_report(
        sample_id=SAMPLE, run_id=RUN_ID,
        output_root=output_root, force_rerun=False, **_DEFAULT_CALL,
    )

    metrics = pd.read_csv(summary_path(output_root, SAMPLE, RUN_ID, "metrics_csv"))
    raw_row = metrics[metrics["h5ad"] == "proseg_raw"].iloc[0]
    assert raw_row["layer"] == "maxpost_counts"

    raw = read_h5ad(raw_p)
    layer = raw.layers["maxpost_counts"]
    expected_median = float(np.median(
        np.asarray(layer.sum(axis=1)).ravel()
        if sp.issparse(layer) else np.asarray(layer).sum(axis=1)
    ))
    x_median = float(np.median(
        np.asarray(raw.X.sum(axis=1)).ravel()
        if sp.issparse(raw.X) else np.asarray(raw.X).sum(axis=1)
    ))
    # Guard: fixtures pick distinct Poisson means (2 vs. 5), so the
    # medians must actually differ. If someone tunes the fixture and
    # this asserts, the invariant test needs new numeric separation.
    assert expected_median != x_median
    assert raw_row["median_counts_per_cell"] == pytest.approx(expected_median)


def test_qc_report_fails_on_missing_maxpost_layer(tmp_path: Path):
    """Fail-loud when proseg_raw has no ``layers['maxpost_counts']``
    — the error message must name the layer + the source h5ad."""
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(
        tmp_path, raw_kwargs={"with_maxpost_layer": False},
    )
    with pytest.raises(SystemExit) as exc:
        run_qc_report(
            sample_id=SAMPLE, run_id=RUN_ID,
            output_root=output_root, force_rerun=False, **_DEFAULT_CALL,
        )
    msg = str(exc.value)
    assert "maxpost" in msg
    assert "proseg_raw" in msg


def test_qc_report_fails_on_missing_purification_status(tmp_path: Path):
    """Fail-loud when proseg_purified.obs has no ``purification_status``
    — the message must name the column + point at SPLIT::purify."""
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(
        tmp_path, purified_kwargs={"with_purification_status": False},
    )
    with pytest.raises(SystemExit) as exc:
        run_qc_report(
            sample_id=SAMPLE, run_id=RUN_ID,
            output_root=output_root, force_rerun=False, **_DEFAULT_CALL,
        )
    msg = str(exc.value)
    assert "purification_status" in msg
    assert "SPLIT" in msg


def test_qc_report_fails_on_missing_umap(tmp_path: Path):
    """Fail-loud when postprocess hasn't laid down a UMAP embedding."""
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(
        tmp_path, purified_kwargs={"with_umap": False},
    )

    with pytest.raises(SystemExit) as exc:
        run_qc_report(
            sample_id=SAMPLE, run_id=RUN_ID,
            output_root=output_root, force_rerun=False, **_DEFAULT_CALL,
        )
    msg = str(exc.value)
    assert "UMAP" in msg
    assert "postprocess" in msg


def test_qc_report_fails_on_missing_hist_obs(tmp_path: Path):
    """Fail-loud when purified has no ``nCount_Proseg`` obs col."""
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(
        tmp_path, purified_kwargs={"with_ncount_proseg": False},
    )
    with pytest.raises(SystemExit) as exc:
        run_qc_report(
            sample_id=SAMPLE, run_id=RUN_ID,
            output_root=output_root, force_rerun=False, **_DEFAULT_CALL,
        )
    msg = str(exc.value)
    assert "nCount_Proseg" in msg
    assert "proseg_purified" in msg


def test_qc_report_sentinel_short_circuits(tmp_path: Path):
    """Second invocation with sentinel present must skip work — the
    plots-dir mtime shouldn't change."""
    from rctd_split._internal.layout import summary_plots_dir
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(tmp_path)
    common = dict(
        sample_id=SAMPLE, run_id=RUN_ID, output_root=output_root,
        **_DEFAULT_CALL,
    )
    run_qc_report(force_rerun=False, **common)
    plot_path = (
        summary_plots_dir(output_root, SAMPLE, RUN_ID)
        / f"{SAMPLE}_umap_first_type.png"
    )
    mtime1 = plot_path.stat().st_mtime_ns

    # Delete the plot then re-run w/o force_rerun; sentinel skip means
    # it stays deleted (proving no work was done).
    plot_path.unlink()
    run_qc_report(force_rerun=False, **common)
    assert not plot_path.exists()

    # force_rerun=True: plot regenerated.
    run_qc_report(force_rerun=True, **common)
    assert plot_path.exists()
    assert plot_path.stat().st_mtime_ns != mtime1


# ---------------------------------------------------------------------------
# RCTD summary
# ---------------------------------------------------------------------------


def test_qc_report_rctd_summary_counts(tmp_path: Path):
    """RCTD spot_class counts + rejected first_type breakdown must match
    the fixture-planted values on ``raw.obs``.

    Fixture (see ``_write_raw``): 3 singlet, 3 doublet_certain,
    2 doublet_uncertain, 2 reject, 2 empty (step-1 dropped).
    Rejected cells' first_type = {tumor, stroma} (one each).
    """
    import pandas as pd

    from rctd_split._internal.layout import summary_path
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(tmp_path)

    run_qc_report(
        sample_id=SAMPLE, run_id=RUN_ID,
        output_root=output_root, force_rerun=False, **_DEFAULT_CALL,
    )

    csv_path = summary_path(output_root, SAMPLE, RUN_ID, "rctd_summary_csv")
    assert csv_path.exists(), csv_path
    df = pd.read_csv(csv_path)

    # Summary section: denominators.
    summary = df[df["section"] == "summary"].set_index("key")
    assert int(summary.loc["n_raw", "count"]) == _N_RAW  # 12
    # 2 empty spot_class cells → step-1-dropped.
    assert int(summary.loc["n_pre_rctd_dropped", "count"]) == 2
    assert int(summary.loc["n_rctd", "count"]) == _N_RAW - 2  # 10
    assert int(summary.loc["n_rejected", "count"]) == 2

    # Spot_class section: exact counts + percentages of RCTD denominator (10).
    sc = df[df["section"] == "spot_class"].set_index("key")
    assert int(sc.loc["singlet", "count"]) == 3
    assert sc.loc["singlet", "pct"] == pytest.approx(30.0)
    assert int(sc.loc["doublet_certain", "count"]) == 3
    assert sc.loc["doublet_certain", "pct"] == pytest.approx(30.0)
    assert int(sc.loc["doublet_uncertain", "count"]) == 2
    assert sc.loc["doublet_uncertain", "pct"] == pytest.approx(20.0)
    assert int(sc.loc["reject", "count"]) == 2
    assert sc.loc["reject", "pct"] == pytest.approx(20.0)

    # Rejected-first_type section: 1 tumor, 1 stroma = 50% each.
    rft = df[df["section"] == "rejected_first_type"].set_index("key")
    assert int(rft.loc["tumor", "count"]) == 1
    assert rft.loc["tumor", "pct"] == pytest.approx(50.0)
    assert int(rft.loc["stroma", "count"]) == 1
    assert rft.loc["stroma", "pct"] == pytest.approx(50.0)


def test_qc_report_rctd_html_shows_percentages(tmp_path: Path):
    """The HTML report must render the RCTD summary rows with counts +
    percentages that agree with the CSV."""
    from rctd_split._internal.layout import summary_path
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(tmp_path)

    run_qc_report(
        sample_id=SAMPLE, run_id=RUN_ID,
        output_root=output_root, force_rerun=False, **_DEFAULT_CALL,
    )

    body = summary_path(output_root, SAMPLE, RUN_ID, "html_report").read_text()
    # Percent strings from the fixture (RCTD denominator = 10).
    assert "30.00%" in body   # singlet + doublet_certain
    assert "20.00%" in body   # doublet_uncertain + reject
    assert "50.00%" in body   # tumor + stroma in rejected


def test_qc_report_fails_on_missing_spot_class(tmp_path: Path):
    """Fail-loud when raw.obs has no spot_class — the RCTD summary
    cannot be computed without it."""
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(
        tmp_path, raw_kwargs={"with_spot_class": False},
    )
    with pytest.raises(SystemExit) as exc:
        run_qc_report(
            sample_id=SAMPLE, run_id=RUN_ID,
            output_root=output_root, force_rerun=False, **_DEFAULT_CALL,
        )
    msg = str(exc.value)
    assert "spot_class" in msg
    assert "proseg_raw" in msg
    assert "writeback_to_step1_raw" in msg


def test_qc_report_fails_on_missing_first_type_on_raw(tmp_path: Path):
    """Fail-loud when raw.obs has no first_type — the rejected-cells
    breakdown cannot be computed without it."""
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(
        tmp_path, raw_kwargs={"with_first_type": False},
    )
    with pytest.raises(SystemExit) as exc:
        run_qc_report(
            sample_id=SAMPLE, run_id=RUN_ID,
            output_root=output_root, force_rerun=False, **_DEFAULT_CALL,
        )
    msg = str(exc.value)
    assert "first_type" in msg
    assert "proseg_raw" in msg


# ---------------------------------------------------------------------------
# Histogram: all cells + threshold line
# ---------------------------------------------------------------------------


def test_qc_report_hist_threshold_from_resolved_config(tmp_path: Path):
    """Threshold values are read from ``config.yaml`` — the
    HTML surface names the source config keys AND the numeric values.
    """
    from rctd_split._internal.layout import summary_path
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(
        tmp_path,
        resolved_config_kwargs={
            "step1_min_counts_cell": 7,
            "step4_postprocess_min_counts": 42,
        },
    )
    run_qc_report(
        sample_id=SAMPLE, run_id=RUN_ID,
        output_root=output_root, force_rerun=False, **_DEFAULT_CALL,
    )

    body = summary_path(output_root, SAMPLE, RUN_ID, "html_report").read_text()
    # Both the config-path anchor and the numeric threshold should render.
    assert "step1.qc_filter.min_counts_cell" in body
    assert "<b>7</b>" in body
    assert "step4.postprocess.qc.min_counts" in body
    assert "<b>42</b>" in body


def test_qc_report_hist_no_threshold_when_config_absent(tmp_path: Path):
    """When ``config.yaml`` is missing entirely, the histograms
    still render and the HTML explicitly says "no dashed threshold
    line drawn" — never falls back to a hard-coded value."""
    from rctd_split._internal.layout import (
        summary_path,
        resolved_config_path,
    )
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(tmp_path)
    # Nuke the resolved_config the fixture defaulted-in above.
    resolved_config_path(output_root, SAMPLE, RUN_ID).unlink()

    run_qc_report(
        sample_id=SAMPLE, run_id=RUN_ID,
        output_root=output_root, force_rerun=False, **_DEFAULT_CALL,
    )
    body = summary_path(output_root, SAMPLE, RUN_ID, "html_report").read_text()
    # Caption for every histogram falls back to the no-line note.
    assert body.count("no dashed threshold line drawn") == 3


def test_qc_report_hist_xenium_never_has_threshold(tmp_path: Path):
    """xenium_ranger.h5ad has no upstream min-counts gate in step-1;
    even with a fully-populated resolved_config the xenium histogram
    caption stays on the "no line" branch."""
    from rctd_split._internal.layout import summary_path
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(tmp_path)
    run_qc_report(
        sample_id=SAMPLE, run_id=RUN_ID,
        output_root=output_root, force_rerun=False, **_DEFAULT_CALL,
    )
    body = summary_path(output_root, SAMPLE, RUN_ID, "html_report").read_text()
    # Locate the xenium histogram block by its <img> alt text and pull
    # a chunk following it — the caption we assert against lives there.
    xenium_anchor = 'alt="log10(xenium_ranger total_counts)"'
    xenium_idx = body.index(xenium_anchor)
    xenium_block = body[xenium_idx: xenium_idx + 400]
    assert "no dashed threshold line drawn" in xenium_block


def test_qc_report_hist_handles_zero_count_cells(tmp_path: Path):
    """A cell with total_counts=0 must not crash the histogram code
    path (the "did not pass the threshold or generally low count"
    branch is exercised), and the resulting PNG is written non-empty."""
    from rctd_split._internal.layout import summary_plots_dir
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(
        tmp_path, raw_kwargs={"include_zero_count_cell": True},
    )
    run_qc_report(
        sample_id=SAMPLE, run_id=RUN_ID,
        output_root=output_root, force_rerun=False, **_DEFAULT_CALL,
    )
    hist_png = (
        summary_plots_dir(output_root, SAMPLE, RUN_ID)
        / f"{SAMPLE}_hist_total_counts_proseg_raw.png"
    )
    assert hist_png.exists()
    assert hist_png.stat().st_size > 0


def test_qc_report_read_hist_thresholds_unit(tmp_path: Path):
    """Direct unit-test of ``_read_hist_thresholds`` — the private
    helper that maps YAML keys onto the three histogram threshold
    slots. Missing file, missing sub-keys, and non-numeric values all
    resolve to None (never a fallback constant)."""
    import yaml

    from rctd_split.stages.qc_report import _read_hist_thresholds

    # 1. Missing file → all None.
    absent = tmp_path / "does_not_exist.yaml"
    got = _read_hist_thresholds(absent)
    assert got == {"raw": None, "xenium": None, "purified": None}

    # 2. Populated file → both numeric thresholds resolved as float;
    #    xenium always None (no config key defined).
    good = tmp_path / "good.yaml"
    good.write_text(yaml.safe_dump({
        "step1": {"qc_filter": {"min_counts_cell": 11}},
        "step4": {"postprocess": {"qc": {"min_counts": 33}}},
    }))
    got = _read_hist_thresholds(good)
    assert got["raw"] == 11.0
    assert got["xenium"] is None
    assert got["purified"] == 33.0

    # 3. Non-numeric config value → None (defensive).
    bad = tmp_path / "bad.yaml"
    bad.write_text(yaml.safe_dump({
        "step1": {"qc_filter": {"min_counts_cell": "not a number"}},
    }))
    got = _read_hist_thresholds(bad)
    assert got["raw"] is None
    assert got["purified"] is None


# ---------------------------------------------------------------------------
# extra_reports (Tracy comment 5322126093, item 4)
# ---------------------------------------------------------------------------


def test_qc_report_extra_reports_absent_no_section(tmp_path: Path):
    """No extra_reports = no "Additional reports" heading in the HTML."""
    from rctd_split._internal.layout import summary_path
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(tmp_path)
    run_qc_report(
        sample_id=SAMPLE, run_id=RUN_ID,
        output_root=output_root, force_rerun=False,
        **_DEFAULT_CALL,
    )
    body = summary_path(output_root, SAMPLE, RUN_ID, "html_report").read_text()
    assert "Additional reports" not in body


def test_qc_report_extra_reports_renders_section(tmp_path: Path):
    """When present, the section renders one heading + iframe + link
    per entry, and the source path appears verbatim in the caption."""
    from rctd_split._internal.layout import summary_path
    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(tmp_path)

    # Put one external report next to the summary/ output; it should
    # be resolved to a RELATIVE path in the href (portable summary/).
    external = tmp_path / "external_diagnostics.html"
    external.write_text("<html><body>hello</body></html>")

    run_qc_report(
        sample_id=SAMPLE, run_id=RUN_ID,
        output_root=output_root, force_rerun=False,
        extra_reports=[
            {"path": str(external),
             "name": "MH-alt diagnostics, run #7"},
            {"path": "/does/not/exist.html",
             "name": "missing on purpose"},
        ],
        **_DEFAULT_CALL,
    )
    body = summary_path(output_root, SAMPLE, RUN_ID, "html_report").read_text()

    assert "Additional reports" in body
    # Display name (with an embedded comma) survives round-trip.
    assert "MH-alt diagnostics, run #7" in body
    assert "missing on purpose" in body
    # Missing file gets an inline "missing at render time" note.
    assert "missing at render time" in body
    # Source paths are surfaced verbatim in the meta line.
    assert "external_diagnostics.html" in body
    assert "/does/not/exist.html" in body
    # An iframe embeds each report.
    assert body.count("<iframe") == 2


def test_qc_report_extra_reports_rejects_malformed(tmp_path: Path):
    """A malformed entry (missing 'path' or 'name') is rejected
    fail-loud before the HTML is written."""
    import pytest

    from rctd_split.stages.qc_report import run_qc_report

    output_root, _, _, _ = _setup_run(tmp_path)
    with pytest.raises(SystemExit):
        run_qc_report(
            sample_id=SAMPLE, run_id=RUN_ID,
            output_root=output_root, force_rerun=False,
            extra_reports=[{"path": "/tmp/foo.html"}],  # missing name
            **_DEFAULT_CALL,
        )
