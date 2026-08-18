"""Stage 9 (terminal): QC HTML report + UMAP plots + count histograms.

Reads (never writes to) the three persisted h5ads under
``spatial_adata/``:

  * ``<S>_xenium_ranger.h5ad`` (step 1; augmented by celltype_writeback)
  * ``<S>_proseg_raw.h5ad``    (step 1; augmented by writeback_to_step1_raw)
  * ``<S>_proseg_purified.h5ad`` (this pipeline; augmented by postprocess)

Emits under ``<run_dir>/summary/``:

  * ``<S>_qc_report.html`` — QC metrics tables + embedded plots + versions.
  * ``<S>_qc_metrics.csv`` — machine-parseable version of the metrics.
  * ``<S>_rctd_summary.csv`` — RCTD spot_class + rejected-first_type table.
  * ``plots/<S>_umap_purification_status.png``
  * ``plots/<S>_umap_leiden_res<r>.png``
  * ``plots/<S>_umap_first_type.png``
  * ``plots/<S>_hist_total_counts_proseg_raw.png``
  * ``plots/<S>_hist_total_counts_xenium_ranger.png``
  * ``plots/<S>_hist_nCount_Proseg_proseg_purified.png``

Reproducibility:

  * All matplotlib figures are drawn from deterministic array orderings
    (obs_names is the canonical order).
  * UMAP #1 is colored by ``proseg_purified.obs['purification_status']``
    (SPLIT-native categorical). Fail-loud when the column is absent
    (settylab/TracyY123-nexus#26 comment 5274781435).
  * The stage does NOT recompute UMAP — it uses the ``X_umap_clipped_norm``
    obsm laid down by ``postprocess`` (fail-loud if missing).
  * Leiden clustering is likewise consumed from ``postprocess`` output —
    the resolution used for the leiden UMAP plot is picked from the
    lowest-numbered ``leiden_*`` obs column and recorded in the plot
    filename + HTML legend.
  * Raw-side QC metrics + histogram are computed on
    ``raw.layers['maxpost_counts']`` (proseg emits multiple assays;
    ``maxpost_counts`` is the argmax-posterior integer matrix Tracy
    tracks — settylab/TracyY123-nexus#26 comment 5274767383;
    canonical layer name confirmed in comment 5274956572). Fail-loud
    when the layer is absent.
  * The three count histograms plot the distribution for ALL cells in
    the source h5ad (no positive-only slice) and overlay a dashed
    vertical line at the upstream filter threshold, read from the
    merged ``config.yaml`` (settylab/TracyY123-nexus#26
    comment 5275064492). Thresholds:
      - proseg_raw histogram: ``step1.qc_filter.min_counts_cell``
      - xenium_ranger histogram: no upstream min-counts gate, no line
      - proseg_purified histogram: ``step4.postprocess.qc.min_counts``
    If the config file or a threshold key is absent, the line is
    omitted — never hard-coded.
  * RCTD spot_class summary + first_type-in-rejected tabulation are
    computed from ``proseg_raw.obs`` (``spot_class`` + ``first_type``
    are folded from ``intermediate/adata/step4_unpurified.h5ad`` onto
    raw by ``writeback_to_step1_raw`` Part 2; ``unpurified.obs`` in
    turn got them from RCTD's ``results_df`` via
    ``SPLIT::run_post_process_RCTD``). Fail-loud when either column is
    absent. Cells with empty ``spot_class`` were dropped by step-1
    QC before RCTD ever saw them and are excluded from the RCTD
    denominator.

Sentinel: marker file under ``intermediate/adata/`` (per pipeline
convention). Idempotent — re-runs are guarded by the sentinel unless
``force_rerun`` is set.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

from rctd_split._internal.compat import sentinel_exists
from rctd_split._internal.layout import (
    intermediate_path,
    summary_path,
    summary_plots_dir,
    resolved_config_path,
    spatial_adata_path,
)
from rctd_split._internal.logging import log


# Canonical RCTD spot_class values (spacexr, post
# SPLIT::run_post_process_RCTD). Order used in the HTML summary table +
# CSV. Any spot_class value observed on raw.obs that isn't in this
# tuple is still tabulated — appended in sorted order after the
# canonical rows so unexpected labels don't get silently dropped.
_SPOT_CLASS_CANONICAL_ORDER: tuple[str, ...] = (
    "singlet",
    "doublet_certain",
    "doublet_uncertain",
    "reject",
)
# Display labels for the HTML table. Tracy's ask (comment 5275064492)
# names four categories: singlet / doublet / rejected / uncertain — map
# them onto RCTD's actual values here.
_SPOT_CLASS_DISPLAY: dict[str, str] = {
    "singlet": "Singlet",
    "doublet_certain": "Doublet (certain)",
    "doublet_uncertain": "Uncertain (doublet_uncertain)",
    "reject": "Rejected",
}
# spot_class values that count as "rejected" for the first_type breakdown.
# Accept both `reject` (canonical spacexr) and `rejected` (defensive
# alias in case a downstream package renames).
_REJECT_VALUES: frozenset[str] = frozenset({"reject", "rejected"})


_PLOT_DPI = 150
_PLOT_FIGSIZE = (6.0, 5.5)
_PLOT_POINT_SIZE = 4.0
_PLOT_ALPHA = 0.7
_PALETTE = "tab20"


# ---------------------------------------------------------------------------
# Metric computation
# ---------------------------------------------------------------------------


def _resolve_matrix(adata, layer: str | None, label: str):
    """Return the expression matrix for `adata` — either ``adata.X`` when
    ``layer is None``, or ``adata.layers[layer]`` (fail-loud on absence).

    Central helper so the QC metrics and the histogram plot share the
    exact same numeric source per h5ad.
    """
    if layer is None:
        return adata.X
    if layer not in adata.layers:
        available = list(adata.layers.keys())
        raise SystemExit(
            f"[qc_report] {label}.h5ad has no layers[{layer!r}]. "
            f"Available layers: {available}. "
            "The canonical proseg argmax-posterior integer count "
            "layer is 'maxpost_counts' (as emitted by step-1 proseg "
            "export). If your h5ad uses that name, drop --qc-raw-layer "
            "or set qc_report.raw_layer: maxpost_counts in your config."
        )
    return adata.layers[layer]


def _h5ad_metrics(adata, label: str, layer: str | None = None) -> dict:
    """Compute per-h5ad QC metrics. `label` is a human-readable tag
    (e.g. 'xenium_ranger') that goes into the metrics row. `layer` is
    the layer name to use for the counts matrix (``None`` = ``.X``)."""
    import numpy as np
    import scipy.sparse as sp

    X = _resolve_matrix(adata, layer, label)
    n_cells, n_genes = adata.shape

    if sp.issparse(X):
        row_sums = np.asarray(X.sum(axis=1)).ravel()
        genes_per_cell = np.asarray((X > 0).sum(axis=1)).ravel()
    else:
        arr = np.asarray(X)
        row_sums = arr.sum(axis=1)
        genes_per_cell = (arr > 0).sum(axis=1)

    median_counts = float(np.median(row_sums)) if n_cells else 0.0
    median_genes = float(np.median(genes_per_cell)) if n_cells else 0.0
    mean_counts = float(np.mean(row_sums)) if n_cells else 0.0

    return {
        "h5ad": label,
        "layer": layer if layer is not None else "X",
        "n_cells": int(n_cells),
        "n_genes": int(n_genes),
        "median_counts_per_cell": median_counts,
        "mean_counts_per_cell": mean_counts,
        "median_genes_per_cell": median_genes,
    }


def _purification_metrics(raw_adata, purified_adata) -> dict:
    """Cell-count invariants tying raw → purified."""
    n_raw = int(raw_adata.n_obs)
    n_pur = int(purified_adata.n_obs)
    n_dropped = n_raw - n_pur
    pct_kept = (100.0 * n_pur / n_raw) if n_raw else 0.0
    return {
        "n_raw": n_raw,
        "n_purified": n_pur,
        "n_dropped_raw_to_purified": n_dropped,
        "pct_purified_of_raw": pct_kept,
    }


def _normalize_spot_series(series):
    """Coerce ``obs['spot_class']`` / ``obs['first_type']`` to a clean
    string Series. Categorical → str; NaN / empty / literal 'nan' →
    "" so downstream masks can filter with a single equality check.
    """
    import pandas as pd

    s = pd.Series(series).astype(object)
    s = s.where(~s.isna(), "")
    s = s.astype(str)
    lower = s.str.lower()
    s = s.where(~lower.isin(("nan", "none", "na")), "")
    return s


def _rctd_summary_metrics(raw_adata) -> dict:
    """Compute RCTD summary metrics from ``proseg_raw.obs``.

    Reads ``raw.obs['spot_class']`` and ``raw.obs['first_type']`` —
    both folded onto raw by ``writeback_to_step1_raw`` Part 2 from
    ``step4_unpurified.obs``. Cells with empty ``spot_class`` were
    dropped by step-1 QC before RCTD saw them and are excluded from
    the RCTD denominator (they are counted separately as
    ``n_pre_rctd_dropped``).

    Returns a dict with:

    - ``n_raw`` — total raw cells.
    - ``n_pre_rctd_dropped`` — cells with empty spot_class
      (step-1 qc-filtered, never reached RCTD).
    - ``n_rctd`` — cells RCTD categorized.
    - ``spot_class_rows`` — list of dicts with keys
      ``class``, ``display_name``, ``count``, ``pct``. Ordered by
      the canonical spacexr order then any unexpected extras.
    - ``n_rejected`` — cells with spot_class in ``_REJECT_VALUES``.
    - ``rejected_first_type_rows`` — list of dicts with keys
      ``first_type``, ``count``, ``pct``. Percentages are of the
      rejected subset. Empty first_type → ``"(none)"``.

    Fail-loud when either source column is absent.
    """
    for col in ("spot_class", "first_type"):
        if col not in raw_adata.obs.columns:
            raise SystemExit(
                f"[qc_report] proseg_raw.h5ad missing obs[{col!r}] — "
                "expected after writeback_to_step1_raw Part 2 folds "
                "step4_unpurified.obs onto raw. Re-run the "
                "writeback_to_step1_raw stage (with all prerequisites)."
            )

    spot = _normalize_spot_series(raw_adata.obs["spot_class"])
    first_type = _normalize_spot_series(raw_adata.obs["first_type"])

    n_raw = int(len(spot))
    rctd_mask = spot.ne("").to_numpy()
    n_rctd = int(rctd_mask.sum())
    n_pre_rctd_dropped = n_raw - n_rctd

    spot_rctd = spot[rctd_mask]
    counts = spot_rctd.value_counts()

    ordered_keys = [k for k in _SPOT_CLASS_CANONICAL_ORDER if k in counts.index]
    extras = sorted(v for v in counts.index if v not in _SPOT_CLASS_CANONICAL_ORDER)
    ordered_keys.extend(extras)

    spot_class_rows = []
    for key in ordered_keys:
        c = int(counts[key])
        p = (100.0 * c / n_rctd) if n_rctd else 0.0
        spot_class_rows.append({
            "class": key,
            "display_name": _SPOT_CLASS_DISPLAY.get(key, key),
            "count": c,
            "pct": p,
        })

    reject_labels = _REJECT_VALUES & set(counts.index)
    if reject_labels:
        rejected_mask = spot.isin(reject_labels).to_numpy()
    else:
        import numpy as np

        rejected_mask = np.zeros(n_raw, dtype=bool)
    n_rejected = int(rejected_mask.sum())

    if n_rejected:
        ft_rejected = first_type[rejected_mask]
        ft_counts = ft_rejected.value_counts()
        rejected_first_type_rows = []
        for label, c in ft_counts.items():
            display = label if label else "(none)"
            p = 100.0 * int(c) / n_rejected
            rejected_first_type_rows.append({
                "first_type": display,
                "count": int(c),
                "pct": p,
            })
    else:
        rejected_first_type_rows = []

    return {
        "n_raw": n_raw,
        "n_pre_rctd_dropped": n_pre_rctd_dropped,
        "n_rctd": n_rctd,
        "spot_class_rows": spot_class_rows,
        "n_rejected": n_rejected,
        "rejected_first_type_rows": rejected_first_type_rows,
    }


def _read_hist_thresholds(resolved_yaml: Path) -> dict:
    """Read the three per-histogram filter thresholds from the merged
    ``config.yaml`` so the dashed vertical line on each
    histogram is data-driven (settylab/TracyY123-nexus#26 comment
    5275064492).

    Missing file / missing key → None (line omitted for that
    histogram). Coerces the config value to ``float`` if possible;
    non-numeric config values → None.

    Returns:
      raw:      step1.qc_filter.min_counts_cell (float | None)
      xenium:   None — xenium_ranger.h5ad has no min-counts gate in
                step-1 (xenium_ranger_to_anndata is pass-through).
      purified: step4.postprocess.qc.min_counts (float | None)
    """
    result = {"raw": None, "xenium": None, "purified": None}
    if not resolved_yaml.exists():
        log(f"[qc_report] config.yaml not found at {resolved_yaml} "
            "— histogram threshold lines will be omitted.")
        return result

    import yaml

    try:
        with open(resolved_yaml) as f:
            merged = yaml.safe_load(f) or {}
    except Exception as exc:  # noqa: BLE001 — never fail-loud on a QC-only read
        log(f"[qc_report] failed to parse {resolved_yaml}: {exc} "
            "— histogram threshold lines will be omitted.")
        return result

    def _get(cfg, path):
        cur = cfg
        for k in path:
            if not isinstance(cur, dict) or k not in cur:
                return None
            cur = cur[k]
        try:
            return float(cur)
        except (TypeError, ValueError):
            return None

    result["raw"] = _get(merged, ("step1", "qc_filter", "min_counts_cell"))
    result["purified"] = _get(merged, ("step4", "postprocess", "qc", "min_counts"))
    return result


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def _configure_matplotlib():
    """Import matplotlib with a headless backend, deterministic settings."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["figure.dpi"] = _PLOT_DPI
    plt.rcParams["savefig.dpi"] = _PLOT_DPI
    plt.rcParams["font.family"] = "DejaVu Sans"
    return plt


def _resolve_umap_obsm(purified_adata) -> str:
    """Return the obsm key holding the UMAP embedding. Fail-loud on
    absence — the qc_report stage depends on postprocess having run."""
    candidates = ("X_umap_clipped_norm", "X_umap")
    for key in candidates:
        if key in purified_adata.obsm:
            return key
    raise SystemExit(
        f"[qc_report] proseg_purified.h5ad is missing a UMAP embedding "
        f"(tried {candidates!r} in obsm). Run the `postprocess` stage first."
    )


def _pick_leiden_column(purified_adata) -> tuple[str, str]:
    """Find a leiden column on purified.obs. Returns (obs_key, resolution_str).

    Picks the lowest-numbered leiden_<r> that's present, so the choice is
    deterministic. Fail-loud on absence."""
    leiden_cols = sorted(
        c for c in purified_adata.obs.columns if c.startswith("leiden_")
    )
    if not leiden_cols:
        raise SystemExit(
            "[qc_report] proseg_purified.h5ad has no leiden_* obs column. "
            "Run the `postprocess` stage first."
        )
    obs_key = leiden_cols[0]
    resolution = obs_key[len("leiden_"):]
    return obs_key, resolution


def _resolve_purification_status(purified_adata, column_name: str):
    """Return a string array aligned with purified.obs_names holding
    the SPLIT purification status categorical.

    ``purification_status`` is a SPLIT-native column emitted by
    ``SPLIT::purify`` and preserved on proseg_purified.h5ad by
    ``mtx_to_h5ad``'s obs-cleanup (allow-list membership). It is a
    categorical (e.g. ``'purified'`` / ``'flagged'`` / ``'discarded'``
    or the sample-specific vocabulary SPLIT emits).

    Fail-loud when the column is missing — do NOT try to reindex from
    raw.obs (raw carries the boolean ``passed_purification`` from
    ``writeback_to_step1_raw``, not the SPLIT categorical).
    """
    import numpy as np

    if column_name not in purified_adata.obs.columns:
        raise SystemExit(
            f"[qc_report] proseg_purified.h5ad has no obs[{column_name!r}]. "
            "This is a SPLIT-native column produced by SPLIT::purify; "
            "check that the split_purify → mtx_to_h5ad chain ran "
            "successfully. (Do not confuse with raw.obs['passed_purification'], "
            "which is the boolean written by writeback_to_step1_raw.)"
        )
    log(f"[qc_report]   using purified.obs[{column_name!r}] "
        "for UMAP #1 coloring")
    return np.asarray(
        purified_adata.obs[column_name].astype(str).to_numpy(),
        dtype=object,
    )


def _compute_row_sums(X):
    """Sum each row of a dense or sparse matrix into a 1D float array."""
    import numpy as np
    import scipy.sparse as sp

    if sp.issparse(X):
        return np.asarray(X.sum(axis=1)).ravel().astype(float)
    return np.asarray(X).sum(axis=1).astype(float)


def _hist_log10_counts(
    plt,
    values,
    *,
    title: str,
    xlabel: str,
    out_path: Path,
    threshold: float | None = None,
    n_bins: int = 50,
):
    """Draw a log10 histogram of counts across ALL cells + an optional
    dashed vertical line at the upstream filter threshold.

    ``values`` is a 1D array-like of raw (non-log) counts. Cells with
    non-positive counts are plotted as a separate leftmost
    "did not pass the threshold or generally low count" bar (rather
    than dropped) so the on-disk plot honors Tracy's ask for the
    distribution of ALL cells (settylab/TracyY123-nexus#26 comment
    5275064492; label rename per comment 5275253425).

    ``threshold``: linear-space cutoff; a dashed vertical line is
    drawn at ``log10(threshold)`` with the raw threshold value in the
    legend. Cells with counts strictly below the threshold are counted
    in the caption so the reader can see how many cells the upstream
    filter dropped. ``None`` → no line, no below-threshold count.

    Deterministic — no RNG in binning.
    """
    import numpy as np

    arr = np.asarray(values, dtype=float)
    n_total = int(arr.size)
    finite_mask = np.isfinite(arr)
    positive_mask = finite_mask & (arr > 0)
    n_zero_or_neg = int((~positive_mask).sum())
    pos_vals = arr[positive_mask]
    log_vals = np.log10(pos_vals) if pos_vals.size else np.array([])

    fig, ax = plt.subplots(figsize=_PLOT_FIGSIZE)

    if log_vals.size:
        # Bin edges span the positive-count log range; extend leftward
        # by one bin's worth of width to hold the 0-count bar (drawn as
        # a separate small bar so the leftmost edge is visually
        # distinct from the log10 distribution).
        min_log, max_log = float(log_vals.min()), float(log_vals.max())
        if max_log == min_log:
            max_log = min_log + 1.0
        span = max_log - min_log
        edges = np.linspace(min_log, max_log, n_bins + 1)
        ax.hist(
            log_vals,
            bins=edges,
            color="#4477aa",
            edgecolor="white",
            linewidth=0.4,
        )
        if n_zero_or_neg:
            # Draw the zero/negative cells as a hatched bar half a
            # bin's width to the LEFT of the log range, so it is
            # visually distinct but still part of the same plot.
            bin_width = span / n_bins
            zero_x = min_log - 1.5 * bin_width
            ax.bar(
                zero_x, n_zero_or_neg, width=bin_width,
                color="#bbbbbb", edgecolor="#666666",
                hatch="//", linewidth=0.6,
                label=(
                    "did not pass the threshold or generally low count "
                    f"(n={n_zero_or_neg})"
                ),
            )
    else:
        ax.text(
            0.5, 0.5, "no positive counts",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=11, color="#888",
        )

    n_below_threshold: int | None = None
    if threshold is not None and float(threshold) > 0:
        thr = float(threshold)
        # Cells that FAIL the upstream gate = counts strictly below
        # threshold (matches sc.pp.filter_cells(min_counts=T) which
        # keeps cells with total >= T).
        n_below_threshold = int(
            (positive_mask & (arr < thr)).sum() + n_zero_or_neg
        )
        log_thr = float(np.log10(thr))
        ax.axvline(
            log_thr,
            color="#cc3311", linestyle="--", linewidth=1.5,
            label=(f"filter threshold = {thr:g} "
                   f"(log10 = {log_thr:.2f}); "
                   f"{n_below_threshold} cells below"),
        )

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(loc="best", fontsize=7, frameon=False)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Cells")
    subtitle_parts = [f"n={n_total} cells (all)"]
    if n_zero_or_neg:
        subtitle_parts.append(
            f"{n_zero_or_neg} did not pass the threshold or generally low count"
        )
    ax.set_title(title + "\n" + "; ".join(subtitle_parts), fontsize=10)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(out_path.name + ".tmp.png")
    fig.savefig(tmp, bbox_inches="tight", dpi=_PLOT_DPI, format="png")
    plt.close(fig)
    os.replace(tmp, out_path)


def _scatter_umap(
    plt,
    umap_xy,
    values,
    *,
    title: str,
    out_path: Path,
    categorical: bool,
):
    """Draw a single 2D-scatter UMAP colored by `values` and save to
    `out_path`. `categorical` picks discrete legend rendering."""
    import numpy as np
    import pandas as pd

    fig, ax = plt.subplots(figsize=_PLOT_FIGSIZE)

    if categorical:
        vals = pd.Series(values).astype(str).fillna("")
        uniq = sorted(vals.unique())
        cmap = plt.get_cmap(_PALETTE, max(len(uniq), 1))
        for i, level in enumerate(uniq):
            mask = (vals == level).to_numpy()
            ax.scatter(
                umap_xy[mask, 0], umap_xy[mask, 1],
                s=_PLOT_POINT_SIZE, alpha=_PLOT_ALPHA,
                color=cmap(i),
                label=str(level) if level else "(empty)",
                linewidths=0,
            )
        # Cap legend at 20 entries to keep the plot readable.
        n_legend = min(len(uniq), 20)
        if n_legend > 0:
            ax.legend(
                loc="center left", bbox_to_anchor=(1.02, 0.5),
                fontsize=7, markerscale=2.0, frameon=False,
                ncol=1 if n_legend <= 15 else 2,
            )
    else:
        arr = np.asarray(values, dtype=float)
        sc = ax.scatter(
            umap_xy[:, 0], umap_xy[:, 1],
            c=arr, s=_PLOT_POINT_SIZE, alpha=_PLOT_ALPHA,
            cmap="viridis", linewidths=0,
        )
        fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.04)

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    # Keep the extension on the tmp path so matplotlib can infer the
    # image format from it — `.tmp` alone raises ValueError.
    tmp = out_path.with_name(out_path.name + ".tmp.png")
    fig.savefig(tmp, bbox_inches="tight", dpi=_PLOT_DPI, format="png")
    plt.close(fig)
    os.replace(tmp, out_path)


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------


_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>QC report — {sample_id}</title>
<style>
  body {{ font-family: -apple-system, "Segoe UI", "Helvetica Neue",
         Arial, sans-serif; max-width: 1100px; margin: 2em auto;
         padding: 0 1em; color: #222; }}
  h1 {{ font-size: 1.6em; border-bottom: 2px solid #333;
        padding-bottom: 0.3em; }}
  h2 {{ font-size: 1.2em; color: #333; margin-top: 1.8em; }}
  table {{ border-collapse: collapse; margin: 0.6em 0; font-size: 0.92em; }}
  th, td {{ border: 1px solid #bbb; padding: 4px 10px; text-align: left; }}
  th {{ background: #f0f0f0; }}
  td.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
  .plot {{ margin: 0.8em 0; }}
  .plot img {{ max-width: 780px; border: 1px solid #ddd; }}
  .caption {{ color: #555; font-size: 0.88em; margin: 0.2em 0 0.5em; }}
  .meta {{ color: #666; font-size: 0.85em; }}
  code {{ background: #f4f4f4; padding: 1px 4px; border-radius: 3px; }}
</style>
</head>
<body>
<h1>QC report — {sample_id}</h1>
<p class="meta">
  Run ID: <code>{run_id}</code><br>
  Generated: {generated_ts}<br>
  Pipeline: rctd-split stage <code>qc_report</code>
</p>

<h2>Per-h5ad metrics</h2>
<table>
  <thead><tr>
    <th>h5ad</th><th>Matrix source</th><th>Cells</th><th>Genes</th>
    <th>Median counts/cell</th><th>Mean counts/cell</th>
    <th>Median genes/cell</th>
  </tr></thead>
  <tbody>
    {h5ad_rows}
  </tbody>
</table>

<h2>Purification funnel</h2>
<table>
  <tbody>
    <tr><th>Raw cells</th><td class="num">{n_raw}</td></tr>
    <tr><th>Purified cells</th><td class="num">{n_purified}</td></tr>
    <tr><th>Dropped (raw → purified)</th><td class="num">{n_dropped}</td></tr>
    <tr><th>% purified of raw</th><td class="num">{pct_purified:.2f}%</td></tr>
  </tbody>
</table>

<h2>RCTD summary</h2>
<p class="meta">
  RCTD categorized <b>{n_rctd}</b> cells (of {n_raw_total} raw;
  {n_pre_rctd_dropped} dropped by step-1 QC before RCTD saw them).
  Source: <code>proseg_raw.obs['spot_class']</code> (folded from
  <code>step4_unpurified.obs</code> by
  <code>writeback_to_step1_raw</code>; originally emitted by
  <code>SPLIT::run_post_process_RCTD</code> from RCTD's
  <code>results_df</code>).
</p>
<table>
  <thead><tr>
    <th>spot_class</th><th>Category</th>
    <th>Cells</th><th>% of RCTD-categorized</th>
  </tr></thead>
  <tbody>
    {rctd_spot_class_rows}
  </tbody>
</table>

<h3>SPLIT-inferred celltype in rejected cells</h3>
<p class="meta">
  Of the <b>{n_rejected}</b> cells RCTD rejected
  (<code>spot_class == 'reject'</code>), the primary cell-type call
  (<code>first_type</code>) still assigned by RCTD breaks down as:
</p>
<table>
  <thead><tr>
    <th>SPLIT-inferred celltype</th>
    <th>Cells</th><th>% of rejected</th>
  </tr></thead>
  <tbody>
    {rctd_rejected_first_type_rows}
  </tbody>
</table>

<h2>Count distributions (log10)</h2>
<div class="plot">
  <img src="{img_hist_raw}" alt="log10(proseg_raw total_counts)">
  <div class="caption">
    <code>np.log10(proseg_raw.obs['total_counts'])</code>, over ALL
    raw cells (no positive-only slice). Total counts are the per-cell
    row sum of <code>raw.layers[{raw_layer!r}]</code>
    (proseg's maxpost expression matrix), matching
    <code>obs['total_counts']</code>. {threshold_caption_raw}
  </div>
</div>
<div class="plot">
  <img src="{img_hist_xenium}" alt="log10(xenium_ranger total_counts)">
  <div class="caption">
    <code>np.log10(xenium_ranger.obs['total_counts'])</code>, over ALL
    xenium_ranger cells. {threshold_caption_xenium}
  </div>
</div>
<div class="plot">
  <img src="{img_hist_purified}" alt="log10(proseg_purified nCount_Proseg)">
  <div class="caption">
    <code>np.log10(proseg_purified.obs['nCount_Proseg'])</code>, over
    ALL purified cells — the Seurat-managed per-cell count column.
    {threshold_caption_purified}
  </div>
</div>

<h2>UMAP: purification_status</h2>
<div class="plot">
  <img src="{img_purification_status}" alt="UMAP by purification_status">
  <div class="caption">
    From <code>proseg_purified.h5ad</code>, obs column
    <code>{purification_status_column}</code> (SPLIT-native categorical
    from <code>SPLIT::purify</code>). Distinct labels: {n_purification_status_levels}.
  </div>
</div>

<h2>UMAP: leiden clusters (resolution = {leiden_resolution})</h2>
<div class="plot">
  <img src="{img_leiden}" alt="UMAP by leiden clusters">
  <div class="caption">
    From <code>proseg_purified.h5ad</code>, obs column
    <code>{leiden_obs_key}</code>. Cluster count: {n_leiden_clusters}.
  </div>
</div>

<h2>UMAP: SPLIT-inferred celltype</h2>
<div class="plot">
  <img src="{img_first_type}" alt="UMAP by SPLIT-inferred celltype">
  <div class="caption">
    From <code>proseg_purified.h5ad</code>, obs column
    <code>first_type</code> (RCTD celltype label). Distinct labels:
    {n_first_type_levels}.
  </div>
</div>

<h2>Provenance</h2>
<table>
  <tbody>
    <tr><th>Python</th><td><code>{python_version}</code></td></tr>
    <tr><th>Package versions</th><td><code>{pkg_versions}</code></td></tr>
    <tr><th>Inputs</th><td><code>{input_paths}</code></td></tr>
  </tbody>
</table>
</body>
</html>
"""


def _render_rctd_spot_class_rows(rctd: dict) -> str:
    rows = rctd["spot_class_rows"]
    if not rows:
        return (
            "<tr><td colspan='4' style='color:#888'>"
            "No cells with a non-empty spot_class on raw.obs "
            "(RCTD did not run for this sample).</td></tr>"
        )
    return "\n    ".join(
        (
            "<tr>"
            f"<td><code>{row['class']}</code></td>"
            f"<td>{row['display_name']}</td>"
            f"<td class='num'>{row['count']}</td>"
            f"<td class='num'>{row['pct']:.2f}%</td>"
            "</tr>"
        )
        for row in rows
    )


def _render_rctd_rejected_first_type_rows(rctd: dict) -> str:
    rows = rctd["rejected_first_type_rows"]
    if not rows:
        return (
            "<tr><td colspan='3' style='color:#888'>"
            "No rejected cells to break down.</td></tr>"
        )
    return "\n    ".join(
        (
            "<tr>"
            f"<td><code>{row['first_type']}</code></td>"
            f"<td class='num'>{row['count']}</td>"
            f"<td class='num'>{row['pct']:.2f}%</td>"
            "</tr>"
        )
        for row in rows
    )


def _threshold_caption(source: str, threshold: float | None) -> str:
    if threshold is None:
        return (
            "<i>No upstream min-counts filter recorded in "
            "<code>config.yaml</code>; no dashed threshold "
            "line drawn.</i>"
        )
    return (
        f"Dashed vertical line at the upstream filter threshold "
        f"<code>{source}</code> = <b>{threshold:g}</b> "
        f"(log10 = {float(_log10(threshold)):.2f})."
    )


def _log10(x: float) -> float:
    import math

    return math.log10(float(x))


def _render_html(
    *,
    out_path: Path,
    sample_id: str,
    run_id: str,
    generated_ts: str,
    h5ad_metrics_rows: list[dict],
    purification: dict,
    rctd_summary: dict,
    img_paths: dict,
    thresholds: dict,
    purification_status_column: str,
    n_purification_status_levels: int,
    leiden_resolution: str,
    leiden_obs_key: str,
    n_leiden_clusters: int,
    n_first_type_levels: int,
    raw_layer: str,
    pkg_versions: str,
    input_paths: str,
) -> None:
    rows_html = "\n    ".join(
        (
            "<tr>"
            f"<td>{row['h5ad']}</td>"
            f"<td><code>{row['layer']}</code></td>"
            f"<td class='num'>{row['n_cells']}</td>"
            f"<td class='num'>{row['n_genes']}</td>"
            f"<td class='num'>{row['median_counts_per_cell']:.2f}</td>"
            f"<td class='num'>{row['mean_counts_per_cell']:.2f}</td>"
            f"<td class='num'>{row['median_genes_per_cell']:.2f}</td>"
            "</tr>"
        )
        for row in h5ad_metrics_rows
    )
    body = _HTML_TEMPLATE.format(
        sample_id=sample_id,
        run_id=run_id,
        generated_ts=generated_ts,
        h5ad_rows=rows_html,
        n_raw=purification["n_raw"],
        n_purified=purification["n_purified"],
        n_dropped=purification["n_dropped_raw_to_purified"],
        pct_purified=purification["pct_purified_of_raw"],
        n_raw_total=rctd_summary["n_raw"],
        n_pre_rctd_dropped=rctd_summary["n_pre_rctd_dropped"],
        n_rctd=rctd_summary["n_rctd"],
        n_rejected=rctd_summary["n_rejected"],
        rctd_spot_class_rows=_render_rctd_spot_class_rows(rctd_summary),
        rctd_rejected_first_type_rows=(
            _render_rctd_rejected_first_type_rows(rctd_summary)
        ),
        threshold_caption_raw=_threshold_caption(
            "step1.qc_filter.min_counts_cell", thresholds["raw"],
        ),
        threshold_caption_xenium=_threshold_caption(
            "step1 (xenium_ranger min-counts)", thresholds["xenium"],
        ),
        threshold_caption_purified=_threshold_caption(
            "step4.postprocess.qc.min_counts", thresholds["purified"],
        ),
        img_hist_raw=img_paths["hist_raw"],
        img_hist_xenium=img_paths["hist_xenium"],
        img_hist_purified=img_paths["hist_purified"],
        img_purification_status=img_paths["purification_status"],
        img_leiden=img_paths["leiden"],
        img_first_type=img_paths["first_type"],
        purification_status_column=purification_status_column,
        n_purification_status_levels=n_purification_status_levels,
        leiden_resolution=leiden_resolution,
        leiden_obs_key=leiden_obs_key,
        n_leiden_clusters=n_leiden_clusters,
        n_first_type_levels=n_first_type_levels,
        raw_layer=raw_layer,
        python_version=sys.version.splitlines()[0],
        pkg_versions=pkg_versions,
        input_paths=input_paths,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    tmp.write_text(body)
    os.replace(tmp, out_path)


def _pkg_versions() -> str:
    parts = []
    for name in ("numpy", "pandas", "scipy", "scanpy", "anndata", "matplotlib"):
        try:
            mod = __import__(name)
            parts.append(f"{name}={getattr(mod, '__version__', '?')}")
        except Exception:
            parts.append(f"{name}=missing")
    return ", ".join(parts)


def _write_metrics_csv(
    csv_path: Path,
    h5ad_metrics_rows: list[dict],
    purification: dict,
) -> None:
    import pandas as pd

    df = pd.DataFrame(h5ad_metrics_rows)
    for k, v in purification.items():
        df[k] = v
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = csv_path.with_suffix(csv_path.suffix + ".tmp")
    df.to_csv(tmp, index=False)
    os.replace(tmp, csv_path)


def _write_rctd_summary_csv(csv_path: Path, rctd: dict) -> None:
    """Machine-parseable sidecar for the RCTD summary. Two logical
    tables concatenated with a ``section`` discriminator column:

      section=spot_class:    class, count, pct  (of RCTD-categorized)
      section=rejected_first_type: first_type, count, pct  (of rejected)

    Plus a leading provenance row (section=summary) with the raw
    denominators so downstream consumers do not need to recompute
    percentages if they want to combine samples.
    """
    import pandas as pd

    rows = [
        {
            "section": "summary",
            "key": "n_raw", "value": rctd["n_raw"],
            "count": rctd["n_raw"], "pct": 100.0,
        },
        {
            "section": "summary",
            "key": "n_pre_rctd_dropped", "value": rctd["n_pre_rctd_dropped"],
            "count": rctd["n_pre_rctd_dropped"],
            "pct": (100.0 * rctd["n_pre_rctd_dropped"] / rctd["n_raw"])
            if rctd["n_raw"] else 0.0,
        },
        {
            "section": "summary",
            "key": "n_rctd", "value": rctd["n_rctd"],
            "count": rctd["n_rctd"],
            "pct": (100.0 * rctd["n_rctd"] / rctd["n_raw"])
            if rctd["n_raw"] else 0.0,
        },
        {
            "section": "summary",
            "key": "n_rejected", "value": rctd["n_rejected"],
            "count": rctd["n_rejected"],
            "pct": (100.0 * rctd["n_rejected"] / rctd["n_rctd"])
            if rctd["n_rctd"] else 0.0,
        },
    ]
    for r in rctd["spot_class_rows"]:
        rows.append({
            "section": "spot_class",
            "key": r["class"], "value": r["display_name"],
            "count": r["count"], "pct": r["pct"],
        })
    for r in rctd["rejected_first_type_rows"]:
        rows.append({
            "section": "rejected_first_type",
            "key": r["first_type"], "value": r["first_type"],
            "count": r["count"], "pct": r["pct"],
        })

    df = pd.DataFrame(rows, columns=["section", "key", "value", "count", "pct"])
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = csv_path.with_suffix(csv_path.suffix + ".tmp")
    df.to_csv(tmp, index=False)
    os.replace(tmp, csv_path)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def _hist_values_from_obs(
    adata, obs_col: str, label: str,
):
    """Look up an obs count column; fail-loud if absent."""
    if obs_col not in adata.obs.columns:
        raise SystemExit(
            f"[qc_report] {label}.h5ad is missing obs[{obs_col!r}] "
            "— cannot draw the log10 count histogram. Check the "
            "upstream stage that populates this column."
        )
    import numpy as np
    return np.asarray(adata.obs[obs_col].to_numpy(), dtype=float)


def run_qc_report(
    sample_id: str,
    run_id: str,
    output_root: Path,
    purification_status_column: str,
    raw_layer: str,
    force_rerun: bool,
) -> Path:
    """Generate QC HTML report + plots. Returns the sentinel path.

    Reads-only w.r.t. the three source h5ads. Writes to
    ``<run_dir>/summary/`` + a sentinel under ``intermediate/adata/``.
    """
    import time

    import anndata as ad

    xenium_p = spatial_adata_path(
        output_root, sample_id, run_id, "xenium_ranger",
    )
    raw_p = spatial_adata_path(
        output_root, sample_id, run_id, "proseg_raw",
    )
    purified_p = spatial_adata_path(
        output_root, sample_id, run_id, "proseg_purified",
    )
    sentinel = intermediate_path(
        output_root, sample_id, run_id, "qc_report_sentinel",
    )
    metrics_csv = summary_path(output_root, sample_id, run_id, "metrics_csv")
    rctd_summary_csv = summary_path(
        output_root, sample_id, run_id, "rctd_summary_csv",
    )
    html_out = summary_path(output_root, sample_id, run_id, "html_report")
    plots_dir = summary_plots_dir(output_root, sample_id, run_id)
    resolved_yaml = resolved_config_path(output_root, sample_id, run_id)

    if sentinel_exists(sentinel, force_rerun):
        log(f"[qc_report] sentinel exists: {sentinel} — skipping "
            "(pass --force-rerun to re-run).")
        return sentinel

    for name, p in (
        ("proseg_purified", purified_p),
        ("proseg_raw", raw_p),
        ("xenium_ranger", xenium_p),
    ):
        if not p.exists():
            raise SystemExit(
                f"[qc_report] required input missing: {name} h5ad not found "
                f"at {p}. Ensure the full pipeline through celltype_writeback "
                "has run."
            )

    log(f"[qc_report] reading {xenium_p}")
    xenium = ad.read_h5ad(xenium_p)
    log(f"[qc_report] reading {raw_p}")
    raw = ad.read_h5ad(raw_p)
    log(f"[qc_report] reading {purified_p}")
    purified = ad.read_h5ad(purified_p)

    log(f"[qc_report] computing per-h5ad metrics "
        f"(raw uses layers[{raw_layer!r}])")
    h5ad_metrics_rows = [
        _h5ad_metrics(xenium, "xenium_ranger"),
        _h5ad_metrics(raw, "proseg_raw", layer=raw_layer),
        _h5ad_metrics(purified, "proseg_purified"),
    ]
    purification = _purification_metrics(raw, purified)
    log(f"[qc_report] purification funnel: "
        f"{purification['n_raw']} raw → "
        f"{purification['n_purified']} purified "
        f"({purification['pct_purified_of_raw']:.2f}%)")

    log("[qc_report] computing RCTD spot_class summary + "
        "first_type-in-rejected tabulation from raw.obs")
    rctd_summary = _rctd_summary_metrics(raw)
    log(f"[qc_report] RCTD summary: {rctd_summary['n_rctd']} categorized "
        f"({rctd_summary['n_pre_rctd_dropped']} step-1-dropped), "
        f"{rctd_summary['n_rejected']} rejected")

    log(f"[qc_report] reading histogram thresholds from {resolved_yaml}")
    thresholds = _read_hist_thresholds(resolved_yaml)
    log(f"[qc_report] thresholds: raw={thresholds['raw']}, "
        f"xenium={thresholds['xenium']}, purified={thresholds['purified']}")

    # Resolve plot inputs (fail-loud on missing).
    umap_key = _resolve_umap_obsm(purified)
    umap_xy = purified.obsm[umap_key]

    purification_status_values = _resolve_purification_status(
        purified, purification_status_column,
    )

    leiden_obs_key, leiden_resolution = _pick_leiden_column(purified)
    leiden_values = purified.obs[leiden_obs_key].astype(str).to_numpy()

    if "first_type" not in purified.obs.columns:
        raise SystemExit(
            "[qc_report] proseg_purified.h5ad is missing obs['first_type'] "
            "— cannot produce the first_type UMAP plot."
        )
    first_type_values = purified.obs["first_type"].astype(str).to_numpy()

    # Count-histogram inputs (fail-loud on missing obs columns).
    hist_raw_vals = _hist_values_from_obs(
        raw, "total_counts", "proseg_raw",
    )
    hist_xenium_vals = _hist_values_from_obs(
        xenium, "total_counts", "xenium_ranger",
    )
    hist_purified_vals = _hist_values_from_obs(
        purified, "nCount_Proseg", "proseg_purified",
    )

    # Also touch layers[raw_layer] to enforce the fail-loud invariant
    # for the histogram source; the layer is guaranteed to exist here
    # (metrics-row assembly above already resolved it).
    _resolve_matrix(raw, raw_layer, "proseg_raw")

    # Emit plots.
    plt = _configure_matplotlib()
    plots_dir.mkdir(parents=True, exist_ok=True)
    img_purification = (
        plots_dir / f"{sample_id}_umap_purification_status.png"
    )
    img_leiden = plots_dir / f"{sample_id}_umap_leiden_res{leiden_resolution}.png"
    img_first_type = plots_dir / f"{sample_id}_umap_first_type.png"
    img_hist_raw = (
        plots_dir / f"{sample_id}_hist_total_counts_proseg_raw.png"
    )
    img_hist_xenium = (
        plots_dir / f"{sample_id}_hist_total_counts_xenium_ranger.png"
    )
    img_hist_purified = (
        plots_dir / f"{sample_id}_hist_nCount_Proseg_proseg_purified.png"
    )

    log(f"[qc_report] writing plot {img_purification}")
    _scatter_umap(
        plt, umap_xy, purification_status_values,
        title=f"{sample_id}: purification_status",
        out_path=img_purification,
        categorical=True,
    )
    log(f"[qc_report] writing plot {img_leiden}")
    _scatter_umap(
        plt, umap_xy, leiden_values,
        title=f"{sample_id}: leiden (resolution={leiden_resolution})",
        out_path=img_leiden,
        categorical=True,
    )
    log(f"[qc_report] writing plot {img_first_type}")
    _scatter_umap(
        plt, umap_xy, first_type_values,
        title=f"{sample_id}: SPLIT-inferred celltype",
        out_path=img_first_type,
        categorical=True,
    )
    log(f"[qc_report] writing histogram {img_hist_raw}")
    _hist_log10_counts(
        plt, hist_raw_vals,
        title=f"{sample_id}: proseg_raw log10(total_counts)",
        xlabel="log10(total_counts)",
        out_path=img_hist_raw,
        threshold=thresholds["raw"],
    )
    log(f"[qc_report] writing histogram {img_hist_xenium}")
    _hist_log10_counts(
        plt, hist_xenium_vals,
        title=f"{sample_id}: xenium_ranger log10(total_counts)",
        xlabel="log10(total_counts)",
        out_path=img_hist_xenium,
        threshold=thresholds["xenium"],
    )
    log(f"[qc_report] writing histogram {img_hist_purified}")
    _hist_log10_counts(
        plt, hist_purified_vals,
        title=f"{sample_id}: proseg_purified log10(nCount_Proseg)",
        xlabel="log10(nCount_Proseg)",
        out_path=img_hist_purified,
        threshold=thresholds["purified"],
    )

    # Metrics CSV (sidecar, machine-parseable).
    _write_metrics_csv(metrics_csv, h5ad_metrics_rows, purification)
    log(f"[qc_report] wrote metrics {metrics_csv}")
    _write_rctd_summary_csv(rctd_summary_csv, rctd_summary)
    log(f"[qc_report] wrote RCTD summary {rctd_summary_csv}")

    # HTML report with relative image paths (so the directory is
    # portable — you can rsync `summary/` anywhere and it renders).
    _render_html(
        out_path=html_out,
        sample_id=sample_id,
        run_id=run_id,
        generated_ts=time.strftime("%Y-%m-%dT%H:%M:%S"),
        h5ad_metrics_rows=h5ad_metrics_rows,
        purification=purification,
        rctd_summary=rctd_summary,
        thresholds=thresholds,
        img_paths={
            "purification_status": f"plots/{img_purification.name}",
            "leiden": f"plots/{img_leiden.name}",
            "first_type": f"plots/{img_first_type.name}",
            "hist_raw": f"plots/{img_hist_raw.name}",
            "hist_xenium": f"plots/{img_hist_xenium.name}",
            "hist_purified": f"plots/{img_hist_purified.name}",
        },
        purification_status_column=purification_status_column,
        n_purification_status_levels=int(
            len(set(purification_status_values))
        ),
        leiden_resolution=leiden_resolution,
        leiden_obs_key=leiden_obs_key,
        n_leiden_clusters=int(len(set(leiden_values))),
        n_first_type_levels=int(len(set(first_type_values))),
        raw_layer=raw_layer,
        pkg_versions=_pkg_versions(),
        input_paths=(
            f"xenium={xenium_p}, raw={raw_p}, purified={purified_p}"
        ),
    )
    log(f"[qc_report] wrote HTML {html_out}")

    sentinel.parent.mkdir(parents=True, exist_ok=True)
    sentinel.write_text(
        f"qc_report complete: {sample_id}\n"
        f"n_raw={purification['n_raw']} "
        f"n_purified={purification['n_purified']} "
        f"leiden_resolution={leiden_resolution} "
        f"raw_layer={raw_layer}\n"
    )
    log(f"[qc_report] wrote sentinel {sentinel}")
    return sentinel
