"""Standalone, iteration-only QC summary — NOT wired into the pipeline
and NOT imported by ``pipeline.py`` or ``qc_report.py``.

This module exists purely so the operator can react to a candidate
layout before any of it is folded into the real ``qc_report`` stage
(settylab/msetty-nexus#41). It never writes to a run directory — it
takes an explicit ``--out`` path — and it only *reads*: the three
existing ``spatial_adata/*.h5ad`` files, the PNGs
``qc_report.run_qc_report`` already produced under ``summary/plots/``,
and the existing ``summary_report.html`` text (to lift its tables /
captions verbatim — see "Superset strategy" below). It does not touch
``qc_report.py``, does not recompute anything ``qc_report`` already
computed, and does not regenerate or modify the existing
``summary_report.html``.

**This is a superset of the existing report, not a replacement — and
its section ORDER matches ``summary_report.html`` exactly.** An
earlier draft led with the spatial panels; the operator corrected that
(settylab/msetty-nexus#41) — the instruction was always "original
report + additions", not a re-architecture. Section list, in the
document's actual top-to-bottom order:

  1. Per-h5ad metrics / Purification funnel / RCTD summary / SPLIT-
     inferred celltype in rejected cells — the four existing tables,
     lifted verbatim, same position as the original (identical
     numbers — same source h5ads, unchanged since the original run).
  2. **[CHANGED LAYOUT]** Count distributions (log10) — same position
     as the original, same three PNGs, laid out in a row instead of
     stacked.
  3. **[CHANGED LAYOUT]** UMAP: purification_status / leiden clusters /
     SPLIT-inferred celltype — same position as the original (right
     after the histograms), same three images + captions, laid out in
     a row instead of one after another. Each panel's original ``<h2>``
     heading text is preserved as a small in-panel title (not deleted,
     just no longer a full-width heading).
  4. Provenance — lifted verbatim, same position (last of the original
     content). Describes how the ORIGINAL pipeline run was produced
     (invocation argv, package versions at run time, resolved config),
     which this script cannot regenerate itself (it's an artifact of
     that run, not of this script's own invocation) — copying it
     forward is the only faithful option.
  5. **[NEW, appended at the end]** Spatial overview — broad leiden
     clusters + broad cell type panels, sized/colored for a two-second
     glance. HE is not available at this pipeline stage, so these
     stand alone.
  6. **[NEW, appended at the end]** QC metrics in space — small spatial
     mirrors of the three count histograms (proseg_raw total_counts,
     xenium_ranger total_counts, proseg_purified nCount_Proseg). No
     cell-type panel here — cell type already has its own broad panel
     above (#5), repeating it here was redundant. These three panels
     clip their COLOR SCALE (not the underlying data — every cell
     still renders) to each panel's own 1st/99th percentile, so a
     handful of extreme cells can't flatten the visible variation;
     each panel says so in its title + colorbar label.

Superset strategy: rather than re-deriving the tables/captions/UMAPs
from the h5ads a second time (risking silent divergence from what
``qc_report`` actually computed), this script parses the ALREADY
-RENDERED ``summary_report.html`` and lifts the relevant HTML
fragments verbatim, rewriting each fragment's ``src="plots/....png"``
references to base64 data URIs. This guarantees byte-identical numbers
and wording to the existing report for every carried-over section, and
means this script only needs to touch the h5ads for the three NEW
spatial panels (obs/obsm only, ``backed="r"`` — no full matrix load).

Data columns used for the NEW panels (state exactly what was read):

  * ``<S>_proseg_purified.h5ad``:
      - ``obsm['spatial']``       -- cell centroid (x, y), xenium frame
      - ``obs['leiden_0.5']``     -- lowest-resolution leiden clustering
        (same column ``qc_report`` picks via ``_pick_leiden_column``)
      - ``obs['first_type']``     -- SPLIT/RCTD-inferred cell type
      - ``obs['nCount_Proseg']``  -- per-cell total counts (Seurat-managed)
  * ``<S>_proseg_raw.h5ad``:      ``obsm['spatial']``, ``obs['total_counts']``
  * ``<S>_xenium_ranger.h5ad``:   ``obsm['spatial']``, ``obs['total_counts']``

All three h5ads carry ``obsm['spatial']`` in the same coordinate frame
(cross-checked: x/y ranges agree within the segmentation-method
cell-count difference), so the small per-histogram spatial mirrors are
directly comparable panel-to-panel.

Orientation: every spatial panel calls ``ax.invert_yaxis()`` (Y only —
X is left alone) so the rendered scatter follows image convention
(origin top-left, y increasing downward), matching how the eventual
H&E overlay will be displayed (settylab/msetty-nexus#41,
``heRegistration/src/hexenium/stages/viz.py`` — ``origin="upper"`` +
``ax.set_ylim(h, 0)``). This was verified EMPIRICALLY, not assumed:
rasterized this run's ``obsm['spatial']`` through the actual
``ax.scatter`` + ``invert_yaxis`` code path at the same resolution as
this run's own native Xenium image
(``he_registration/register/*/masks/morphology_focus_0000.png`` — the
un-warped Xenium morphology image, in the same frame as
``obsm['spatial']``, with no H&E-registration transform in the way)
and swept all four X/Y flip combinations, scoring each by
intersection-over-union against the morphology image's tissue mask.
``invert_y=True, invert_x=False`` (what's implemented) wins outright
(IoU 0.773 vs 0.689 for no inversion, 0.60-0.62 for the two
X-inverted variants) — a decisive, unambiguous margin, confirmed by
a pixel-level overlay showing scatter and morphology-image tissue
boundary lining up edge-for-edge. Do not "fix" this again without
re-running that check — the plausible-sounding alternative (add an
inversion because a naive obsm['spatial']-is-image-row-order
assumption suggests one is missing) was tested FIRST and empirically
LOST; a first pass at this check used a raw-numpy-array rasterization
helper that silently used the opposite row/y convention from
matplotlib's default and produced a misleading result — the fix was
to test the *actual* ``ax.scatter``/``invert_yaxis`` rendering path
directly against ground truth, not reason about conventions on paper.

Point size / alpha were tuned by visual inspection against this run's
146k-cell purified population (settylab/msetty-nexus#36's lesson: a
plot that renders sub-pixel is not "broad", it's invisible) — broad
panels use ``s=2.5, alpha=1.0`` (unchanged, not touched by the small
QC-panel size tuning below). The smaller per-metric panels use
``s=0.4, alpha=1.0`` — reduced from the earlier ``s=1.0`` per operator
request (settylab/msetty-nexus#41); these panels are already narrower
than the broad ones, so a smaller marker keeps them from
over-saturating into flat color blocks. Both settings show clear
tissue-scale domains (e.g. the liver rim, the fibrotic band) at the
rendered size.

Usage::

    python -m rctd_split.stages.qc_report_spatial_v2 \\
        --sample-id MH8_S24-158053-01_BR_2453 \\
        --run-dir /path/to/MH8_S24-158053-01_BR_2453_env_config_v1 \\
        --out /path/to/out.html
"""
from __future__ import annotations

import argparse
import base64
import html
import re
import time
from pathlib import Path


_PLOT_DPI = 150
_BROAD_FIG_WIDTH = 7.0
_SMALL_FIG_WIDTH = 3.4
_BROAD_POINT_SIZE = 2.5
_SMALL_POINT_SIZE = 0.4
_POINT_ALPHA = 1.0
_CLIP_PERCENTILES = (1.0, 99.0)


def _configure_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["figure.dpi"] = _PLOT_DPI
    plt.rcParams["savefig.dpi"] = _PLOT_DPI
    plt.rcParams["font.family"] = "DejaVu Sans"
    return plt


def _png_to_data_uri(path: Path) -> str:
    raw = path.read_bytes()
    b64 = base64.b64encode(raw).decode("ascii")
    return f"data:image/png;base64,{b64}"


def _fig_to_data_uri(fig) -> str:
    import io

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=_PLOT_DPI, bbox_inches="tight")
    b64 = base64.b64encode(buf.getvalue()).decode("ascii")
    return f"data:image/png;base64,{b64}"


def _spatial_extent(xy):
    x_span = float(xy[:, 0].max() - xy[:, 0].min())
    y_span = float(xy[:, 1].max() - xy[:, 1].min())
    return x_span, y_span


def _broad_spatial_categorical(
    plt, xy, values, *, color_map, title, legend_title,
):
    import numpy as np

    x_span, y_span = _spatial_extent(xy)
    aspect = y_span / x_span if x_span else 1.0
    fig_h = _BROAD_FIG_WIDTH * aspect
    fig, ax = plt.subplots(figsize=(_BROAD_FIG_WIDTH, fig_h))

    vals = np.asarray(values, dtype=object).astype(str)
    uniq = sorted(set(vals))
    for level in uniq:
        mask = vals == level
        color = color_map.get(level, "#BBBBBB") if color_map else None
        ax.scatter(
            xy[mask, 0], xy[mask, 1],
            s=_BROAD_POINT_SIZE, alpha=_POINT_ALPHA,
            color=color, linewidths=0, label=level,
        )
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
    ax.set_title(title, fontsize=13)
    n_legend = min(len(uniq), 20)
    ax.legend(
        loc="upper left", bbox_to_anchor=(1.01, 1.0),
        fontsize=9, markerscale=4.0, frameon=False,
        title=legend_title, ncol=1 if n_legend <= 20 else 2,
    )
    uri = _fig_to_data_uri(fig)
    plt.close(fig)
    return uri


def _small_spatial_continuous(
    plt, xy, values, *, title, clip_percentiles=_CLIP_PERCENTILES,
):
    """Spatial scatter colored by a continuous metric (log10 scale).

    The COLOR SCALE (``vmin``/``vmax``) is clipped to the
    ``clip_percentiles`` of THIS panel's own values — computed
    independently per panel, never shared across metrics, and computed
    only over cells with a finite value (``np.nanpercentile`` — some
    metrics, e.g. ``proseg_raw.obs['total_counts']``, carry NaN for
    cells that never reached the stage that populates them; a plain
    ``np.percentile`` propagates those NaNs into ``vmin``/``vmax`` and
    silently degenerates the whole color scale to a single flat color,
    caught by visual inspection while building this panel). Clipping
    the color scale never drops a cell from the scatter: every
    finite-value cell is plotted, points beyond the clip range simply
    saturate at the end color. Cells with no finite value at all
    (distinct from "clipped out of range") still get a marker, in a
    neutral gray, so the tissue outline stays complete. The title and
    colorbar label both say the scale is clipped so a reader can't
    mistake it for the true range.
    """
    import numpy as np

    x_span, y_span = _spatial_extent(xy)
    aspect = y_span / x_span if x_span else 1.0
    fig_h = _SMALL_FIG_WIDTH * aspect
    fig, ax = plt.subplots(figsize=(_SMALL_FIG_WIDTH, fig_h))

    arr = np.asarray(values, dtype=float)
    log_arr = np.log10(np.clip(arr, 1.0, None))
    finite_mask = np.isfinite(log_arr)

    if finite_mask.any():
        ax.scatter(
            xy[~finite_mask, 0], xy[~finite_mask, 1],
            s=_SMALL_POINT_SIZE, alpha=_POINT_ALPHA,
            color="#dddddd", linewidths=0,
        )
        lo, hi = clip_percentiles
        vmin, vmax = np.nanpercentile(log_arr[finite_mask], [lo, hi])
        if not (vmax > vmin):
            vmin = vmax = None  # degenerate distribution — fall back to auto-scale
        sc = ax.scatter(
            xy[finite_mask, 0], xy[finite_mask, 1], c=log_arr[finite_mask],
            s=_SMALL_POINT_SIZE, alpha=_POINT_ALPHA,
            cmap="viridis", linewidths=0,
            vmin=vmin, vmax=vmax,
        )
        fig.colorbar(
            sc, ax=ax, fraction=0.045, pad=0.04,
            label=f"log10 (color clipped {lo:g}\N{EN DASH}{hi:g}%ile)",
        )
        n_missing = int((~finite_mask).sum())
        subtitle = f"\n(color clipped to {lo:g}\N{EN DASH}{hi:g}th %ile"
        subtitle += f"; {n_missing} cells no value)" if n_missing else ")"
    else:
        ax.scatter(
            xy[:, 0], xy[:, 1], s=_SMALL_POINT_SIZE, alpha=_POINT_ALPHA,
            color="#dddddd", linewidths=0,
        )
        subtitle = "\n(no finite values)"

    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title + subtitle, fontsize=8)
    uri = _fig_to_data_uri(fig)
    plt.close(fig)
    return uri


# ---------------------------------------------------------------------------
# Lifting sections verbatim from the existing rendered summary_report.html
# ---------------------------------------------------------------------------


def _slice_between(text: str, start_marker: str, end_marker: str | None) -> str:
    start = text.index(start_marker)
    if end_marker is None:
        return text[start:]
    end = text.index(end_marker, start)
    return text[start:end]


_IMG_SRC_RE = re.compile(r'src="plots/([^"]+\.png)"')


_UMAP_HEADINGS: tuple[str, ...] = (
    "UMAP: purification_status",
    "UMAP: leiden clusters (resolution = 0.5)",
    "UMAP: SPLIT-inferred celltype",
)
_UMAP_ROW_HEIGHT_PX = 380


def _build_umap_row(text: str, plots_dir: Path) -> str:
    """Reassemble the three individual UMAP sections (each its own
    ``<h2>`` + ``<div class="plot">`` + caption in the existing report)
    into one horizontal row — same position in the document as the
    existing report's UMAP sections, just laid out side by side
    instead of stacked. Each section's original ``<h2>`` text is kept,
    demoted from a page-level heading into a small in-panel title
    (content preserved verbatim, only the layout changes).

    Sets explicit pixel ``width``/``height`` attributes on each ``img``
    (read from the actual PNG via Pillow) rather than relying on
    CSS ``height: Npx; width: auto`` inside a ``display: table-cell``
    column. The latter was tried first and broke: the leiden panel (a
    much narrower source PNG paired with the row's longest title text)
    rendered squashed into a slim column while its siblings rendered
    full width — an auto-table-layout column-width quirk, caught by
    rendering and visually comparing all three panels side by side,
    not assumed correct from a clean run. Explicit pixel dimensions
    sidestep the browser's column-width guess entirely.
    """
    from PIL import Image

    markers = [f"<h2>{h}</h2>" for h in _UMAP_HEADINGS] + ["<h2>Provenance</h2>"]
    panels = []
    for i, heading_text in enumerate(_UMAP_HEADINGS):
        chunk = _slice_between(text, markers[i], markers[i + 1])
        rest = chunk.split(markers[i], 1)[1]  # drop the heading tag itself
        m = _IMG_SRC_RE.search(rest)
        if m is None:
            raise SystemExit(
                f"[qc_report_spatial_v2] no plots/*.png <img> found under "
                f"heading {heading_text!r} in the existing report — cannot "
                "lift this panel verbatim."
            )
        png_path = plots_dir / m.group(1)
        if not png_path.exists():
            raise SystemExit(
                f"[qc_report_spatial_v2] existing report references "
                f"{png_path} but it does not exist."
            )
        with Image.open(png_path) as im:
            w_px, h_px = im.size
        disp_w = round(_UMAP_ROW_HEIGHT_PX * w_px / h_px)
        data_uri = _png_to_data_uri(png_path)
        rest = _IMG_SRC_RE.sub(
            f'src="{data_uri}" width="{disp_w}" height="{_UMAP_ROW_HEIGHT_PX}"',
            rest, count=1,
        )
        title_html = f'<div class="umap-title">{html.escape(heading_text)}</div>'
        rest = rest.replace('<div class="plot">', f'<div class="plot">{title_html}', 1)
        panels.append(rest.strip())
    return '<div class="row umap">\n' + "\n".join(panels) + "\n</div>"


def _inline_plot_images(fragment: str, plots_dir: Path) -> str:
    """Rewrite every ``src="plots/<name>.png"`` in ``fragment`` to a
    base64 data URI, reading the PNG from ``plots_dir``. Fail-loud if a
    referenced PNG is missing — a silently-dropped image would ship a
    broken report."""

    def _sub(m: re.Match) -> str:
        png_path = plots_dir / m.group(1)
        if not png_path.exists():
            raise SystemExit(
                f"[qc_report_spatial_v2] existing report references "
                f"{png_path} but it does not exist — cannot lift this "
                "section verbatim."
            )
        return f'src="{_png_to_data_uri(png_path)}"'

    return _IMG_SRC_RE.sub(_sub, fragment)


def _lift_existing_sections(existing_html_path: Path, plots_dir: Path) -> dict:
    """Parse the existing ``summary_report.html`` and return the four
    carried-over fragments (image-inlined where relevant), keyed for
    the new template. Fail-loud if the expected section markers are
    not found — the existing report's structure changed and this
    extraction needs re-checking rather than silently emitting an
    incomplete page."""
    if not existing_html_path.exists():
        raise SystemExit(
            f"[qc_report_spatial_v2] existing report not found at "
            f"{existing_html_path} — run the `qc_report` pipeline stage "
            "first (this script only adds to its output, never "
            "regenerates it)."
        )
    text = existing_html_path.read_text()

    required_markers = (
        "<h2>Per-h5ad metrics</h2>",
        "<h2>Count distributions (log10)</h2>",
        *(f"<h2>{h}</h2>" for h in _UMAP_HEADINGS),
        "<h2>Provenance</h2>",
        "</body>",
    )
    for marker in required_markers:
        if marker not in text:
            raise SystemExit(
                f"[qc_report_spatial_v2] expected marker {marker!r} not "
                f"found in {existing_html_path} — the existing report's "
                "structure has changed; re-check the section-lifting "
                "logic in qc_report_spatial_v2.py before proceeding."
            )

    tabular_blob = _slice_between(
        text, "<h2>Per-h5ad metrics</h2>", "<h2>Count distributions (log10)</h2>",
    )
    # Slice starting AFTER the heading — our own template supplies its own
    # "Count distributions (log10)" <h2> immediately before this blob, so
    # keeping the lifted heading here would duplicate it.
    hist_blob = _slice_between(
        text, "<h2>Count distributions (log10)</h2>", "<h2>UMAP: purification_status</h2>",
    ).split("<h2>Count distributions (log10)</h2>", 1)[1]
    umap_row = _build_umap_row(text, plots_dir)
    provenance_blob = _slice_between(text, "<h2>Provenance</h2>", "</body>")

    return {
        "tabular": tabular_blob,
        "hist": _inline_plot_images(hist_blob, plots_dir),
        "umap_row": umap_row,
        "provenance": provenance_blob,
    }


_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>QC summary (v2) — {sample_id}</title>
<style>
  body {{ font-family: -apple-system, "Segoe UI", "Helvetica Neue",
         Arial, sans-serif; max-width: 1500px; margin: 2em auto;
         padding: 0 1em; color: #222; }}
  h1 {{ font-size: 1.6em; border-bottom: 2px solid #333;
        padding-bottom: 0.3em; }}
  h2 {{ font-size: 1.2em; color: #333; margin-top: 1.8em; }}
  h3 {{ font-size: 1.05em; color: #333; margin-top: 1.4em; }}
  table {{ border-collapse: collapse; margin: 0.6em 0; font-size: 0.92em; }}
  th, td {{ border: 1px solid #bbb; padding: 4px 10px; text-align: left; }}
  th {{ background: #f0f0f0; }}
  td.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
  .meta {{ color: #666; font-size: 0.85em; }}
  .row {{ display: table; table-layout: auto; margin: 0.8em 0;
          border-spacing: 1em 0; border-collapse: separate; }}
  .row .plot {{ display: table-cell; vertical-align: top; margin: 0; }}
  .plot img {{ max-width: 100%; border: 1px solid #ddd; }}
  .row.hist .plot img {{ height: 320px; width: auto; }}
  .row.broad .plot img {{ height: 640px; width: auto; }}
  .row.small .plot img {{ height: 340px; width: auto; }}
  .row.umap .plot {{ max-width: none; }}
  .row.umap .plot img {{ max-width: none; }}
  .umap-title {{ font-weight: 600; font-size: 0.95em; margin: 0 0 0.3em; }}
  .caption {{ color: #555; font-size: 0.85em; margin: 0.3em 0 0; }}
  code {{ background: #f4f4f4; padding: 1px 4px; border-radius: 3px; }}
  .banner {{ background: #fff8e1; border: 1px solid #e6c200;
             padding: 0.6em 1em; border-radius: 4px; font-size: 0.9em; }}
  details summary {{ cursor: pointer; }}
  pre {{ background: #f8f8f8; padding: 0.8em; overflow-x: auto;
         font-size: 0.82em; }}
</style>
</head>
<body>
<h1>QC summary (v2) — {sample_id}</h1>
<p class="banner">
  Draft layout for review, generated outside the pipeline
  (<code>qc_report_spatial_v2.py</code>, not wired into any stage).
  <b>Superset</b> of <code>{sample_id}_summary_report.html</code>, same
  section order as that file: tables, histograms (now a row), UMAPs
  (now a row), Provenance — then two new spatial sections appended at
  the end. That file is untouched. Tracking: settylab/msetty-nexus#41.
</p>
<p class="meta">
  Run: <code>{run_id}</code><br>
  Generated: {generated_ts}
</p>

{tabular_blob}

<h2>Count distributions (log10)</h2>
<div class="row hist">
{hist_blob}
</div>

{umap_row}

{provenance_blob}

<h2>Spatial overview — leiden clusters &amp; cell type</h2>
<p class="meta">
  From <code>{sample_id}_proseg_purified.h5ad</code>:
  <code>obsm['spatial']</code> + <code>obs['{leiden_col}']</code> /
  <code>obs['first_type']</code>. HE is not available at this stage of
  the pipeline, so these stand alone rather than overlaying anything.
</p>
<div class="row broad">
  <div class="plot"><img src="{img_broad_leiden}" alt="Spatial leiden clusters">
    <div class="caption">Leiden clusters ({leiden_col}, {n_leiden} clusters)</div></div>
  <div class="plot"><img src="{img_broad_celltype}" alt="Spatial cell type">
    <div class="caption">SPLIT-inferred cell type ({n_celltype} types)</div></div>
</div>

<h2>QC metrics in space</h2>
<p class="meta">
  Same information as the count histograms above, plotted spatially
  instead of as a distribution (cell type has its own broad panel
  above, so it isn't repeated here). Color scale only is clipped to
  each panel's own 1st&ndash;99th percentile (independently per
  metric) &mdash; every cell is still plotted, points beyond the clip
  range simply saturate at the end color.
</p>
<div class="row small">
  <div class="plot"><img src="{img_small_raw}" alt="Spatial proseg_raw total_counts, color clipped 1st-99th percentile">
    <div class="caption">proseg_raw: total_counts</div></div>
  <div class="plot"><img src="{img_small_xenium}" alt="Spatial xenium_ranger total_counts, color clipped 1st-99th percentile">
    <div class="caption">xenium_ranger: total_counts</div></div>
  <div class="plot"><img src="{img_small_purified}" alt="Spatial proseg_purified nCount_Proseg, color clipped 1st-99th percentile">
    <div class="caption">proseg_purified: nCount_Proseg</div></div>
</div>
</body>
</html>
"""


def build_report(
    *,
    sample_id: str,
    run_id: str,
    xenium_h5ad: Path,
    raw_h5ad: Path,
    purified_h5ad: Path,
    plots_dir: Path,
    existing_html_path: Path,
    out_path: Path,
) -> None:
    import anndata as ad

    from rctd_split._internal.palette import celltype_color_map

    plt = _configure_matplotlib()

    print(f"[qc_report_spatial_v2] reading {purified_h5ad}")
    purified = ad.read_h5ad(purified_h5ad, backed="r")
    xy_pur = purified.obsm["spatial"]

    leiden_cols = sorted(
        c for c in purified.obs.columns if c.startswith("leiden_")
    )
    if not leiden_cols:
        raise SystemExit(
            "[qc_report_spatial_v2] proseg_purified.h5ad has no leiden_* "
            "obs column."
        )
    leiden_col = leiden_cols[0]
    leiden_vals = purified.obs[leiden_col].astype(str).to_numpy()

    if "first_type" not in purified.obs.columns:
        raise SystemExit(
            "[qc_report_spatial_v2] proseg_purified.h5ad missing "
            "obs['first_type']."
        )
    celltype_vals = purified.obs["first_type"].astype(str).to_numpy()

    if "nCount_Proseg" not in purified.obs.columns:
        raise SystemExit(
            "[qc_report_spatial_v2] proseg_purified.h5ad missing "
            "obs['nCount_Proseg']."
        )
    ncount_pur = purified.obs["nCount_Proseg"].to_numpy()

    celltype_cmap = celltype_color_map(celltype_vals)

    leiden_uniq = sorted(set(leiden_vals), key=lambda v: (len(v), v))
    leiden_tab = plt.get_cmap("tab10" if len(leiden_uniq) <= 10 else "tab20")
    leiden_cmap = {
        lvl: leiden_tab(i / max(len(leiden_uniq) - 1, 1))
        for i, lvl in enumerate(leiden_uniq)
    }

    print("[qc_report_spatial_v2] rendering broad leiden panel")
    img_broad_leiden = _broad_spatial_categorical(
        plt, xy_pur, leiden_vals,
        color_map=leiden_cmap,
        title=f"{sample_id}: leiden clusters ({leiden_col})",
        legend_title="leiden",
    )
    print("[qc_report_spatial_v2] rendering broad cell-type panel")
    img_broad_celltype = _broad_spatial_categorical(
        plt, xy_pur, celltype_vals,
        color_map=celltype_cmap,
        title=f"{sample_id}: SPLIT-inferred cell type",
        legend_title="cell type",
    )

    print("[qc_report_spatial_v2] rendering small purified nCount_Proseg panel (clipped 1st-99th %ile)")
    img_small_purified = _small_spatial_continuous(
        plt, xy_pur, ncount_pur, title="purified: nCount_Proseg",
    )

    print(f"[qc_report_spatial_v2] reading {raw_h5ad}")
    raw = ad.read_h5ad(raw_h5ad, backed="r")
    if "total_counts" not in raw.obs.columns:
        raise SystemExit(
            "[qc_report_spatial_v2] proseg_raw.h5ad missing "
            "obs['total_counts']."
        )
    print("[qc_report_spatial_v2] rendering small proseg_raw total_counts panel (clipped 1st-99th %ile)")
    img_small_raw = _small_spatial_continuous(
        plt, raw.obsm["spatial"], raw.obs["total_counts"].to_numpy(),
        title="proseg_raw: total_counts",
    )

    print(f"[qc_report_spatial_v2] reading {xenium_h5ad}")
    xenium = ad.read_h5ad(xenium_h5ad, backed="r")
    if "total_counts" not in xenium.obs.columns:
        raise SystemExit(
            "[qc_report_spatial_v2] xenium_ranger.h5ad missing "
            "obs['total_counts']."
        )
    print("[qc_report_spatial_v2] rendering small xenium_ranger total_counts panel (clipped 1st-99th %ile)")
    img_small_xenium = _small_spatial_continuous(
        plt, xenium.obsm["spatial"], xenium.obs["total_counts"].to_numpy(),
        title="xenium_ranger: total_counts",
    )

    print(f"[qc_report_spatial_v2] lifting carried-over sections from {existing_html_path}")
    lifted = _lift_existing_sections(existing_html_path, plots_dir)

    body = _HTML_TEMPLATE.format(
        sample_id=html.escape(sample_id),
        run_id=html.escape(run_id),
        generated_ts=time.strftime("%Y-%m-%dT%H:%M:%S"),
        leiden_col=html.escape(leiden_col),
        n_leiden=len(leiden_uniq),
        n_celltype=len(set(celltype_vals)),
        img_broad_leiden=img_broad_leiden,
        img_broad_celltype=img_broad_celltype,
        img_small_raw=img_small_raw,
        img_small_xenium=img_small_xenium,
        img_small_purified=img_small_purified,
        hist_blob=lifted["hist"],
        umap_row=lifted["umap_row"],
        tabular_blob=lifted["tabular"],
        provenance_blob=lifted["provenance"],
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(body)
    print(f"[qc_report_spatial_v2] wrote {out_path}")


def _main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--run-dir", required=True, type=Path,
                    help="Run directory, e.g. .../<S>_env_config_v1")
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args(argv)

    sample_id = args.sample_id
    run_dir = args.run_dir
    spatial_dir = run_dir / "spatial_adata"
    summary_dir = run_dir / "summary"
    plots_dir = summary_dir / "plots"

    build_report(
        sample_id=sample_id,
        run_id=run_dir.name,
        xenium_h5ad=spatial_dir / f"{sample_id}_xenium_ranger.h5ad",
        raw_h5ad=spatial_dir / f"{sample_id}_proseg_raw.h5ad",
        purified_h5ad=spatial_dir / f"{sample_id}_proseg_purified.h5ad",
        plots_dir=plots_dir,
        existing_html_path=summary_dir / f"{sample_id}_summary_report.html",
        out_path=args.out,
    )


if __name__ == "__main__":
    _main()
