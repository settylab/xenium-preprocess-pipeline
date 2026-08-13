"""Stage 1: proseg output → raw AnnData.

Adapted from
    the internal proseg-to-anndata reference script
(the original 40-line CLI script). Behaviour preserved for the
expected-counts matrix; dual-matrix support added 2026-07-10 (user request
request, (internal issue review).

- Reads a proseg expected-counts matrix (parquet or csv; auto-detected
  from extension). Rows = cells, columns = genes.
- Optionally ALSO reads a proseg maxpost-counts matrix (same shape) and
  stores it as `adata.layers["maxpost_counts"]`. When the maxpost glob
  is null or the file is missing, the maxpost layer is not populated —
  reverting to the original single-matrix behaviour for backward
  compatibility.
- The expected matrix is also stored as `adata.layers["expected_counts"]`
  to give downstream stages a stable, layer-scoped handle.
- `.X` mirrors ONE of the two layers, chosen by `x_source`
  (`maxpost_counts` — the default, integer MAP counts, user request
  2026-07-23; or `expected_counts` — the pre-2026-07-23 fractional
  default). When `x_source="maxpost_counts"` and the maxpost matrix is
  not loaded, `.X` silently falls back to expected_counts (with a WARN
  log) — preserves back-compat for configs that set
  `maxpost_matrix_glob: null`.
- Reads a proseg cell metadata table (parquet or csv). Must carry
  `centroid_x`, `centroid_y`, and `cell` columns (names configurable).
- Stacks the centroids into `adata.obsm["spatial"]`.
- Uses the cell metadata as `adata.obs`, indexed by `cell` (renamed
  to `cell_d` to match the original script — Anndata's index name
  cannot collide with a regular column).
- Converts each count matrix to a CSR sparse matrix (float — the
  proseg-to-anndata.py original doesn't dtype-cast on this side; the
  int cast happens later, in split_prep, when the exported mtx needs
  integers).
- Validates that the two matrices share cells + genes (exact
  dimensions AND column order) when both are loaded. Fails loud on
  mismatch — a silent gene-order shift here would silently corrupt
  every downstream layer-scoped result.
- Writes `<output_root>/<sample_id>/<sample_id>_<run_id>/spatial_adata/<sample_id>_proseg_raw.h5ad`.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

from xenium_preprocess._internal.compat import sentinel_exists
from xenium_preprocess._internal.layout import (
    atomic_write_h5ad,
    spatial_adata_dir,
    spatial_adata_path,
)
from xenium_preprocess._internal.logging import log


def _read_by_ext(path: Path) -> pd.DataFrame:
    """Read a proseg tabular output; auto-select parquet vs csv."""
    ext = os.path.splitext(str(path))[1].lower()
    if ext == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path)


def _resolve_input(
    label: str,
    explicit: Path | None,
    proseg_dir: Path,
    glob: str,
) -> Path:
    """CLI-explicit path wins; otherwise glob under proseg_dir."""
    if explicit is not None:
        if not explicit.exists():
            raise SystemExit(f"{label}: explicit path does not exist: {explicit}")
        return explicit
    matches = sorted(proseg_dir.glob(glob))
    if not matches:
        raise SystemExit(
            f"{label}: no files matching {glob!r} under {proseg_dir}. "
            f"Pass an explicit path via CLI (e.g. --{label.replace('_', '-')}) "
            f"or update config.proseg_to_anndata.{label}_glob."
        )
    if len(matches) > 1:
        log(f"[proseg_to_anndata] WARN: multiple {label} matches; "
            f"using the first ({matches[0]}). All matches: {[str(m) for m in matches]}")
    return matches[0]


def _resolve_optional(
    label: str,
    explicit: Path | None,
    proseg_dir: Path,
    glob: str | None,
) -> Path | None:
    """Like _resolve_input but returns None when the glob is null or
    matches nothing (with a WARN in the missing case). Used for the
    optional maxpost matrix.
    """
    if explicit is not None:
        if not explicit.exists():
            raise SystemExit(f"{label}: explicit path does not exist: {explicit}")
        return explicit
    if glob is None:
        log(f"[proseg_to_anndata] {label}: glob is null — not loaded.")
        return None
    matches = sorted(proseg_dir.glob(glob))
    if not matches:
        log(f"[proseg_to_anndata] WARN: {label}: no files matching {glob!r} "
            f"under {proseg_dir} — not loaded. Downstream stages that require "
            f"this layer will fall back to `.X` (single-matrix mode).")
        return None
    if len(matches) > 1:
        log(f"[proseg_to_anndata] WARN: multiple {label} matches; "
            f"using the first ({matches[0]}). All matches: {[str(m) for m in matches]}")
    return matches[0]


VALID_X_SOURCES = ("maxpost_counts", "expected_counts")


def run_proseg_to_anndata(
    sample_id: str,
    run_id: str,
    proseg_dir: Path,
    output_root: Path,
    count_matrix_path: Path | None,
    cell_metadata_path: Path | None,
    count_matrix_glob: str,
    cell_metadata_glob: str,
    centroid_x_col: str,
    centroid_y_col: str,
    cell_id_col: str,
    force_rerun: bool,
    maxpost_matrix_path: Path | None = None,
    maxpost_matrix_glob: str | None = None,
    proseg_run_script: Path | None = None,
    x_source: str = "maxpost_counts",
) -> Path:
    """Read the proseg outputs and write a raw h5ad.

    Returns the h5ad path.
    """
    # anndata + scipy imported inside the run function so a bare
    # `python -c "import xenium_preprocess.stages.proseg_to_anndata"`
    # succeeds in an env without them.
    import anndata
    import shutil
    from scipy.sparse import csr_matrix

    if x_source not in VALID_X_SOURCES:
        raise SystemExit(
            f"[proseg_to_anndata] x_source={x_source!r} is not valid. "
            f"Valid values: {VALID_X_SOURCES}."
        )

    out_h5ad = spatial_adata_path(output_root, sample_id, run_id, "proseg_raw")
    out_dir = out_h5ad.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    if sentinel_exists(out_h5ad, force_rerun):
        log(f"[proseg_to_anndata] sentinel exists: {out_h5ad} — skipping "
            f"(pass --force-rerun to re-run).")
        return out_h5ad

    # Optional provenance: copy the proseg-run shell script into the
    # output tree so the origin of the proseg outputs is audit-trailed.
    if proseg_run_script is not None:
        script_src = Path(proseg_run_script)
        if not script_src.exists():
            log(f"[proseg_to_anndata] WARN: proseg_run_script does not "
                f"exist: {script_src} — skipping provenance copy.")
        else:
            prov_dir = spatial_adata_dir(output_root, sample_id, run_id) / "provenance"
            prov_dir.mkdir(parents=True, exist_ok=True)
            prov_dst = prov_dir / script_src.name
            shutil.copy2(script_src, prov_dst)
            log(f"[proseg_to_anndata] provenance: copied {script_src} -> "
                f"{prov_dst}")

    cm_path = _resolve_input(
        "count_matrix", count_matrix_path, proseg_dir, count_matrix_glob
    )
    md_path = _resolve_input(
        "cell_metadata", cell_metadata_path, proseg_dir, cell_metadata_glob
    )
    mp_path = _resolve_optional(
        "maxpost_matrix", maxpost_matrix_path, proseg_dir, maxpost_matrix_glob
    )
    log(f"[proseg_to_anndata] count_matrix:   {cm_path}")
    log(f"[proseg_to_anndata] cell_metadata:  {md_path}")
    log(f"[proseg_to_anndata] maxpost_matrix: {mp_path}")

    count_matrix = _read_by_ext(cm_path)
    cell_metadata = _read_by_ext(md_path)

    for col in (centroid_x_col, centroid_y_col, cell_id_col):
        if col not in cell_metadata.columns:
            raise SystemExit(
                f"cell metadata missing required column {col!r}. "
                f"Present columns: {list(cell_metadata.columns)}. "
                f"Override the column name via config.proseg_to_anndata."
            )

    obsm = {
        "spatial": np.stack(
            [
                cell_metadata[centroid_x_col].to_numpy(),
                cell_metadata[centroid_y_col].to_numpy(),
            ],
            axis=1,
        )
    }
    var = pd.DataFrame(index=count_matrix.columns)
    # Match the original: use `cell` as obs and reindex to a named
    # Index("cell_d") so anndata doesn't complain about the index name
    # colliding with the raw "cell" column.
    obs = cell_metadata.set_index(
        pd.Index(cell_metadata[cell_id_col], str, True, "cell_d")
    )
    expected_csr = csr_matrix(count_matrix.to_numpy())

    layers: dict = {"expected_counts": expected_csr}

    # Load + validate the maxpost matrix if a path was resolved.
    if mp_path is not None:
        maxpost_matrix = _read_by_ext(mp_path)
        # Dimensional match: same n_cells × n_genes.
        if maxpost_matrix.shape != count_matrix.shape:
            raise SystemExit(
                f"[proseg_to_anndata] maxpost matrix shape mismatch: "
                f"expected={count_matrix.shape}, maxpost={maxpost_matrix.shape}. "
                f"Refusing to build the h5ad — the two proseg outputs are "
                f"supposed to describe the same cells × genes."
            )
        # Column-order match: both matrices must carry the same gene
        # names in the same order. A silent column reshuffle would mean
        # `.X[:, i]` and `layers['maxpost_counts'][:, i]` describe
        # different genes.
        if list(maxpost_matrix.columns) != list(count_matrix.columns):
            raise SystemExit(
                f"[proseg_to_anndata] maxpost matrix gene order does not "
                f"match expected matrix. First 5 expected: "
                f"{list(count_matrix.columns[:5])}; first 5 maxpost: "
                f"{list(maxpost_matrix.columns[:5])}."
            )
        # Row-order match: proseg writes both matrices with cells in
        # the same order (matching `cell-metadata.csv`). If the two
        # matrices carry an identifying column, cross-check it — but
        # since proseg count matrices are gene-column-only (no
        # explicit cell-id column), an index-order check by row
        # position is all we can do here. Log a note.
        log(f"[proseg_to_anndata] maxpost matrix shape + gene order match "
            f"expected matrix (both {maxpost_matrix.shape}).")
        layers["maxpost_counts"] = csr_matrix(maxpost_matrix.to_numpy())

    # Pick which layer `.X` mirrors. `maxpost_counts` (the 2026-07-23
    # default) needs the maxpost layer to be loaded; if it isn't (e.g.
    # a pre-dual-mode user YAML with `maxpost_matrix_glob: null`),
    # silently fall back to expected with a WARN log — keeps old configs
    # from breaking. `expected_counts` always works.
    if x_source == "maxpost_counts":
        if "maxpost_counts" in layers:
            x_csr = layers["maxpost_counts"]
            log("[proseg_to_anndata] `.X` <- layers['maxpost_counts'] "
                "(integer MAP counts).")
        else:
            log("[proseg_to_anndata] WARN: x_source='maxpost_counts' but the "
                "maxpost layer was not loaded (maxpost_matrix_glob null or the "
                "file was missing). Falling back to `.X` <- expected_counts. "
                "Set x_source='expected_counts' explicitly to silence this "
                "warning if you intended single-matrix mode.")
            x_csr = expected_csr
    else:  # "expected_counts" (validated at top of function)
        x_csr = expected_csr
        log("[proseg_to_anndata] `.X` <- layers['expected_counts'] "
            "(fractional expected counts, pre-2026-07-23 default).")

    adata = anndata.AnnData(X=x_csr, obs=obs, var=var, obsm=obsm, layers=layers)
    log(f"[proseg_to_anndata] adata: {adata.n_obs} cells × {adata.n_vars} genes; "
        f"layers={list(adata.layers.keys())}")

    atomic_write_h5ad(adata, out_h5ad)
    log(f"[proseg_to_anndata] wrote {out_h5ad}")
    return out_h5ad
