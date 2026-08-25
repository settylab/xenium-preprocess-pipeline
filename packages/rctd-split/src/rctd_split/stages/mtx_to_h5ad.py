"""mtx_to_h5ad: 10X-style mtx bundles → AnnData .h5ad files.

Pure Python (no R). Reads each of the two intermediate mtx bundles
written by ``export_mtx`` and builds an AnnData:

    rows = cells (barcodes.tsv[.gz])
    cols = genes (features.tsv[.gz])
    X    = counts.mtx[.gz], TRANSPOSED from genes × cells to cells × genes.
    obs  = metadata.csv (cell-indexed), reindexed to match barcodes.
    obsm['spatial'] = (n_obs, 2) float array from spatial_coords.csv[.gz].

Writes:

  * ``<S>_proseg_purified.h5ad`` under ``spatial_adata/`` — the
    persisted purified adata.
  * ``<S>_unpurified.h5ad`` under ``intermediate/adata/`` — the
    intermediate unpurified adata, consumed by ``filter_status`` and
    ``writeback_to_raw``. Dropped from persisted final outputs.
"""
from __future__ import annotations

from pathlib import Path

from rctd_split._internal.compat import sentinel_exists
from rctd_split._internal.layout import (
    atomic_write_h5ad,
    intermediate_path,
    spatial_adata_path,
)
from rctd_split._internal.logging import log
from rctd_split._internal.split_native import filter_purified_obs


def _pick_ext(mtx_dir: Path, sample_id: str, variant: str) -> str:
    """Return `.gz` if the counts .mtx is gzipped in this bundle, else ``."""
    gz = mtx_dir / f"{sample_id}_{variant}_counts.mtx.gz"
    plain = mtx_dir / f"{sample_id}_{variant}_counts.mtx"
    if gz.exists():
        return ".gz"
    if plain.exists():
        return ""
    raise SystemExit(
        f"[mtx_to_h5ad] no counts.mtx found in {mtx_dir} for variant {variant!r}"
    )


def _build_adata(
    sample_id: str,
    mtx_dir: Path,
    variant: str,
    unpurified_mtx_dir: Path | None = None,
) -> "anndata.AnnData":
    """Assemble one AnnData from a mtx bundle.

    When ``variant == "purified"`` and ``unpurified_mtx_dir`` is given, the
    corresponding unpurified counts (subset to the purified cells) are attached
    as ``.layers['unpurified_counts']``.
    """
    import anndata as ad  # noqa: F401 — imported so the type hint above resolves
    import pandas as pd
    import scanpy as sc

    ext = _pick_ext(mtx_dir, sample_id, variant)
    stem = mtx_dir / f"{sample_id}_{variant}"

    mtx_path      = Path(f"{stem}_counts.mtx{ext}")
    feats_path    = Path(f"{stem}_features.tsv{ext}")
    barcodes_path = Path(f"{stem}_barcodes.tsv{ext}")
    meta_path     = Path(f"{stem}_metadata.csv")
    coords_path   = Path(f"{stem}_spatial_coords.csv{ext}")

    log(f"[mtx_to_h5ad]   reading {mtx_path}")
    adata = sc.read_mtx(mtx_path).T

    features = pd.read_csv(feats_path, sep="\t", header=None, dtype=str)
    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None, dtype=str)
    if features.shape[1] != 1:
        raise SystemExit(
            f"[mtx_to_h5ad] features file has {features.shape[1]} columns; "
            f"expected single-column form: {feats_path}"
        )
    if barcodes.shape[1] != 1:
        raise SystemExit(
            f"[mtx_to_h5ad] barcodes file has {barcodes.shape[1]} columns; "
            f"expected single-column form: {barcodes_path}"
        )

    gene_names = features.iloc[:, 0].astype(str).tolist()
    cell_ids   = barcodes.iloc[:, 0].astype(str).tolist()

    if adata.shape != (len(cell_ids), len(gene_names)):
        raise SystemExit(
            f"[mtx_to_h5ad] mtx shape {adata.shape} does not match "
            f"(cells={len(cell_ids)}, genes={len(gene_names)})"
        )

    adata.obs_names = cell_ids
    adata.var_names = gene_names

    log(f"[mtx_to_h5ad]   reading {meta_path}")
    meta = pd.read_csv(meta_path, index_col=0)
    meta.index = meta.index.astype(str)
    meta = meta.reindex(cell_ids)
    adata.obs = meta

    # Object columns with NaN + str/bool mix trip h5py's vlen-string writer.
    # Pandas 3.x reads CSV string columns as `pd.StringDtype` (arrow-backed)
    # by default (`future.infer_string=True`), NOT object dtype — the naive
    # `dtype == object` check silently skips them and leaves NaN sentinels
    # in place, which then reappear as literal `"nan"` strings after the
    # h5ad round-trip. Catch both dtypes.
    for c in adata.obs.columns:
        col_dtype = adata.obs[c].dtype
        if col_dtype == object or isinstance(col_dtype, pd.StringDtype):
            adata.obs[c] = adata.obs[c].fillna("").astype(str)

    if coords_path.exists():
        log(f"[mtx_to_h5ad]   reading {coords_path}")
        coords = pd.read_csv(coords_path, index_col=0)
        coords.index = coords.index.astype(str)
        needed = [c for c in ("x", "y") if c in coords.columns]
        if len(needed) != 2:
            log(f"[mtx_to_h5ad]   WARN: spatial_coords has cols "
                f"{list(coords.columns)} — expected 'x' and 'y'; "
                f"skipping obsm['spatial'].")
        else:
            coords_aligned = coords.reindex(cell_ids)[needed].to_numpy()
            adata.obsm["spatial"] = coords_aligned
            log(f"[mtx_to_h5ad]   attached obsm['spatial'] shape={coords_aligned.shape}")
    else:
        log(f"[mtx_to_h5ad]   WARN: {coords_path} missing; "
            "skipping obsm['spatial'].")

    if variant == "purified" and unpurified_mtx_dir is not None:
        _attach_unpurified_layer(
            adata, sample_id, unpurified_mtx_dir, gene_names,
        )

    layer_keys = [k for k in adata.layers.keys() if k is not None]
    log(f"[mtx_to_h5ad]   built {variant} AnnData: "
        f"{adata.n_obs} cells × {adata.n_vars} genes; "
        f".obs cols: {list(adata.obs.columns)}"
        + (f"; layers: {layer_keys}" if layer_keys else ""))
    return adata


def _attach_unpurified_layer(
    adata: "anndata.AnnData",
    sample_id: str,
    unpurified_mtx_dir: Path,
    gene_names: list[str],
) -> None:
    """Load unpurified counts for ``adata``'s cells and attach as
    ``.layers['unpurified_counts']``."""
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import scipy.sparse as sp

    stem = unpurified_mtx_dir / f"{sample_id}_unpurified"
    gz_mtx = Path(f"{stem}_counts.mtx.gz")
    plain_mtx = Path(f"{stem}_counts.mtx")
    if not gz_mtx.exists() and not plain_mtx.exists():
        log(f"[mtx_to_h5ad]   WARN: no unpurified counts.mtx in "
            f"{unpurified_mtx_dir}; skipping layers['unpurified_counts'].")
        return

    unp_ext = _pick_ext(unpurified_mtx_dir, sample_id, "unpurified")
    unp_mtx_path      = Path(f"{stem}_counts.mtx{unp_ext}")
    unp_feats_path    = Path(f"{stem}_features.tsv{unp_ext}")
    unp_barcodes_path = Path(f"{stem}_barcodes.tsv{unp_ext}")

    log(f"[mtx_to_h5ad]   reading unpurified counts for layer: {unp_mtx_path}")
    unp_adata = sc.read_mtx(unp_mtx_path).T
    unp_features = pd.read_csv(
        unp_feats_path, sep="\t", header=None, dtype=str,
    ).iloc[:, 0].astype(str).tolist()
    unp_barcodes = pd.read_csv(
        unp_barcodes_path, sep="\t", header=None, dtype=str,
    ).iloc[:, 0].astype(str).tolist()

    if unp_adata.shape != (len(unp_barcodes), len(unp_features)):
        raise SystemExit(
            f"[mtx_to_h5ad] unpurified mtx shape {unp_adata.shape} does not "
            f"match (cells={len(unp_barcodes)}, genes={len(unp_features)})"
        )

    unp_gene_to_idx = {g: i for i, g in enumerate(unp_features)}
    shared_genes = [g for g in gene_names if g in unp_gene_to_idx]
    if not shared_genes:
        log(f"[mtx_to_h5ad]   WARN: purified and unpurified bundles share "
            f"no gene names; skipping layers['unpurified_counts'].")
        return
    n_missing = len(gene_names) - len(shared_genes)
    if n_missing:
        log(f"[mtx_to_h5ad]   {n_missing}/{len(gene_names)} purified genes "
            "absent from unpurified bundle; zero-filled in the layer.")

    unp_adata.obs_names = unp_barcodes
    missing_cells = [c for c in adata.obs_names if c not in set(unp_barcodes)]
    if missing_cells:
        raise SystemExit(
            f"[mtx_to_h5ad] {len(missing_cells)} purified cells absent from the "
            f"unpurified bundle (first: {missing_cells[:5]}); cannot attach "
            "layers['unpurified_counts']."
        )
    unp_col_idx = [unp_gene_to_idx[g] for g in shared_genes]
    unp_shared = unp_adata[list(adata.obs_names), unp_col_idx]

    purified_gene_to_col = {g: i for i, g in enumerate(gene_names)}
    target_cols = np.asarray(
        [purified_gene_to_col[g] for g in shared_genes], dtype=int,
    )
    x_unp = unp_shared.X
    if sp.issparse(x_unp):
        x_coo = x_unp.tocoo()
        layer = sp.csr_matrix(
            (x_coo.data, (x_coo.row, target_cols[x_coo.col])),
            shape=(adata.n_obs, adata.n_vars),
            dtype=x_unp.dtype,
        )
    else:
        layer = np.zeros((adata.n_obs, adata.n_vars), dtype=x_unp.dtype)
        layer[:, target_cols] = np.asarray(x_unp)
    adata.layers["unpurified_counts"] = layer
    log(f"[mtx_to_h5ad]   attached layers['unpurified_counts'] shape={layer.shape}")


def run_mtx_to_h5ad(
    sample_id: str,
    run_id: str,
    unpurified_mtx_dir: Path,
    purified_mtx_dir: Path,
    output_root: Path,
    h5ad_compression: str | None,
    force_rerun: bool,
) -> tuple[Path, Path]:
    """Read both mtx bundles and write matching .h5ad files.

    Returns ``(unpurified_h5ad, purified_h5ad)`` where the purified path
    is under ``spatial_adata/<S>_proseg_purified.h5ad`` and the
    unpurified path is under ``intermediate/adata/<S>_unpurified.h5ad``.
    """
    unpurified_h5ad = intermediate_path(
        output_root, sample_id, run_id, "unpurified_h5ad",
    )
    purified_h5ad = spatial_adata_path(
        output_root, sample_id, run_id, "proseg_purified",
    )
    unpurified_h5ad.parent.mkdir(parents=True, exist_ok=True)
    purified_h5ad.parent.mkdir(parents=True, exist_ok=True)

    # Sentinel = purified.h5ad (second write).
    if sentinel_exists(purified_h5ad, force_rerun) and unpurified_h5ad.exists():
        log(f"[mtx_to_h5ad] sentinels exist: {unpurified_h5ad}, {purified_h5ad} "
            "— skipping (pass --force-rerun to re-run).")
        return unpurified_h5ad, purified_h5ad

    log("[mtx_to_h5ad] --- unpurified variant ---")
    adata_unp = _build_adata(sample_id, unpurified_mtx_dir, "unpurified")
    write_kwargs = {"compression": h5ad_compression} if h5ad_compression else {}
    atomic_write_h5ad(adata_unp, unpurified_h5ad, **write_kwargs)
    log(f"[mtx_to_h5ad] wrote {unpurified_h5ad}")

    log("[mtx_to_h5ad] --- purified variant ---")
    adata_pur = _build_adata(
        sample_id, purified_mtx_dir, "purified",
        unpurified_mtx_dir=unpurified_mtx_dir,
    )
    # SPLIT-native obs cleanup: preserve_meta_from_unpurified (R side)
    # copies proseg_raw-side metadata (e.g. `total_counts`) onto the
    # purified Seurat, which is then persisted into purified.obs by
    # `_build_adata`. the user wants purified.obs to carry ONLY SPLIT
    # output. Filter here at the h5ad boundary so the persisted file
    # never contains foreign columns.
    dropped = filter_purified_obs(adata_pur)
    if dropped:
        log(f"[mtx_to_h5ad]   dropped {len(dropped)} non-SPLIT-native "
            f"obs cols from purified: {dropped}")
    atomic_write_h5ad(adata_pur, purified_h5ad, **write_kwargs)
    log(f"[mtx_to_h5ad] wrote {purified_h5ad}")

    return unpurified_h5ad, purified_h5ad
