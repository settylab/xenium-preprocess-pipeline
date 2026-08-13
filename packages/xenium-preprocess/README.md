# xenium-preprocess

**Step 1 of the Xenium spatial-data preprocessing pipeline: proseg output → raw AnnData + xenium-ranger AnnData with a NN-mapped proseg cell_id column joining the two.**

At a glance (default stages):

- **`proseg_to_anndata`** — reads a proseg output directory (BOTH the expected-counts and maxpost-counts matrices + cell metadata) and writes `<sample>_proseg_raw.h5ad` with `.X` = expected-counts, `.layers['expected_counts']` mirroring `.X`, `.layers['maxpost_counts']` for the integer MAP assignments, `.obs` = cell metadata, and `.obsm['spatial']` = centroids.
- **`qc_filter`** — annotates the raw proseg h5ad with `.obs['qc_filtered']` (bool). Strictly additive — no cells removed.
- **`xenium_ranger_to_anndata`** — reads the xenium-ranger bundle (`cell_feature_matrix.h5` + `cells.csv.gz`) and writes `<sample>_xenium_ranger.h5ad` with `.X` = raw counts, `.obsm['spatial']` = centroids, `.obs` carrying the standard xenium QC columns.
- **`enrich_xenium_id`** — for each xenium cell, finds the nearest proseg cell by centroid and adds a `proseg_cell_id_nn` column (plus `proseg_id_nn_distance` + `proseg_id_nn_note`) to `<sample>_xenium_ranger.h5ad`. Strictly additive — no existing xenium obs columns mutated. Reversed 2026-08-11 from the previous "xenium id onto proseg h5ad" direction (internal issue review).

Legacy (opt-in via `--stages ... split_prep rctd_prep`, no longer default):

- **`preprocess`** — QC + PCA/UMAP/Leiden/celltype inference. Deprecated; retained only for reversibility.
- **`split_prep`** — 10x-style `{counts.mtx.gz, features.tsv.gz, barcodes.tsv.gz}` triple + sidecars. Dropped from defaults 2026-08-11 (user request: remove all 10x-bundle output).
- **`rctd_prep`** — R + SpatialExperiment → RCTD test object RDS. Consumed the `split_prep/` triple, so it went out with it.

Every stage is idempotent — sentinel-file resume — and every run snapshots its resolved configuration to disk.

Table of contents: [What it does](#what-it-does) · [Installation](#installation) · [Quickstart](#quickstart) · [Outputs](#outputs) · [Configuration](#configuration) · [HPC (Slurm)](#hpc-usage-slurm) · [Citation](#citation) · [License](#license) · [Acknowledgements](#acknowledgements)

## What it does

Default stage order (`proseg_to_anndata → qc_filter → xenium_ranger_to_anndata → enrich_xenium_id`). Restrict a run to a subset with `--stages`; each stage is guarded by a sentinel file, so re-running with the same `--sample-id` and `--output-root` picks up where the last run left off. Nuke sentinels with `--force-rerun`. Opt into the legacy `preprocess` / `split_prep` / `rctd_prep` stages by naming them explicitly under `--stages`.

### `proseg_to_anndata` — proseg output → raw AnnData

Adapted from the internal proseg-to-anndata reference script, extended 2026-07-10 for dual-matrix support (user request, (internal issue review). Reads BOTH proseg count matrices — expected-counts (continuous) and maxpost-counts (integer) — plus the cell metadata table, auto-detecting `.parquet` vs `.csv` from the extension. The centroids (`centroid_x`, `centroid_y`) are stacked into `adata.obsm['spatial']`; the metadata rows become `adata.obs` (indexed by `cell` → renamed to `cell_d` to sidestep anndata's index-name collision rule); counts are stored as CSR sparse matrices.

Layers:

- `adata.X` = `adata.layers['expected_counts']` — continuous expected transcript counts (unchanged downstream contract for anything that reads `.X`).
- `adata.layers['maxpost_counts']` — integer maximum-a-posteriori assignments (fed to RCTD in stage 4).

The two matrices are validated for exact cells × genes match and identical gene column order before the h5ad is written; a shape or column-order mismatch fails loud. Setting `--maxpost-matrix` or `proseg_to_anndata.maxpost_matrix_glob: null` reverts to a single-matrix load (expected-counts only, no `maxpost_counts` layer).

**Writes:** `<output_root>/<sample_id>/proseg_adata/raw.h5ad`.

### `preprocess` — QC → normalization → clustering → celltype inference

Adapted verbatim from the internal preprocessing reference notebook, using the routines that live in the internal source the internal Xenium data-processing routines. Those routines have been ported into `xenium_preprocess.stages.preprocess` (each with `# Adapted from …` in the docstring) so the pipeline does not depend on `XeniumDataProcessing` being on `sys.path` in the runtime env.

The pipeline stops *before* the manual `label_clusters` fixes in the reference notebook — those are curation and areuser-driven. The output is deliberately named `adata_unpurified.h5ad`.

**Steps** (each parameterised in `config/default.yaml`):

1. `remove_control_probes` — drop genes whose name starts with `Neg` or `Unassigned`. (Shared across both dual-matrix passes.)
2. `sc.pp.calculate_qc_metrics` with `percent_top=(10, 20, 50, 150)` + `sc.pp.filter_cells(min_counts=10)`. (Shared.)
3. `daniel_approach_hvg_pca` — clipped-log normalization (`X /= X.sum(axis=1); log(clip(X, min_prop, 1))`), stored as `X_clipped{suffix}`; PCA on that; KNN via `umap.umap_.nearest_neighbors`; UMAP via `umap.UMAP(precomputed_knn=…)`.
4. `leiden_clustering_via_knn_graph` — build an unweighted igraph from `obsm['knn_indices{suffix}']` and run `leidenalg.find_partition` with `RBConfigurationVertexPartition` at `resolution=0.4`.
5. **Optional** celltype inference (fires when `--global-non-tumor-json` is set): score each cluster against every marker-gene group in the JSON via `sc.tl.score_genes`; top-scoring module per cluster; low-scoring → `unknown_maybe_tumor`; duplicate-winner clusters merged (`B_T`, `Myeloid_T`, …); low-count → `low count`. 1:1 with `generate_dict_for_celltype_annotation_global` in the reference notebook. Runs INDEPENDENTLY per pass — labels may differ between the expected and maxpost passes.

**Dual-matrix mode (default).** Steps 3–5 run twice — once with `source_layer=expected_counts`, `suffix=_expected`, and once with `source_layer=maxpost_counts`, `suffix=_maxpost`. New adata slots after dual-mode preprocess:

- `adata.obsm`: `X_pca_expected`, `X_pca_maxpost`, `X_umap_expected`, `X_umap_maxpost`, `knn_indices_expected`, `knn_indices_maxpost`, `knn_dists_expected`, `knn_dists_maxpost`.
- `adata.layers`: `X_clipped_expected`, `X_clipped_maxpost` (plus the two source layers `expected_counts` and `maxpost_counts` unchanged). The dense float32 `counts_expected` / `counts_maxpost` layers were dropped 2026-07-30 (internal issue review) — they were per-pass copies of the source layers, ~10 GiB overhead with no downstream consumer.
- `adata.obs`: `leiden_0.4_expected`, `leiden_0.4_maxpost`, `global_level1_celltype_expected`, `global_level1_celltype_maxpost` (when celltype JSON given).

Set `preprocess.dual_matrix_mode: false` (or `--dual-matrix-mode false`) to revert to a single pass on `.X` (writing unsuffixed keys). Combined with `proseg_to_anndata.maxpost_matrix_glob: null`, this reverts to the pre-2026-07-10 single-matrix pipeline.

**Writes:** `preprocessed/adata_unpurified.h5ad`, and per pass `preprocessed/plots/umap_leiden{_pass}.png`, `preprocessed/plots/umap_celltype{_pass}.png`, `preprocessed/plots/umap_tumor_markers_<TumorType>{_pass}.png` (in single-mode the plots are named `umap_leiden.png`, `umap_global_level1_celltype.png`, `umap_tumor_markers_<TumorType>.png` for backward compat).

### `split_prep` — 10x-style mtx triple + sidecars (LEGACY, opt-in only)

Dropped from `DEFAULT_STAGES` on 2026-08-11 (user request, (internal issue review): the 10x-bundle is a legacy artifact of the pre-reorg architecture. Kept in `VALID_STAGES` — opt in with `--stages ... split_prep ...` if you need the bundle for a legacy R/Seurat workflow. When enabled, the stage still materialises `adata.layers['maxpost_counts']`, casts to `int32`, transposes to `genes × cells`, and writes:

- `<sample>{suffix}_counts.mtx.gz` — MatrixMarket, gzipped.
- `<sample>{suffix}_features.tsv.gz` — one gene name per line (single-column form; matches Seurat `ReadMtx(feature.column=1)`).
- `<sample>{suffix}_barcodes.tsv.gz` — one cell id per line.
- `<sample>{suffix}_metadata.csv` — `adata.obs`, cell-indexed.
- `<sample>{suffix}_spatial_coords.csv.gz` — `(x, y)` per cell.

The suffix defaults to `_lowCountThreshold` (matches the `proseg_<sample>_lowCountThreshold_*` convention); override with `--split-name-suffix`.

**Writes (only when opted in):** `split_prep/` with the five files above.

### `rctd_prep` — SPLIT triple → RCTD test object (RDS) (LEGACY, opt-in only)

Dropped from `DEFAULT_STAGES` on 2026-08-11 (user request, (internal issue review) together with `split_prep`, which it consumes. Kept in `VALID_STAGES` — opt in with `--stages ... split_prep rctd_prep ...`. When enabled, Python shells out to `Rscript src/xenium_preprocess/r/rctd_prep.R`, threading through the paths + args. The R script:

1. Reads the mtx triple via `Seurat::ReadMtx(feature.column=1, cell.column=1)`.
2. Builds a `SpatialExperiment` from the counts + metadata + spatial coords.
3. Converts to a Seurat object with a `Proseg` assay.
4. Attaches the spatial coordinates as a `DimReduc` keyed `ST_` (columns renamed to `ST_1`, `ST_2`).
5. `saveRDS(seu, "test_object.rds")`.

That RDS is the RCTD "test object" — downstream `spacexr::create.RCTD(spatial_seurat, reference)` calls take it as the query side. Since RCTD's `require_int=TRUE` rejects the continuous expected-counts matrix, `rctd_prep.source_layer` defaults to `maxpost_counts` (threaded into `split_prep.layer` by `pipeline.py`).

The second half of the reference Rmd (building an RCTD `Reference` from a scRNA mtx triple) is out of scope for step 1; it will land in a later step.

**Writes:** `rctd_prep/test_object.rds`.

## Installation

`xenium-preprocess` targets Python 3.10+ (validated on 3.11). The R stage depends on `Seurat`, `Matrix`, `SpatialExperiment`, `readr`, and `spacexr` — the `fhR/4.4.1-foss-2023b` R module (Lmod) carries all of these.

```bash
# 1. Install micromamba (skip if you already have it)
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# 2. Create the env
micromamba env create -f environment.yml
micromamba activate xeniumPreprocess

# 3. Editable install
pip install -e /path/to/xenium-preprocess

# 4. Verify
xenium-preprocess --help
xenium-preprocess run --help
```

## Quickstart

```bash
micromamba activate xeniumPreprocess
ml fhR/4.4.1-foss-2023b   # for the R stage; skip if you already have Rscript with Seurat + spacexr

xenium-preprocess run \
    --sample-id       SAMPLE1 \
    --proseg-dir      /data/SAMPLE1/proseg \
    --output-root     /data/xenium_preprocess_runs
```

To run the celltype-inference sub-step of the `preprocess` stage, pass the two marker-JSON paths and the tumor type:

```bash
xenium-preprocess run \
    --sample-id       SAMPLE1 \
    --proseg-dir      /data/SAMPLE1/proseg \
    --output-root     /data/xenium_preprocess_runs \
    --global-non-tumor-json /data/markers/markers_global_no_tumor_cells_level1.json \
    --global-tumor-json     /data/markers/markers_global_tumor_cells_level1.json \
    --tumor-type            Prostate
```

Skip a stage set with `--stages`:

```bash
# Regenerate only the RCTD test object (assumes SPLIT triple exists).
xenium-preprocess run --sample-id SAMPLE1 --proseg-dir /… --output-root /… \
    --stages rctd_prep
```

## Outputs

Under `<output_root>/<sample_id>/<sample_id>_<run_id>/` (default stages):

```
spatial_adata/<sample>_proseg_raw.h5ad         — proseg_to_anndata (+ expected_counts, maxpost_counts layers);
                                                 qc_filter adds .obs['qc_filtered']
spatial_adata/<sample>_xenium_ranger.h5ad      — xenium_ranger_to_anndata; enrich_xenium_id then adds
                                                 .obs['proseg_cell_id_nn' + 'proseg_id_nn_distance'
                                                     + 'proseg_id_nn_note']
resolved_config.yaml                           — snapshot of the resolved config (merged across steps)
logs/step1.log                                 — the step-1 log
```

Opting into the legacy `preprocess` / `split_prep` / `rctd_prep` stages via `--stages` additionally writes (respectively): `legacy_preprocess/adata_unpurified.h5ad` + UMAP plots; `split_prep/<sample>{suffix}_*` (mtx triple + metadata + spatial_coords); `rctd/<sample>_test_object.rds`.

## Configuration

Every knob lives in `config/default.yaml`. Override precedence: CLI flag > user YAML (`--config …`) > default YAML.

The four required-at-runtime keys (`sample_id`, `proseg_dir`, `output_root`) are `null` in the default YAML — the pipeline refuses to run until they're supplied. Everything else has a sensible default matched to the reference notebooks.

Full CLI reference: `xenium-preprocess run --help`. Semantics of each knob: `docs/usage.md`.

## HPC usage (Slurm)

```bash
sbatch scripts/submit.slurm.sh SAMPLE1 /data/SAMPLE1/proseg \
    --output-root /data/xenium_preprocess_runs
```

The submit wrapper carries the Slurm-tee/pipefail deadlock fix (see `scripts/submit.slurm.sh:70` for the comment block). Under sbatch it uses a plain `exec >> "$LOG_FILE" 2>&1`; interactive `bash submit.slurm.sh` still gets the tee for live console echo. Skipping this split re-introduces the silent-freeze deadlock the H&E pipeline hit on 2026-07.

## Citation

If you use `xenium-preprocess` in published work, cite it via the metadata in [`CITATION.cff`](CITATION.cff). Also cite the upstream tools this pipeline glues together:

- **proseg** — Jones et al., "Cell segmentation of Xenium in-situ transcriptomics data with proseg."
- **scanpy** — Wolf, F.A. *et al.* Genome Biology 19, 15 (2018).
- **leidenalg** — Traag, V.A., Waltman, L. & van Eck, N.J. Sci. Rep. 9, 5233 (2019).
- **UMAP** — McInnes, L. & Healy, J. arXiv:1802.03426 (2018).
- **Seurat** — Hao *et al.* Cell 184, 3573–3587.e29 (2021).
- **RCTD / spacexr** — Cable *et al.* Nat. Biotechnol. 40, 517–526 (2022).

## License

MIT — see [`LICENSE`](LICENSE). See the "Third-party attribution" block in `LICENSE` for the internal SPLIT/Proseg workflow this package is a port of.

## Acknowledgements

Adapted from an internal SPLIT/Proseg workflow. See `docs/methods.md` for a per-stage attribution to the source scripts and notebooks.
