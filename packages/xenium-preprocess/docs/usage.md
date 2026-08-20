# Usage

## Command-line reference

```
xenium-preprocess run [OPTIONS]
```

The three required-at-runtime flags:

| Flag | Meaning |
|---|---|
| `--sample-id`   | Sample identifier (e.g. `SAMPLE1`). Used to name the output subdirectory. |
| `--proseg-dir`  | Directory holding the proseg count matrix + cell metadata. |
| `--output-root` | Root output directory. Per-sample results land at `<output_root>/<sample_id>/`. |

Every other flag has a default in `config/default.yaml`. Override precedence: CLI flag > user YAML (`--config user.yaml`) > default YAML.

### Stage selection + resume

Default stages (2026-08-11 onward): `proseg_to_anndata`, `qc_filter`, `xenium_ranger_to_anndata`, `enrich_xenium_id`. The legacy `preprocess`, `split_prep`, and `rctd_prep` stages remain in `VALID_STAGES` but are no longer default — opt in explicitly.

```bash
xenium-preprocess run --stages proseg_to_anndata qc_filter    # partial pipeline
xenium-preprocess run --stages split_prep rctd_prep           # legacy 10x-bundle → RDS
xenium-preprocess run --force-rerun                            # nuke sentinels + re-run
```

Each stage has a sentinel file at its expected output path. Re-running with the same `--sample-id` and `--output-root` picks up where the last run left off. `--force-rerun` overrides.

### proseg_to_anndata flags

By default the pipeline globs the proseg directory for the count matrix + cell metadata files. Override with:

- `--count-matrix /path/to/expected-counts.parquet` — explicit path to the expected-counts matrix (skips the glob).
- `--maxpost-matrix /path/to/maxpost_counts.parquet` — explicit path to the maxpost integer matrix (skips the glob). Feeds `adata.layers['maxpost_counts']`.
- `--cell-metadata /path/to/cell_metadata.parquet` — same.
- `--x-source {maxpost_counts,expected_counts}` — which layer `.X` mirrors. Default `maxpost_counts` (integer MAP counts — the 2026-07-23 default flip; QC and single-pass preprocess then read integer counts). Pass `expected_counts` to opt back into the pre-2026-07-23 fractional layout. Both matrices are stored as `.layers[...]` regardless.
- `--proseg-run-script /path/to/proseg_run_<sample>.sh` — optional audit trail. When set, the `proseg_to_anndata` stage copies the referenced shell script verbatim into `<run-dir>/spatial_adata/provenance/<script-name>` so a reader of the run folder can trace the raw h5ad back to the exact proseg invocation (CLI flags, model params, seed) that produced its counts. Not the same as `summary/<sample>_summary_report.html`'s "Provenance" section, which captures rctd-split's own invocation + config + package versions. Omitted when the flag is unset — no empty `provenance/` directory is created.

Set `proseg_to_anndata.maxpost_matrix_glob: null` in your user YAML to disable maxpost loading entirely (reverts to the pre-2026-07-10 single-matrix behaviour). With the default `x_source: maxpost_counts` in that setup, `.X` silently falls back to `expected_counts` with a WARN log — no need to also flip `x_source`.

Update the `centroid_x_col` / `centroid_y_col` / `cell_id_col` keys under `proseg_to_anndata` in your user YAML if the proseg version you're running writes different column names.

### enrich_xenium_id flags

Runs after `xenium_ranger_to_anndata` and adds three NN-mapping columns to `<sample>_xenium_ranger.h5ad`. Direction (reversed 2026-08-11, (internal issue review): FOR EACH xenium cell → find nearest proseg cell (by centroid distance) → write the proseg cell id onto the xenium h5ad's obs.

The stage takes no CLI flags — both source h5ads (`<sample>_proseg_raw.h5ad` + `<sample>_xenium_ranger.h5ad`) are located from the pipeline's canonical `spatial_adata/` layout, so nothing extra needs threading. Tune column names / distance threshold via `config.enrich_xenium_id` in a user YAML:

```yaml
enrich_xenium_id:
  enabled: true                       # set false to skip
  distance_threshold: null            # cells beyond N units get marked `unassigned:distance_over_threshold`
  nn_id_col: proseg_cell_id_nn        # xenium.obs column to receive the proseg id
  distance_col: proseg_id_nn_distance
  note_col: proseg_id_nn_note
  write_back_h5ad: true               # false → sidecar `_xenium_enriched.h5ad`
```

### preprocess flags (LEGACY, opt-in only)

- `--dual-matrix-mode true` — default. Runs two independent preprocess passes on the expected + maxpost layers. Set false to revert to a single pass on `.X`.
- `--leiden-resolution 0.4` — resolution passed to `leidenalg.find_partition` (matches the reference notebook).
- `--min-counts-cell 10` — floor for `sc.pp.filter_cells`.
- `--n-neighbors 15` — neighbours for KNN + UMAP.
- `--min-prop 0.001` — clip floor for the clipped-log normalization.
- `--positive-x false` — set true to shift `X_clipped` positive.
- `--random-state 0` — seed for KNN + UMAP + Leiden.
- `--global-non-tumor-json /path/to/markers_global_no_tumor_cells_level1.json` — enables celltype inference (runs independently per pass in dual mode).
- `--global-tumor-json /path/to/markers_global_tumor_cells_level1.json` — for the tumor-marker overlay UMAP.
- `--tumor-type Prostate` — key into the tumor-marker JSON.

Without the marker-JSON flags, the `preprocess` stage still writes the unpurified adata (with Leiden clusters + UMAP + PCA populated) plus `umap_leiden{_pass}.png`, and just skips the celltype-inference sub-step.

### split_prep flags (LEGACY, opt-in only)

Dropped from `DEFAULT_STAGES` on 2026-08-11 (user request, (internal issue review). Opt in with `--stages ... split_prep ...`.

- `--split-layer maxpost_counts` — which adata layer to export as the mtx. Default `maxpost_counts` (the integer proseg output — needed by RCTD's `require_int=TRUE`). Overridden by `rctd_prep.source_layer` when the `rctd_prep` stage is in the run.
- `--split-name-suffix _lowCountThreshold` — file-name suffix (matches the `proseg_<sample>_lowCountThreshold_*` convention).

### rctd_prep flags (LEGACY, opt-in only)

Dropped from `DEFAULT_STAGES` on 2026-08-11 together with `split_prep` (which it consumes). Opt in with `--stages ... split_prep rctd_prep ...`.

- `--rscript-bin Rscript` — path to the R interpreter. Default `Rscript` (PATH lookup); pass an absolute path to skip PATH.

## User YAML

If you have a lot of overrides, park them in a user YAML and load it in one flag:

```yaml
# my_run.yaml
preprocess:
  leiden_resolution: 0.5
  n_neighbors: 20
  global_non_tumor_json: /data/markers/markers_global_no_tumor_cells_level1.json
  global_tumor_json: /data/markers/markers_global_tumor_cells_level1.json
  tumor_type: Prostate
split_prep:
  name_suffix: _proseg_v2
```

```bash
xenium-preprocess run \
    --sample-id SAMPLE1 \
    --proseg-dir /data/SAMPLE1/proseg \
    --output-root /data/xenium_preprocess_runs \
    --config my_run.yaml
```

CLI flags always win over the user YAML.

## Output layout

See the [Outputs](../README.md#outputs) section of the README. Every run also writes a `config.yaml` at the run root, capturing the exact merged config that the run used — useful when reproducing a specific result months later.
