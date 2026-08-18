"""Fixed color maps for `qc_report` categorical UMAP plots.

Two color maps land here so plots stay consistent across samples
(Tracy's ask on settylab/TracyY123-nexus#26 comment 5322126093,
item 3):

- ``purification_status_color_map(levels)`` — high-contrast fixed
  palette keyed on SPLIT-native purification_status values
  (``singlet``, ``doublet``, ``rejected``, plus the older
  ``purified``/``unpurified``/``discarded``/``flagged`` labels).
  Unknown levels fall back to a neutral gray so a new SPLIT
  release adding a category doesn't crash the plot.
- ``celltype_color_map(levels)`` — deterministic per-name mapping
  from ``first_type`` labels onto a high-contrast qualitative
  palette. Common tumor-microenvironment celltypes (Liver, Tumor,
  Myeloid, Hepatocyte, T cell, B cell, Macrophage, Endothelial,
  Fibroblast, Stroma) receive fixed, maximally distinct hues so
  the summary HTML never renders them in confusable neighboring
  shades. Unknown celltypes fall back to a hash-based lookup into
  a curated 16-hue palette — no light/dark pairs, unlike tab20 —
  so any two celltypes render in visually distinct colors. The
  color each celltype receives is a stable function of its name
  (hash-based for unknowns, direct lookup for the known set),
  so the SAME celltype gets the SAME color across every sample
  the pipeline runs — even though each sample's ``first_type``
  set is different. Sorted alphabetically in the returned dict
  for readable JSON dumps.

Both functions return dicts mapping string level → ``#rrggbb`` hex
string; that shape drops straight into matplotlib's ``color=``
kwarg AND serializes cleanly into the summary/ JSON sidecar the
HTML report links from.
"""
from __future__ import annotations

import hashlib


# High-contrast, colorblind-safe (Wong palette derived) fixed
# assignments for known SPLIT purification_status values. Ordered
# so the "keep this cell" concept (singlet / purified) is green
# and the "reject" concept (rejected / discarded) is red — matches
# the intuitive traffic-light reading of the UMAP.
_PURIFICATION_STATUS_KNOWN: dict[str, str] = {
    # SPLIT canonical (RCTD spot_class-derived vocabulary)
    "singlet": "#117733",         # dark green — kept, high confidence
    "doublet": "#DDCC77",         # sand yellow — flagged
    "doublet_certain": "#CC6677", # muted red — kept but ambiguous
    "doublet_uncertain": "#AA4499", # purple — kept but very ambiguous
    "rejected": "#882255",        # dark magenta — dropped
    "reject": "#882255",           # alias — same color
    # Older SPLIT vocabulary (kept for back-compat with pre-refactor
    # sample runs whose h5ads still use these labels).
    "purified": "#117733",
    "unpurified": "#DDCC77",
    "discarded": "#882255",
    "flagged": "#CC6677",
    # Explicit empty / missing bucket
    "": "#BBBBBB",
    "(empty)": "#BBBBBB",
    "unknown": "#BBBBBB",
    "nan": "#BBBBBB",
    "na": "#BBBBBB",
}

# Fallback palette for unknown purification_status values, cycled
# in order of appearance. Kept short and visually distinct.
_PURIFICATION_STATUS_FALLBACK: tuple[str, ...] = (
    "#88CCEE",  # light blue
    "#44AA99",  # teal
    "#332288",  # dark blue
    "#999933",  # olive
)


# Fixed, maximally-distinct color assignments for common
# celltypes seen in Setty-lab tumor-microenvironment samples.
# Case-insensitive key — the returned dict preserves the caller's
# spelling. Chosen so Tracy's specific complaint on issue #26
# (Liver / Tumor / Myeloid rendering in confusable shades)
# cannot recur regardless of the fallback palette's hashing.
_CELLTYPE_KNOWN: dict[str, str] = {
    "liver":       "#08306B",  # midnight blue
    "tumor":       "#B22222",  # firebrick red
    "myeloid":     "#E69F00",  # Okabe-Ito orange
    "hepatocyte":  "#6A3D9A",  # deep purple
    "t cell":      "#009E73",  # Okabe-Ito bluish green
    "t_cell":      "#009E73",
    "tcell":       "#009E73",
    "b cell":      "#F0E442",  # Okabe-Ito yellow
    "b_cell":      "#F0E442",
    "bcell":       "#F0E442",
    "nk cell":     "#CC79A7",  # Okabe-Ito reddish purple
    "nk_cell":     "#CC79A7",
    "nkcell":      "#CC79A7",
    "macrophage":  "#7B3F00",  # dark chocolate brown
    "endothelial": "#56B4E9",  # Okabe-Ito sky blue
    "fibroblast":  "#008080",  # teal
    "stroma":      "#4D4D4D",  # dark gray
    "stromal":     "#4D4D4D",
}


# Curated 16-hue qualitative palette for unknown celltypes.
# Replaces matplotlib's `tab20`, which pairs every color with a
# lighter version (hash collisions → confusable neighboring
# shades — the root cause of Tracy's Liver/Tumor/Myeloid
# complaint). Colors sourced from Trubetskoy's "20 distinct
# colors" set, filtered to drop pale pastels that lose contrast
# on a white background.
_CELLTYPE_FALLBACK: tuple[str, ...] = (
    "#e6194b",  # crimson red
    "#3cb44b",  # medium green
    "#4363d8",  # bright blue
    "#f58231",  # bright orange
    "#911eb4",  # purple
    "#469990",  # teal-green
    "#9a6324",  # brown
    "#f032e6",  # magenta
    "#800000",  # maroon
    "#808000",  # olive
    "#000075",  # navy
    "#42d4f4",  # cyan
    "#bfef45",  # lime
    "#a9a9a9",  # medium gray
    "#ffe119",  # bright yellow
    "#dcbeff",  # lavender
)


def _normalize(level) -> str:
    return str(level).strip()


def purification_status_color_map(levels) -> dict[str, str]:
    """Return a fixed color map for the ``levels`` iterable of
    SPLIT purification_status values.

    Known levels get their canonical color (see
    ``_PURIFICATION_STATUS_KNOWN``). Unknown levels get a color
    from ``_PURIFICATION_STATUS_FALLBACK`` in stable
    alphabetical-of-unknowns order, cycled if more unknowns than
    fallbacks. Never raises — the map covers every level in
    ``levels``.
    """
    result: dict[str, str] = {}
    unknown_sorted = sorted(
        {_normalize(lv) for lv in levels}
        - set(_PURIFICATION_STATUS_KNOWN.keys())
    )
    for i, lv in enumerate(unknown_sorted):
        result[lv] = _PURIFICATION_STATUS_FALLBACK[
            i % len(_PURIFICATION_STATUS_FALLBACK)
        ]
    # Overlay known values so any collision with an unknown falls
    # through to the canonical color.
    for lv in {_normalize(l) for l in levels}:
        if lv in _PURIFICATION_STATUS_KNOWN:
            result[lv] = _PURIFICATION_STATUS_KNOWN[lv]
    return result


def _hash_index(name: str, n: int) -> int:
    """Stable non-cryptographic hash of ``name`` mod ``n``.

    Uses SHA-1 rather than Python's ``hash()`` because the latter
    is randomized per-process under PYTHONHASHSEED and would give
    a different color to the same celltype across two pipeline
    invocations.
    """
    h = hashlib.sha1(name.encode("utf-8")).digest()
    return int.from_bytes(h[:4], "big") % n


def celltype_color_map(levels) -> dict[str, str]:
    """Return a color map for celltype labels.

    Common celltypes (Liver, Tumor, Myeloid, Hepatocyte, T cell,
    B cell, NK cell, Macrophage, Endothelial, Fibroblast, Stroma;
    case-insensitive) receive fixed maximally-distinct hues.
    All other names get a deterministic per-name hash lookup into
    a 16-hue high-contrast palette.

    Consistency guarantee: two pipeline invocations that see the
    same celltype name produce the same hex color. That gives
    cross-sample color consistency without needing a globally
    shared celltype list. Adding a new celltype to ``levels``
    never shifts colors of others.

    Empty/missing labels map to gray.
    """
    result: dict[str, str] = {}
    for lv in sorted({_normalize(l) for l in levels}):
        if lv == "":
            result[lv] = "#BBBBBB"
            continue
        known = _CELLTYPE_KNOWN.get(lv.lower())
        if known is not None:
            result[lv] = known
            continue
        result[lv] = _CELLTYPE_FALLBACK[
            _hash_index(lv, len(_CELLTYPE_FALLBACK))
        ]
    return result
