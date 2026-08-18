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
  from ``first_type`` labels onto the tab20 palette. The color
  each celltype receives is a stable function of its name string
  (hash-based mod 20), so the SAME celltype gets the SAME color
  across every sample the pipeline runs — even though each
  sample's ``first_type`` set is different. Sorted alphabetically
  in the returned dict for readable JSON dumps.

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


# tab20 hex codes (RGB tuples converted). Matches matplotlib's
# 'tab20' qualitative cmap so downstream notebooks that render
# celltype legends with `plt.get_cmap('tab20')` see the same
# colors.
_TAB20_HEX: tuple[str, ...] = (
    "#1f77b4", "#aec7e8",
    "#ff7f0e", "#ffbb78",
    "#2ca02c", "#98df8a",
    "#d62728", "#ff9896",
    "#9467bd", "#c5b0d5",
    "#8c564b", "#c49c94",
    "#e377c2", "#f7b6d2",
    "#7f7f7f", "#c7c7c7",
    "#bcbd22", "#dbdb8d",
    "#17becf", "#9edae5",
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
    """Return a color map for celltype labels using a deterministic
    per-name hash → tab20 lookup.

    Consistency guarantee: two pipeline invocations that see the
    same celltype name produce the same hex color. That gives
    cross-sample color consistency without needing a globally
    shared celltype list.

    Empty/missing labels map to gray. Occasional colour collisions
    are possible (two celltypes hashing to the same tab20 slot)
    but rare with typical sample sizes.
    """
    result: dict[str, str] = {}
    for lv in sorted({_normalize(l) for l in levels}):
        if lv == "":
            result[lv] = "#BBBBBB"
            continue
        result[lv] = _TAB20_HEX[_hash_index(lv, len(_TAB20_HEX))]
    return result
