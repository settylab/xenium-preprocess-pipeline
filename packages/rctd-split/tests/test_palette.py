"""Unit tests for `rctd_split._internal.palette`.

Cross-sample consistency invariants:

  * `purification_status_color_map` returns a canonical color for
    every known level regardless of input ordering / capitalization
    of neighbors.
  * `celltype_color_map` returns the SAME hex color for the SAME
    celltype name across two independent invocations — this is the
    load-bearing property Tracy asked for (color continuity across
    per-sample runs on
    settylab/TracyY123-nexus#26 comment 5322126093, item 3).
"""
from __future__ import annotations

from rctd_split._internal.palette import (
    celltype_color_map,
    purification_status_color_map,
)


def _is_hex(s: str) -> bool:
    return isinstance(s, str) and s.startswith("#") and len(s) == 7


class TestPurificationStatusColorMap:
    def test_known_singlet_reject_canonical_colors(self):
        m = purification_status_color_map(
            ["singlet", "reject", "doublet_certain"]
        )
        # traffic-light intent: singlet green, reject dark magenta
        assert m["singlet"] == "#117733"
        assert m["reject"] == "#882255"
        # doublet_certain gets the pre-declared canonical color
        assert m["doublet_certain"] == "#CC6677"

    def test_legacy_labels_mapped(self):
        m = purification_status_color_map(["purified", "discarded"])
        assert m["purified"] == "#117733"
        assert m["discarded"] == "#882255"

    def test_unknown_level_gets_fallback_hex(self):
        m = purification_status_color_map(["singlet", "reject", "brand_new"])
        assert "brand_new" in m
        assert _is_hex(m["brand_new"])
        assert m["brand_new"] not in ("#117733", "#882255")

    def test_input_order_does_not_matter(self):
        m1 = purification_status_color_map(["singlet", "reject"])
        m2 = purification_status_color_map(["reject", "singlet"])
        assert m1 == m2


class TestCelltypeColorMap:
    def test_returns_hex_for_each_level(self):
        m = celltype_color_map(["T cell", "B cell", "Macrophage"])
        assert set(m.keys()) == {"T cell", "B cell", "Macrophage"}
        for h in m.values():
            assert _is_hex(h)

    def test_same_name_gets_same_color_across_calls(self):
        """The core cross-sample consistency property."""
        m1 = celltype_color_map(["T cell", "B cell", "NK cell"])
        m2 = celltype_color_map(["T cell", "Fibroblast"])
        assert m1["T cell"] == m2["T cell"]

    def test_deterministic_across_neighbors(self):
        """Adding an unrelated celltype must not shift others."""
        base = celltype_color_map(["Alpha", "Beta"])
        extended = celltype_color_map(["Alpha", "Beta", "Gamma"])
        assert base["Alpha"] == extended["Alpha"]
        assert base["Beta"] == extended["Beta"]

    def test_empty_string_is_gray(self):
        m = celltype_color_map(["T cell", ""])
        assert m[""] == "#BBBBBB"
