"""Shared NN core tests.

Verifies `nn_map_celltype` on a controlled fixture:
  * k=1, min_dist matches the sklearn-only baseline.
  * unmatched_policy='nearest_label' (default) passes labels through
    even beyond distance_threshold.
  * unmatched_policy='mark_unassigned' fills 'Unassigned' beyond threshold.
"""
from __future__ import annotations

import pytest

pytest.importorskip("sklearn")


def test_nn_basic_k1():
    import numpy as np
    from rctd_split._internal.nn import nn_map_celltype

    ref_xy = np.array([[0.0, 0.0], [10.0, 0.0], [0.0, 10.0]])
    ref_labels = np.array(["A", "B", "C"], dtype=object)
    query_xy = np.array([[0.1, 0.1], [9.9, 0.0], [0.0, 9.9], [5.0, 5.0]])
    out = nn_map_celltype(
        reference_xy=ref_xy,
        reference_labels=ref_labels,
        query_xy=query_xy,
    )
    assert list(out["labels"]) == ["A", "B", "C", "A"]  # (5,5) closer to (0,0)
    assert out["distances"].shape == (4,)
    assert out["reference_idx"].tolist() == [0, 1, 2, 0]


def test_nn_nearest_label_passes_threshold():
    """Default policy: distance_threshold is IGNORED when policy is
    'nearest_label' — labels pass through unchanged."""
    import numpy as np
    from rctd_split._internal.nn import nn_map_celltype

    ref_xy = np.array([[0.0, 0.0]])
    ref_labels = np.array(["A"], dtype=object)
    query_xy = np.array([[1000.0, 1000.0]])
    out = nn_map_celltype(
        reference_xy=ref_xy,
        reference_labels=ref_labels,
        query_xy=query_xy,
        distance_threshold=1.0,
        unmatched_policy="nearest_label",
    )
    assert out["labels"][0] == "A"


def test_nn_mark_unassigned_over_threshold():
    import numpy as np
    from rctd_split._internal.nn import nn_map_celltype

    ref_xy = np.array([[0.0, 0.0]])
    ref_labels = np.array(["A"], dtype=object)
    query_xy = np.array([[1000.0, 1000.0], [0.1, 0.1]])
    out = nn_map_celltype(
        reference_xy=ref_xy,
        reference_labels=ref_labels,
        query_xy=query_xy,
        distance_threshold=1.0,
        unmatched_policy="mark_unassigned",
    )
    assert list(out["labels"]) == ["Unassigned", "A"]


def test_nn_shape_validation():
    import numpy as np
    import pytest as _pt
    from rctd_split._internal.nn import nn_map_celltype

    with _pt.raises(ValueError):
        nn_map_celltype(
            reference_xy=np.zeros((3, 3)),  # bad shape
            reference_labels=np.array(["A", "B", "C"], dtype=object),
            query_xy=np.zeros((2, 2)),
        )


def test_nn_k_ge_1():
    import numpy as np
    import pytest as _pt
    from rctd_split._internal.nn import nn_map_celltype

    with _pt.raises(ValueError):
        nn_map_celltype(
            reference_xy=np.zeros((3, 2)),
            reference_labels=np.array(["A", "B", "C"], dtype=object),
            query_xy=np.zeros((2, 2)),
            k=0,
        )
