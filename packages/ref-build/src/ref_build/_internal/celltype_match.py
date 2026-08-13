"""Fuzzy celltype label matching.

Ported from user request 2026-07-10 (internal issue review) on
(internal issue review) — the census stage was doing strict-equality
matching (`obs[celltype_col] == X`), which silently dropped composite labels
like `B/Plasma_T/NK_rbc` and `unknown_maybe_Fibroblast`, leading to
`missing_no_donor` decisions for celltypes the data DID contain.

Three match rules — a cell with label `L` matches expected celltype `X` if
ANY of:

  1. **Exact**       — `L == X`. Always active.
  2. **Unknown-maybe** — `L == f"unknown_maybe_{X}"`. Active when
                          `include_unknown_maybe=True`.
  3. **Composite delimited** — `X` appears in `L` as a token surrounded by
                                 `_` / `/` / string boundary. Regex:
                                 ``(?:^|[_/])re.escape(X)(?:$|[_/])``.
                                 Active when `fuzzy=True`.

`re.escape` is important — expected celltypes may themselves contain `/`
(e.g. `T/NK`), which would otherwise be interpreted as a regex metachar.

Backward-compat: `fuzzy=False` + `include_unknown_maybe=False` reduces to
strict equality — the pre-2026-07-10 behavior.
"""
from __future__ import annotations

import re
from collections.abc import Callable, Iterable


def build_matcher(
    celltype: str,
    *,
    fuzzy: bool = True,
    include_unknown_maybe: bool = True,
) -> Callable[[str], bool]:
    """Return a predicate ``label -> bool`` for `celltype` under the given rules.

    See the module docstring for rule semantics. Non-string labels (e.g.
    NaN from pandas' object columns) never match.
    """
    exact = celltype
    unknown_maybe = f"unknown_maybe_{celltype}"
    composite_re = (
        re.compile(rf"(?:^|[_/]){re.escape(celltype)}(?:$|[_/])")
        if fuzzy
        else None
    )

    def _match(label: object) -> bool:
        if not isinstance(label, str):
            return False
        if label == exact:
            return True
        if include_unknown_maybe and label == unknown_maybe:
            return True
        if composite_re is not None and composite_re.search(label):
            return True
        return False

    return _match


def matching_labels(
    celltype: str,
    labels: Iterable[object],
    *,
    fuzzy: bool = True,
    include_unknown_maybe: bool = True,
) -> list[str]:
    """Return the ordered, deduped subset of `labels` matching `celltype`.

    First-seen order is preserved — useful for the census's `matched_labels`
    audit column, which lists actual on-cell labels that were rolled up under
    each expected celltype `X`.
    """
    predicate = build_matcher(
        celltype, fuzzy=fuzzy, include_unknown_maybe=include_unknown_maybe,
    )
    seen: set[str] = set()
    kept: list[str] = []
    for L in labels:
        if not isinstance(L, str):
            continue
        if L in seen:
            continue
        seen.add(L)
        if predicate(L):
            kept.append(L)
    return kept
