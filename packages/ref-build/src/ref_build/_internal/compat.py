"""Small compatibility shims used across stages.

Kept minimal on purpose; if something grows past ~30 lines it probably
wants its own module.
"""
from __future__ import annotations

from pathlib import Path


def sentinel_exists(path: Path, force_rerun: bool) -> bool:
    """Return True if `path` already exists and `force_rerun` is False.

    Every stage guards its work behind a call to this helper so that
    re-running the pipeline picks up where the last run left off.
    """
    return path.exists() and not force_rerun


_LEGACY_NULL_PATCHED = False


def patch_legacy_null_encoding() -> None:
    """Register a reader for legacy `IOSpec("null", "0.1.0")` scalars.

    Older scanpy/anndata releases wrote scalar ``None`` values (notably
    ``uns/log1p/base`` produced by :func:`scanpy.pp.log1p` when the base
    was left unspecified) as an h5py Dataset with attrs
    ``encoding-type='null'`` / ``encoding-version='0.1.0'``. Current
    anndata (0.11.x) has no registered reader for that spec, so
    :func:`anndata.read_h5ad` raises ``IORegistryError`` when it
    encounters such a field. Reproduced on an internal input
    (``/uns/log1p/base``); see (internal issue review).

    This shim registers a ``None``-returning reader for the missing
    spec on both ``h5py.Dataset`` and ``h5py.Group`` (and their zarr
    counterparts when zarr is available). Idempotent — a second call
    is a no-op.
    """
    global _LEGACY_NULL_PATCHED
    if _LEGACY_NULL_PATCHED:
        return

    from anndata._io.specs.registry import _REGISTRY, IOSpec

    def _read_legacy_null(elem, *, _reader):
        return None

    spec = IOSpec("null", "0.1.0")
    src_types: list[type] = []
    try:
        import h5py
        src_types.extend([h5py.Dataset, h5py.Group])
    except ImportError:
        pass
    try:
        import zarr
        src_types.extend([zarr.Array, zarr.Group])
    except ImportError:
        pass

    for src_type in src_types:
        _REGISTRY.register_read(src_type, spec)(_read_legacy_null)

    _LEGACY_NULL_PATCHED = True
