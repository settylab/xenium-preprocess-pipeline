"""Unified `refbuild_<sample>_<tumor>/` layout + legacy-symlink tests
for ref-build (TracyY123-nexus#14 comment 5028672016)."""
from __future__ import annotations

import os
from pathlib import Path


def test_refbuild_dir_name_shape(tmp_path: Path):
    from ref_build._internal.layout import refbuild_dir

    root = tmp_path
    assert refbuild_dir(root, "MH7", "Prostate") == root / "refbuild_MH7_Prostate"
    assert refbuild_dir(root, "MH7", None) == root / "refbuild_MH7"
    # Non-safe characters (spaces, slashes) get slug-safe replacements so
    # the folder name never contains a path separator or shell surprise.
    assert refbuild_dir(root, "MH 7", "Liver Benign") == (
        root / "refbuild_MH_7_Liver_Benign"
    )


def test_new_and_legacy_paths_have_expected_shape(tmp_path: Path):
    from ref_build._internal.layout import legacy_h5ad_path, new_h5ad_path

    root = tmp_path
    assert new_h5ad_path(root, "MH7", "Prostate", "loaded_concat") == (
        root / "refbuild_MH7_Prostate" / "loaded_concat.h5ad"
    )
    assert new_h5ad_path(root, "MH7", "Prostate", "assembled_reference") == (
        root / "refbuild_MH7_Prostate" / "assembled_reference.h5ad"
    )
    assert legacy_h5ad_path(root, "MH7", "loaded_concat") == (
        root / "MH7" / "loaded" / "concat.h5ad"
    )
    assert legacy_h5ad_path(root, "MH7", "assembled_reference") == (
        root / "MH7" / "assembled" / "reference.h5ad"
    )


def test_install_legacy_symlink_relative_and_round_trips(tmp_path: Path):
    from ref_build._internal.layout import install_legacy_symlink

    new_path = tmp_path / "refbuild_MH7_Prostate" / "loaded_concat.h5ad"
    new_path.parent.mkdir(parents=True, exist_ok=True)
    new_path.write_bytes(b"NEW")

    legacy_path = tmp_path / "MH7" / "loaded" / "concat.h5ad"
    install_legacy_symlink(new_path, legacy_path, enabled=True)

    assert legacy_path.is_symlink()
    assert legacy_path.resolve() == new_path.resolve()
    assert not os.path.isabs(os.readlink(legacy_path))


def test_install_legacy_symlink_disabled_is_noop(tmp_path: Path):
    from ref_build._internal.layout import install_legacy_symlink

    new_path = tmp_path / "refbuild_MH7" / "loaded_concat.h5ad"
    new_path.parent.mkdir(parents=True, exist_ok=True)
    new_path.write_bytes(b"X")

    legacy_path = tmp_path / "MH7" / "loaded" / "concat.h5ad"
    install_legacy_symlink(new_path, legacy_path, enabled=False)
    assert not legacy_path.exists()


def test_install_legacy_symlink_preserves_pre_refactor_file(tmp_path: Path):
    from ref_build._internal.layout import install_legacy_symlink

    new_path = tmp_path / "refbuild_MH7" / "loaded_concat.h5ad"
    new_path.parent.mkdir(parents=True, exist_ok=True)
    new_path.write_bytes(b"NEW")

    legacy_path = tmp_path / "MH7" / "loaded" / "concat.h5ad"
    legacy_path.parent.mkdir(parents=True, exist_ok=True)
    legacy_path.write_bytes(b"OLD")

    install_legacy_symlink(new_path, legacy_path, enabled=True)
    assert legacy_path.is_symlink()
    bak = legacy_path.with_suffix(legacy_path.suffix + ".pre_relayout_bak")
    assert bak.exists() and bak.read_bytes() == b"OLD"
