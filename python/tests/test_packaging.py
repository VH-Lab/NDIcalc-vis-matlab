"""Checks that the installed ndicalc_vis_matlab carries all of ndi_common."""

from __future__ import annotations

import os
from importlib.metadata import entry_points
from pathlib import Path

import pytest

import ndicalc_vis_matlab

REPO_ROOT = Path(__file__).resolve().parents[2]
KEYS = ("database_documents", "schema_documents")


def _repo_json_files() -> list[Path]:
    files = []
    for key in KEYS:
        root = REPO_ROOT / "ndi_common" / key
        files.extend(p.relative_to(REPO_ROOT / "ndi_common") for p in root.rglob("*.json"))
    return sorted(files)


def test_toolbox_dir_contains_ndi_common():
    toolbox = ndicalc_vis_matlab.toolbox_dir()
    assert isinstance(toolbox, Path)
    assert toolbox.is_dir()
    for key in KEYS:
        assert (toolbox / "ndi_common" / key).is_dir()


def test_every_json_file_is_installed():
    expected = _repo_json_files()
    assert expected, f"no JSON files found under {REPO_ROOT / 'ndi_common'}"
    installed = ndicalc_vis_matlab.toolbox_dir() / "ndi_common"
    missing = [str(rel) for rel in expected if not (installed / rel).is_file()]
    assert not missing, f"missing from the installed package: {missing}"
    differ = [
        str(rel)
        for rel in expected
        if (installed / rel).read_bytes() != (REPO_ROOT / "ndi_common" / rel).read_bytes()
    ]
    assert not differ, f"installed copies differ from the repository: {differ}"


def test_entry_point_is_discoverable():
    matches = [ep for ep in entry_points(group="ndi.calc") if ep.name == "NDIcalc-vis-matlab"]
    assert len(matches) == 1
    ep = matches[0]
    assert ep.value == "ndicalc_vis_matlab:toolbox_dir"
    func = ep.load()
    assert callable(func)
    assert func() == ndicalc_vis_matlab.toolbox_dir()


def test_wheel_install_uses_packaged_copy():
    # Set by CI when it installs the built wheel: there toolbox_dir() must be
    # the installed package, not the repository it fell back to.
    if os.environ.get("NDICALC_EXPECT_WHEEL") != "1":
        pytest.skip("NDICALC_EXPECT_WHEEL is not set")
    package_dir = Path(ndicalc_vis_matlab.__file__).resolve().parent
    assert ndicalc_vis_matlab.toolbox_dir() == package_dir
    assert REPO_ROOT not in package_dir.parents
