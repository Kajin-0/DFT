"""Workspace containment: no writes may escape the repository."""
import os

import pytest

from mct_dft.runner import ensure_within_workspace, workspace_root


def test_inside_paths_ok(tmp_path, monkeypatch):
    root = tmp_path / "repo"
    root.mkdir()
    monkeypatch.setenv("MCT_DFT_ROOT", str(root))
    assert ensure_within_workspace("calculations/cdte/scf.in", root) == (
        root / "calculations/cdte/scf.in"
    )
    assert ensure_within_workspace(root / "results") == root / "results"


def test_escape_rejected(tmp_path):
    root = tmp_path / "repo"
    root.mkdir()
    for bad in ("../../x", "/etc/passwd", "/tmp/x", str(tmp_path / "other")):
        with pytest.raises(ValueError):
            ensure_within_workspace(bad, root)


def test_workspace_root_matches_repo():
    root = workspace_root()
    assert (root / "pyproject.toml").exists()
