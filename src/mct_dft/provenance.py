"""Provenance capture: environment, versions, hashes, run records.

Every production result in results/ must be traceable back to (at minimum)
QE version, pseudopotential SHA-256, input SHA-256, MPI/process count and the
git commit at runtime.
"""
from __future__ import annotations

import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path

from .pseudopotentials import sha256_file

_TRACKED_PY_PKGS = ["numpy", "scipy", "matplotlib", "ase", "icet", "spglib",
                    "seekpath", "pyyaml", "pytest"]


def package_versions() -> dict[str, str]:
    out = {}
    for name in _TRACKED_PY_PKGS:
        try:
            out[name] = importlib_metadata.version(name)
        except importlib_metadata.PackageNotFoundError:
            out[name] = None
    return out


def git_commit(root: Path | None = None) -> str | None:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=str(root or Path.cwd()),
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except Exception:
        return None


def hardware_info() -> dict:
    info = {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
    }
    try:
        info["cpu_logical_count"] = os.cpu_count()
    except Exception:
        pass
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemTotal"):
                    info["mem_total_kb"] = int(line.split()[1])
                if line.startswith("MemAvailable"):
                    info["mem_available_kb"] = int(line.split()[1])
    except OSError:
        pass
    return info


def executable_head(executable: str, env=None) -> str | None:
    """First output lines of a QE binary (contains the version banner)."""
    try:
        proc = subprocess.run([executable], input="", capture_output=True,
                              text=True, timeout=20, env=env)
        head = "\n".join((proc.stdout + proc.stderr).splitlines()[:4])
        return head.strip() or None
    except Exception:
        return None


def write_provenance(record: dict, out_path: str | Path) -> Path:
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    record = dict(record)
    record["written_utc"] = datetime.now(timezone.utc).isoformat()
    with open(out, "w") as fh:
        json.dump(record, fh, indent=2, sort_keys=True, default=str)
    return out


__all__ = [
    "package_versions", "git_commit", "hardware_info",
    "executable_head", "write_provenance", "sha256_file",
]
