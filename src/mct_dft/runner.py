"""Safe process execution and workspace path containment.

Every path this project writes must resolve inside WORKSPACE_ROOT. The guard
below is enforced by all runners and input writers.
"""
from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path


def workspace_root() -> Path:
    """Repository root. Overridable for tests via MCT_DFT_ROOT."""
    env = os.environ.get("MCT_DFT_ROOT")
    if env:
        return Path(env).resolve()
    # src/mct_dft/runner.py -> repo root is 3 levels up
    return Path(__file__).resolve().parents[2]


def ensure_within_workspace(path: os.PathLike | str, root: Path | None = None) -> Path:
    """Resolve *path* and assert it lives inside the workspace.

    Rejects absolute paths outside the repo and any '..' escape attempts.
    Raises ValueError on violation; returns the resolved Path on success.
    """
    root = root or workspace_root()
    p = Path(path)
    if not p.is_absolute():
        p = root / p
    try:
        resolved = p.resolve()
    except (OSError, RuntimeError) as exc:  # broken symlink loops etc.
        raise ValueError(f"cannot resolve path {path!r}: {exc}") from exc
    if resolved != root and root not in resolved.parents:
        raise ValueError(
            f"refusing path outside workspace: {resolved} (root={root})"
        )
    return resolved


def run_command(
    argv: list[str],
    cwd: os.PathLike | str,
    stdout_path: os.PathLike | str,
    stderr_path: os.PathLike | str | None = None,
    env: dict | None = None,
    timeout_s: float | None = None,
    root: Path | None = None,
) -> dict:
    """Run argv (no shell) capturing stdout/stderr to files inside workspace.

    Returns a small record with timing and exit status. The command itself is
    never string-interpolated into a shell.
    """
    root = root or workspace_root()
    cwd = ensure_within_workspace(cwd, root)
    stdout_path = ensure_within_workspace(stdout_path, root)
    stderr_path = ensure_within_workspace(stderr_path or Path(str(stdout_path) + ".err"), root)
    stdout_path.parent.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    with open(stdout_path, "w") as out, open(stderr_path, "w") as err:
        proc = subprocess.run(
            list(argv), cwd=str(cwd), stdout=out, stderr=err,
            env=env, timeout=timeout_s, check=False,
        )
    dt = time.time() - t0
    return {
        "argv": list(argv),
        "cwd": str(cwd),
        "returncode": proc.returncode,
        "wall_time_s": dt,
        "stdout": str(stdout_path),
        "stderr": str(stderr_path),
    }
