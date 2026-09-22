#!/usr/bin/env bash
# Process-local environment for this project. Source it; nothing touches the
# user's shell configuration or any system directory.
#   source scripts/env.sh
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export WORKSPACE_ROOT="$REPO_ROOT"
export MCT_DFT_ROOT="$REPO_ROOT"
export TMPDIR="$REPO_ROOT/.tmp"
export PIP_CACHE_DIR="$REPO_ROOT/.cache/pip"
export XDG_CACHE_HOME="$REPO_ROOT/.cache/xdg"
export MPLCONFIGDIR="$REPO_ROOT/.cache/matplotlib"
export MAMBA_ROOT_PREFIX="$REPO_ROOT/.local/micromamba"
export QE_ENV="$REPO_ROOT/.local/qe-env"
export PATH="$QE_ENV/bin:$REPO_ROOT/.venv/bin:$PATH"
export VIRTUAL_ENV="$REPO_ROOT/.venv"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
# pseudopotential directory (downloaded UPFs live here, not committed)
export PSEUDO_DIR="$REPO_ROOT/pseudopotentials"
export QE_PSEUDO="$PSEUDO_DIR"
mkdir -p "$TMPDIR" "$PIP_CACHE_DIR" "$XDG_CACHE_HOME" "$MPLCONFIGDIR"
