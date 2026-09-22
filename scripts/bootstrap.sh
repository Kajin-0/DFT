#!/usr/bin/env bash
# One-time, fully project-local bootstrap:
#   1. python .venv + mct_dft package
#   2. micromamba static binary (if absent)
#   3. conda-forge Quantum ESPRESSO environment under .local/qe-env
# Nothing is written outside the repository; no sudo, no system changes.
set -euo pipefail
cd "$(dirname "$0")/.."
REPO_ROOT="$PWD"

mkdir -p .tmp .cache/pip .cache/xdg .cache/matplotlib .local work

if [ ! -x .venv/bin/python ]; then
  echo "[bootstrap] creating python venv"
  python3 -m venv .venv
fi
echo "[bootstrap] installing python package"
PIP_CACHE_DIR="$REPO_ROOT/.cache/pip" .venv/bin/pip install -q --upgrade pip
PIP_CACHE_DIR="$REPO_ROOT/.cache/pip" .venv/bin/pip install -q -e ".[alloy,dev]"

if [ ! -x .local-bin/micromamba ]; then
  echo "[bootstrap] downloading micromamba (static, project-local)"
  mkdir -p .local-bin .local/micromamba
  curl -sL https://micro.mamba.pm/api/micromamba/linux-64/latest -o .tmp/micromamba.tar.bz2
  tar -xjf .tmp/micromamba.tar.bz2 -C .local-bin --strip-components=1 bin/micromamba
  rm -f .tmp/micromamba.tar.bz2
fi

if [ ! -x .local/qe-env/bin/pw.x ]; then
  echo "[bootstrap] installing Quantum ESPRESSO (conda-forge) into .local/qe-env"
  MAMBA_ROOT_PREFIX="$REPO_ROOT/.local/micromamba" \
    .local-bin/micromamba create -y -p "$REPO_ROOT/.local/qe-env" -c conda-forge qe
fi

echo "[bootstrap] done. Now: source scripts/env.sh && scripts/check_environment.sh"
