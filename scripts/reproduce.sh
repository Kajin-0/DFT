#!/usr/bin/env bash
# End-to-end reproduction. Assumes: scripts/bootstrap.sh has been run once,
# scripts/env.sh active. Only run stages that are computationally affordable,
# in the order used to produce the committed results.
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/env.sh

python -m pytest                      # fast unit tests

# --- Phase A: benchmark -----------------------------------------------------
.venv/bin/python -m mct_dft.cli empirical --x 0.199 --temps 0 77 300
.venv/bin/python -m mct_dft.cli empirical --x 0.200 --temps 0 77 300

# --- Phase B/C: convergence (moderate cost) --------------------------------
# CdTe
.venv/bin/python scripts/run_convergence.py cdte --scan ecut \
  --values 40 50 60 70 80 90 100 110 120 --fixed-kgrid 8 8 8 --np 6
.venv/bin/python scripts/run_convergence.py cdte --scan ecutrho \
  --values 480 560 640 720 --fixed-ecut 80 --fixed-kgrid 8 8 8 --np 6
.venv/bin/python scripts/run_convergence.py cdte --scan kgrid \
  --values 4 6 8 10 12 --fixed-ecut 80 --np 6
# HgTe
.venv/bin/python scripts/run_convergence.py hgte --scan ecut \
  --values 40 50 60 70 80 90 100 110 120 --fixed-kgrid 8 8 8 --np 6
.venv/bin/python scripts/run_convergence.py hgte --scan ecutrho \
  --values 480 560 640 720 800 --fixed-ecut 80 --fixed-kgrid 8 8 8 --np 6
.venv/bin/python scripts/run_convergence.py hgte --scan kgrid \
  --values 4 6 8 10 12 --fixed-ecut 80 --np 6

# --- Phase B/C: production (longer) ----------------------------------------
.venv/bin/python scripts/run_production.py cdte --np 6
.venv/bin/python scripts/run_production.py hgte --np 6

echo "reproduce.sh finished. See results/ and reports/."
