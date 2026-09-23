#!/usr/bin/env bash
# Remaining production queue, sequential, ~6 MPI ranks each.
set -uo pipefail
cd "$(dirname "$0")/.."
source scripts/env.sh
PY=.venv/bin/python

set -x
# --- CdTe bands: both treatments, fixed nbnd + full-acc SOC ----------------
$PY scripts/run_stage.py cdte bands --ecutwfc 120 --ecutrho 720 --kgrid 10 10 10 \
  --np 6 --lattice 6.623044 --tag production_bands_nosoc --force \
  --charge-from calculations/cdte/production_scf_nosoc
$PY scripts/run_stage.py cdte bands --ecutwfc 120 --ecutrho 720 --kgrid 10 10 10 \
  --np 6 --soc --lattice 6.626419 --tag production_bands_soc --force \
  --charge-from calculations/cdte/production_scf_soc

# --- CdTe SOC DOS (tetrahedra, 14^3) ---------------------------------------
$PY scripts/run_stage.py cdte dos --ecutwfc 120 --ecutrho 720 --kgrid 14 14 14 \
  --np 6 --soc --lattice 6.626419 --tag production_dos_soc --force \
  --charge-from calculations/cdte/production_scf_soc

# --- HgTe full production (relax, scf, bands, dos; SOC + no-SOC) -----------
$PY scripts/run_production.py hgte --np 6

echo QUEUE1_DONE
