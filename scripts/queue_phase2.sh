#!/usr/bin/env bash
# Phase 2 queue (run AFTER queue4/alloy completes): optics + PDOS + masses.
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/env.sh
PY=.venv/bin/python
set -x

# 1. CdTe optics via epsilon.x (NC SG15 branch, no SOC -- epsilon.x does not
#    support USPP/PAW): scf at relaxed lattice, then eps on 6^3 shifted grid
$PY scripts/run_stage.py cdte scf --ecutwfc 60 --ecutrho 240 --kgrid 6 6 6 \
  --np 6 --lattice 6.623044 --tag nc_test_scf --pp-role nc
$PY scripts/run_stage.py cdte eps --ecutwfc 60 --ecutrho 240 --kgrid 6 6 6 \
  --nbnd 50 --lattice 6.623044 --tag optics_nc --np 6 --pp-role nc \
  --charge-from calculations/cdte/nc_test_scf

# 2. PDOS (orbital character) for CdTe, both treatments
$PY scripts/run_stage.py cdte pdos --ecutwfc 120 --ecutrho 720 --kgrid 10 10 10 \
  --np 6 --lattice 6.623044 --tag production_pdos_nosoc \
  --charge-from calculations/cdte/production_scf_nosoc
$PY scripts/run_stage.py cdte pdos --ecutwfc 120 --ecutrho 720 --kgrid 10 10 10 \
  --np 6 --soc --lattice 6.626419 --tag production_pdos_soc \
  --charge-from calculations/cdte/production_scf_soc

# 3. SOC effective masses (electron CB; holes) for CdTe at the SOC lattice
$PY scripts/fit_effective_mass.py cdte --band cb --soc --ecut 120 \
  --ecutrho 720 --lattice 6.626419 --np 6 \
  --charge-from calculations/cdte/production_scf_soc || true
$PY scripts/fit_effective_mass.py cdte --band vb --soc --ecut 120 \
  --ecutrho 720 --lattice 6.626419 --np 6 \
  --charge-from calculations/cdte/production_scf_soc || true

echo QUEUE5_DONE
