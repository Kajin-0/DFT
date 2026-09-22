# HgCdTe first-principles DFT workflow (Quantum ESPRESSO)

Goal: trace the physics chain

```
atomic arrangement → electronic structure → bandgap/topology
→ effective masses / DOS → optical properties → ideal detector quantities
```

for Hg(1-x)Cd(x)Te (MCT) photodetector materials, starting from the
zincblende endpoints **CdTe** and **HgTe** and building an explicit-atom
**Hg0.8Cd0.2Te** SQS alloy model (x = 0.20).

## Status: what has actually been calculated

See `reports/results.md` (auto-updated numbers, all traceable to raw QE
output under `calculations/`). Nothing in this repository is handed in from
literature as a "calculation"; every derived number carries provenance
(input SHA-256, PP SHA-256, QE version, MPI ranks, git commit).

* ✅ Phase A: repository, package, tests, empirical benchmark (Hansen
  Eg(x,T)), environment audit, PP fetch + SOC validation
* ⏳ Phase B/C/D: CdTe / HgTe convergence, structure, bands, DOS, masses
* ⏳ Phase E/F: x=0.20 SQS alloy
* ⏳ Phase G: optics + ideal detector chain

## Quick start (fully self-contained, no sudo, no system changes)

```bash
git clone git@github.com:Kajin-0/DFT.git && cd DFT
bash scripts/bootstrap.sh          # .venv + micromamba + conda-forge QE
source scripts/env.sh              # process-local env vars (tmp/cache/paths)
scripts/check_environment.sh       # writes results/provenance/environment.json
python scripts/fetch_pseudopotentials.py
python -m pytest                   # fast, QE-free unit tests
```

All artifacts stay inside the repository (`.venv/`, `.local/`, `.cache/`,
`.tmp/`, `work/`). Nothing touches `/etc`, `/usr`, your shell rc files, or
system packages.

## Reproduce a calculation

```bash
source scripts/env.sh
# CdTe SCF at chosen converged settings (example)
python scripts/run_stage.py cdte scf  --ecutwfc 60 --ecutrho 480 --kgrid 10 10 10
python scripts/run_stage.py cdte bands --soc --ecutwfc 60 --ecutrho 480
python scripts/run_stage.py cdte dos   --soc --ecutwfc 60 --ecutrho 480 --kgrid 14 14 14
```

See `reports/methodology.md` for the full reproducible sequence, including
the converged cutoffs/k-grids selected by the convergence studies.

## What is first-principles vs empirical?

| Output | Type |
|---|---|
| Empirical Eg(x,T) and λc (Hansen) | **EMPIRICAL BENCHMARK** (`mct_dft.detector`) |
| PBE(+SOC) structures, gaps, masses, DOS | first-principles DFT (QE) |
| Optical ε(ω) → n, κ, α, R, η, ideal Rλ | DFT + analytic optics |
| Carrier lifetimes, mobility, D*, noise | **not yet first-principles** (future work; never silently substituted) |

Key definition issue for narrow/negative gap materials: `λc = 1.23984/Eg` is
returned **only for a positive, normal semiconductor gap**. Inverted or
overlapping bands (realistic for HgTe and maybe the alloy) are reported as
such — the code will not fabricate a cutoff wavelength.

## Why HgCdTe is hard

* Very small band gaps (~0.06–0.15 eV at x≈0.2) — convergence must be ~meV.
* HgTe is band-inverted/semimetallic; spin-orbit controls the ordering.
* Semilocal DFT (PBE) band-gap error may exceed the physical gap itself:
  discrepancies vs the empirical Hansen gap are expected and are reported,
  not tuned away.
* The T=0 static-lattice DFT gap is compared against the T=0 empirical
  value; finite-T values (77/300 K) are detector context only.

## Repository layout

```
config/        YAML configs (defaults + per-system)
src/mct_dft/   Python package (constants, structures, QE I/O, parser,
               bands, effective mass, optics, detector, provenance, CLI)
scripts/       bootstrap / env / audit / PP fetch / run_stage / reproduce
calculations/  per-system run dirs (inputs committed; raw scratch ignored)
pseudopotentials/  manifest.json (+UFP files, git-ignored, sha256-pinned)
results/       tables, figures, provenance
reports/       methodology.md, results.md
tests/         QE-free unit tests (pytest)
work/          transient scratch (ignored)
```

## Pseudopotentials

PSlibrary 1.0.0, PBE, rrkjus ultrasoft, one matched family for Hg/Cd/Te:
fully-relativistic `has_so` files for SOC production plus the scalar-
relativistic twins as SOC-off controls. Verified from UPF headers, SHA-256
pinned in `pseudopotentials/manifest.json`. Never selected by filename.
