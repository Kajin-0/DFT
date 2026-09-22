# Methodology

_Last updated: Phase A + start of Phase B. Numbers quoted here must trace to
raw QE outputs under `calculations/`; if a claim cannot be traced, it doesn't
belong here._

## 1. Scope

First-principles DFT (Quantum ESPRESSO 7.5, plane waves + ultrasoft
pseudopotentials) of the Hg(1-x)Cd(x)Te photodetector materials chain:

1. zincblende CdTe primitive cell — convergence, structure, bands, DOS,
   effective masses, optics;
2. zincblende HgTe — same, with explicit band-topology care (inversion);
3. Hg0.8Cd0.2Te explicit SQS supercell (cation-sublattice alloying only);
4. derived optical + ideal photovoltaic figures of merit.

## 2. Software environment (fully project-local)

* Python 3.13 venv in `.venv/` — numpy, scipy, matplotlib, ase, spglib,
  seekpath, icet, pyyaml, pytest. No global installs.
* Quantum ESPRESSO 7.5 from conda-forge, installed under `.local/qe-env`
  via a project-local static micromamba binary. MPI: Open MPI 5.0.10.
* Machine: 8-core AMD EPYC 9354P VPS, 31 GB RAM, ~200 GB free disk.
* MPI ranks: ≤ 6 (= 75% of 8 logical cores), 1 OMP thread per rank.
  Never oversubscribed.
* Exact versions/hashes: `results/provenance/environment.json` and each
  run's `provenance.json`.

## 3. Structures

* Primitive zincblende cells built with ASE (`ase.build.bulk`,
  `crystalstructure='zincblende'`): cation at (0,0,0), anion at
  (1/4,1/4,1/4), space group F-43m (No. 216) verified with spglib.
  Structural checks (stoichiometry, volume, 4-fold NN coordination,
  NN distance -> implied a_conv) are stored in every run's `result.json`.
* Starting lattice constants are literature values (CdTe 6.482 A,
  HgTe 6.453 A). The DFT equilibrium lattice parameter is obtained from
  a vc-relax/EOS step with converged settings and recorded separately;
  the starting value is never reported as a result.

## 4. Pseudopotentials

PSlibrary 1.0.0 (rrkjus ultrasoft), **PBE for all elements**
(no functional mixing), one generator family:

| role | elements | files | purpose |
|------|----------|-------|---------|
| `soc`   | Hg, Cd, Te | `*.rel-pbe-*` fully relativistic, `has_so="T"` | production SOC (`noncolin=.true.`, `lspinorb=.true.`) |
| `nosoc` | Hg, Cd, Te | `*.pbe-*` scalar-relativistic twins | SOC-off control |

Note: QE >= 7 refuses fully relativistic PPs with `lspinorb=.false.`
(`average_pp` error) — hence the paired scalar-relativistic files from the
*same family* for an honest SOC comparison. Headers (not filenames) are
parsed to verify relativistic="full" and has_so; SHA-256 pins recorded in
`pseudopotentials/manifest.json`. UPF files are git-ignored; retrieval is
reproducible via `scripts/fetch_pseudopotentials.py`.

## 5. Convergence policy

* Successive-step total-energy criterion: **< 1 meV/atom**.
* Gap stability criterion: **< 5 meV** (uniform-grid HO/LU proxy from NSCF
  eigenvalue listings).
* ecutrho = 6 x ecutwfc baseline for rrkjUS potentials; the chosen pair is
  then used consistently for structure/bands/DOS. Any deviation is
  documented per system in `reports/results.md`.
* k-grids scanned 4x4x4 .. 12x12x12 for primitive cells; chosen at the
  smallest mesh meeting both criteria.
* `conv_thr = 1e-10 Ry` everywhere; `occupations='fixed'` for SCF,
  tetrahedra for NSCF/DOS.

## 6. Electronic structure protocol (per endpoint)

1. ecut convergence scan (SCF + NSCF) at fixed k-grid;
2. k-grid convergence at fixed ecut;
3. vc-relax at final settings -> DFT lattice constant;
4. dense NSCF + dos.x (+ projwfc.x for PDOS);
5. `bands` along SeeK-path FCC-zincblende path;
6. no-SOC control (scalar-relativistic PPs) for the SOC comparison;
7. directional effective masses from dense Gamma-centred segments.

## 7. Spin-orbit coupling

Production SOC runs require `noncolin=.true.`, `lspinorb=.true.` **and**
UPFs validated as fully relativistic with spin-orbit projectors; the Python
validator (`mct_dft.pseudopotentials.validate_soc_set`) refuses otherwise.
`Delta_SO` and gap changes are computed as SOC-minus-control at identical
cutoffs/grids/lattices.

## 8. Band-edge bookkeeping

`mct_dft.bands.band_edges` never enforces a positive gap. If the computed
ordering is overlapping/inverted (VBM > CBM across k), the signed gap is
reported and status marked `inverted/overlap`; `cutoff_um_from_gap` returns
None instead of a wavelength. This is essential for HgTe and possibly the
alloy.

## 9. Effective masses

Parabolic fits `E = E0 + hbar^2 k^2 / 2m*` on signed, Cartesian-scaled
segments around Gamma in [100], [110], [111], at window half-widths
0.01/0.02/0.03 1/A; R^2 reported; non-parabolic cases are flagged, not
hidden. All unit conversions are centralised in `mct_dft/constants.py`
(CODATA via scipy).

## 10. Alloy (Phase E; planned)

* Te anion sublattice untouched; Hg/Cd occupy the FCC **cation** sublattice
  only, with exactly 4 Cd of 20 cations (x = 0.20 exact).
* SQS generated with icet; correlation metrics + random seed recorded.
* Model A: ideal unrelaxed lattice with Vegard-interpolated a; Model B:
  cell/internal-coordinate relaxed. Compared directly.
* Primitive-cell band pictures do not survive supercell band folding;
  alloy results emphasise DOS, band-edge eigenvalues and spectral gap,
  not decorated high-symmetry spaghetti plots.

## 11. Physics of detector quantities

Where computed: eps(omega) -> n+ik -> alpha, R, ideal absorptance with a
documented thickness (default 10 um), ideal PV responsivity
`R = eta q lambda / h c` with eta = (1-R)(1-exp(-alpha d)). These are
**ideal** numbers assuming perfect collection and unity gain; they are not
measured HgCdTe detector performance. Lifetimes/mobility/noise/D* are
explicitly out of scope of the baseline and labelled as future work
(see README "What remains").

## 12. Reproducibility contract

Every production result row carries: git commit, QE version, python package
versions, PP sha256 set, input sha256, MPI ranks, run wall time, parser
verdict. Failed/convergence-NOT-achieved/OOM runs are never parsed into
results tables (see `mct_dft.parser` failure markers and tests).
