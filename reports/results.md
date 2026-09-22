# Results (auto-maintained; every number traces to `calculations/` raw output)

## Empirical benchmark (computed by `mct_dft.detector`, not DFT)

Hansen-type Eg(x,T) = -0.302 + 1.93x - 0.81x^2 + 0.832x^3 + 5.35e-4(1-2x)T:

| x     | T (K) | Eg (eV)    | lambda_c (um) |
|-------|-------|------------|---------------|
| 0.199 |   0   | 0.05655    | 21.9248       |
| 0.199 |  77   | 0.08135    | 15.2410       |
| 0.199 | 300   | 0.15317    | 8.0945        |
| 0.200 |   0   | 0.058256   | 21.28265      |
| 0.200 |  77   | 0.082973   | 14.94272      |
| 0.200 | 300   | 0.154556   | 8.02196       |

Unit tests pin these to printed precision (tests/test_detector.py).

> Reminder: static-lattice ground-state DFT compares against the **T=0 K**
> value; 77/300 K are detector context, not DFT validation targets.

## Phase B/C/D status

_Pending convergence scans — will be filled from results/tables/*.csv and
calculations/*/result.json only._

## Environment (from results/provenance/environment.json)

QE 7.5 (conda-forge, project-local), Open MPI 5.0.10, Python 3.13.7,
8-core EPYC 9354P, 31 GB RAM. PPs: PSlibrary 1.0.0 PBE rrkjus,
fully-relativistic/soc + scalar-relativistic controls (manifest.json).
