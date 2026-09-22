# AGENTS.md — working rules for this repository

## Mission
Reproducible first-principles DFT (Quantum ESPRESSO) workflow for HgCdTe
photodetector materials. Scientific honesty beats output volume.

## Hard rules
- Containment: everything under the repo root. Never write `/tmp`, `/etc`,
  `/usr`, `$HOME` config. Project-local caches: `.tmp/`, `.cache/`,
  `.local/`, `.venv/`, `work/`. Set env via `source scripts/env.sh`.
- No `sudo`, no `apt`, no global pip/conda. QE comes from the project-local
  micromamba env at `.local/qe-env`.
- Do NOT commit: `.venv/`, `.cache/`, `.tmp/`, `.local*/`, `work/*`,
  QE `out/` scratch dirs, `*.UPF`, `results/raw/` (see `.gitignore`).
- Never fabricate a DFT number. Parsed results must carry provenance
  (git commit, QE version, PP sha256, input sha256, MPI ranks).
- Gap/cutoff logic: `λc` only for a positive normal-semiconductor gap;
  inverted/overlapping band orderings must be reported honestly (NaN/None).
- Empirical Eg(x,T) (Hansen) is a benchmark, never fitted to, and T≠0 K
  values are context, not static-lattice DFT comparisons.
- Ordinary git only (`git status/add/commit/push`, `git pull --ff-only`).
  No `gh`, no force-push.

## Commands that work
```bash
source scripts/env.sh
python -m pytest            # 38 fast tests, no QE
scripts/check_environment.sh
python scripts/fetch_pseudopotentials.py
python scripts/run_stage.py cdte scf --ecutwfc 60 --ecutrho 480 --kgrid 10 10 10
```

## Machine (recorded in results/provenance/environment.json)
8-core AMD EPYC 9354P VPS, 31 GB RAM, ~200 GB free disk. QE 7.5 (conda-
forge), Open MPI 5.0.10. Default MPI ranks: 6 (≤75% of cores).

## Tests
`python -m pytest` runs the QE-free suite. Tests requiring real QE runs are
marked `@pytest.mark.qe` and excluded by default.
