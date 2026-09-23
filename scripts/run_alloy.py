#!/usr/bin/env python3
"""Phase E/F driver for the Hg0.8Cd0.2Te SQS alloy.

Model A: ideal zincblende-site SQS (unrelaxed), lattice from Vegard mixing
         of the *DFT-relaxed* endpoint lattice constants.
Model B: internal-coordinate relaxation (fixed cell) + SCF.

Runs SCF with and without SOC at the documented alloy cutoffs
(80/480 Ry — endpoint tests showed gaps converged to <0.1 meV there; the
120/720 values are only needed for sub-meV total-energy convergence).
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sh(cmd: list[str]) -> None:
    print("+", " ".join(map(str, cmd)), flush=True)
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        print(p.stdout[-3000:])
        print(p.stderr[-2000:], file=sys.stderr)
        raise SystemExit(f"failed: {' '.join(map(str, cmd))}")


def endpoint_a(system: str, soc: bool) -> float:
    tag = "relax_soc" if soc else "relax_nosoc"
    res = json.loads((ROOT / "calculations" / system / tag
                      / "result.json").read_text())
    import numpy as np
    cell = np.asarray(res["pw"]["final_cell_ang"], float)
    return float(np.linalg.norm(cell[0]) * np.sqrt(2.0))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--struct", required=True,
                    help="repo-relative path to the SQS structure (extxyz)")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--np", type=int, default=6)
    ap.add_argument("--ecut", type=float, default=80.0)
    ap.add_argument("--ecutrho", type=float, default=480.0)
    ap.add_argument("--kgrid", type=int, nargs=3, default=[2, 2, 2])
    ap.add_argument("--real-space-paw", action="store_true")
    ap.add_argument("--relax", action="store_true",
                    help="Model B: ionic relaxation before the SOC SCF")
    args = ap.parse_args()

    e, rho = args.ecut, args.ecutrho
    kg = [str(k) for k in args.kgrid]
    s = args.struct

    def stage(stage_, soc, tag, extra=None, check_energy=True):
        cmd = [sys.executable, str(ROOT / "scripts/run_stage.py"), "hgcdte",
               stage_, "--ecutwfc", str(e), "--ecutrho", str(rho),
               "--kgrid", *kg, "--np", str(args.np), "--struct", s,
               "--tag", tag]
        if soc:
            cmd.append("--soc")
        if extra:
            cmd += extra
        sh(cmd)
        if check_energy:
            res = json.loads((ROOT / "calculations/hgcdte" / tag
                              / "result.json").read_text())
            if not res["pw"]["ok"]:
                raise SystemExit(f"stage {tag} failed validation")

    stage("scf", False, f"{args.tag}_scf_nosoc")
    stage("scf", True, f"{args.tag}_scf_soc")
    if args.relax:
        stage("relax", True, f"{args.tag}_relax_soc")
    return 0


if __name__ == "__main__":
    main()
