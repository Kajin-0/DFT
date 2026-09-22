#!/usr/bin/env python3
"""Production pipeline for one endpoint system (CdTe or HgTe).

Sequence (all at converged cutoffs from config/systems/<sys>.yaml):
  1. vc-relax (no SOC)        -> a_nosoc
  2. vc-relax (SOC)           -> a_soc
  3. per treatment {nosoc, soc}: scf at own relaxed lattice (k_scf)
  4. bands along SeekPath (from scf charge)
  5. nscf dense (k_dos, tetrahedra) + dos.x + projwfc.x
Stopping is explicit and every stage is provenance-tagged.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]


def sh(cmd: list[str]) -> None:
    print("+", " ".join(map(str, cmd)), flush=True)
    p = subprocess.run(cmd, capture_output=True, text=True)
    (ROOT / ".tmp/last_production.log").write_text(p.stdout + "\n" + p.stderr)
    if p.returncode != 0:
        print(p.stdout[-2000:])
        print(p.stderr[-2000:], file=sys.stderr)
        raise SystemExit(f"failed: {cmd}")
    # surface the compact result dict
    tail = p.stdout.strip().splitlines()
    for line in reversed(tail):
        if line.strip().startswith("{"):
            try:
                print(json.dumps(json.loads(line), indent=1))
            except Exception:
                pass
            break


def relaxed_lattice(system: str, tag: str) -> float:
    """Extract final lattice constant from a vc-relax run dir result.json."""
    res = json.loads(
        (ROOT / "calculations" / system / tag / "result.json").read_text())
    cell = np.asarray(res["pw"]["final_cell_ang"], dtype=float)
    # zincblende lattice parameter from primitive cell: |a1|*sqrt(2)
    a = float(np.linalg.norm(cell[0]) * np.sqrt(2.0))
    return a


def stage(system, stage, *, soc, ecut, rho, kgrid, np_, lat=None,
          tag=None, charge_from=None):
    cmd = [sys.executable, str(ROOT / "scripts/run_stage.py"), system, stage,
           "--ecutwfc", str(ecut), "--ecutrho", str(rho),
           "--kgrid", *[str(k) for k in kgrid], "--np", str(np_)]
    if soc:
        cmd.append("--soc")
    if lat is not None:
        cmd += ["--lattice", f"{lat:.6f}"]
    if tag:
        cmd += ["--tag", tag]
    if charge_from:
        cmd += ["--charge-from", charge_from]
    sh(cmd)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("system", choices=["cdte", "hgte"])
    ap.add_argument("--np", type=int, default=6)
    ap.add_argument("--skip-relax", action="store_true")
    args = ap.parse_args()

    cfg = yaml.safe_load(
        (ROOT / "config/systems" / f"{args.system}.yaml").read_text())
    cs = cfg["converged_settings"]
    ecut, rho = cs["ecutwfc_ry"], cs["ecutrho_ry"]
    k_scf, k_dos = cs["kgrid_scf"], cs["kgrid_dos"]
    a0 = cfg["lattice_start_ang"]

    record = {"system": args.system, "ecutwfc_ry": ecut, "ecutrho_ry": rho,
              "kgrid_scf": k_scf, "kgrid_dos": k_dos}

    for soc in (False, True):
        name = "soc" if soc else "nosoc"
        if not args.skip_relax:
            stage(args.system, "vc-relax", soc=soc, ecut=ecut, rho=rho,
                  kgrid=k_scf, np_=args.np, lat=a0, tag=f"relax_{name}")
        a = relaxed_lattice(args.system, f"relax_{name}")
        record[f"a_relaxed_{name}_ang"] = a
        print(f"[{args.system}] relaxed lattice ({name}): {a:.5f} A")

        stage(args.system, "scf", soc=soc, ecut=ecut, rho=rho, kgrid=k_scf,
              np_=args.np, lat=a, tag=f"production_scf_{name}")
        charge = str(Path("calculations") / args.system
                     / f"production_scf_{name}")
        stage(args.system, "bands", soc=soc, ecut=ecut, rho=rho, kgrid=k_scf,
              np_=args.np, lat=a, tag=f"production_bands_{name}",
              charge_from=charge)
        stage(args.system, "dos", soc=soc, ecut=ecut, rho=rho, kgrid=k_dos,
              np_=args.np, lat=a, tag=f"production_dos_{name}",
              charge_from=charge)

    out = ROOT / "results/provenance" / f"{args.system}_production.json"
    out.write_text(json.dumps(record, indent=2))
    print(json.dumps(record, indent=2))
    return 0


if __name__ == "__main__":
    main()
