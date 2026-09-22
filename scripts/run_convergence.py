#!/usr/bin/env python3
"""Cutoff- and k-grid convergence studies for one system.

For each parameter point:
  * SCF run     -> total energy  (meV/atom stability criterion)
  * NSCF run    -> highest occupied / lowest unoccupied level on the same
                   uniform grid (gap convergence proxy)

Soft saving: every point lands in its own run directory with full
provenance; the aggregated CSVs go to results/tables/ and the plots to
results/figures/.
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

RHO_RATIO_DEFAULT = 6.0   # ecutrho = ratio * ecutwfc (rrkj ultrasoft, PSlibrary)


def run_stage(stage: str, system: str, ecut: float, ecutrho: float,
              kgrid: tuple[int, int, int], tag: str, np_ranks: int,
              charge_from: str | None = None) -> dict:
    cmd = [sys.executable, str(ROOT / "scripts/run_stage.py"), system, stage,
           "--ecutwfc", str(ecut), "--ecutrho", str(ecutrho),
           "--kgrid", *[str(k) for k in kgrid], "--np", str(np_ranks),
           "--tag", tag]
    if charge_from:
        cmd += ["--charge-from", charge_from]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"stage failed: {tag}\n{p.stdout}\n{p.stderr}")
    res = json.loads((ROOT / "calculations" / system / tag / "result.json").read_text())
    return {
        "ecutwfc_ry": ecut, "ecutrho_ry": ecutrho,
        "kgrid": "x".join(map(str, kgrid)),
        "total_energy_ry": res["pw"]["total_energy_ry"],
        "highest_occupied_ev": res["pw"]["highest_occupied_ev"],
        "lowest_unoccupied_ev": res["pw"]["lowest_unoccupied_ev"],
        "gap_ev": res["pw"].get("gap_ev"),
        "wall_s": sum(r["wall_time_s"] for r in res["runs"]),
        "ok": res["pw"]["ok"],
    }


def scan(system: str, scan_type: str, values: list[float], fixed_ecut: float,
         fixed_kgrid: tuple[int, int, int], np_ranks: int) -> list[dict]:
    rows = []
    for v in values:
        if scan_type == "ecut":
            ecut, rho, kg = v, RHO_RATIO_DEFAULT * v, fixed_kgrid
            tag = f"conv_scf_e{v:g}"
        elif scan_type == "ecutrho":
            ecut, rho, kg = fixed_ecut, v, fixed_kgrid
            tag = f"conv_scf_rho{v:g}"
        else:
            kv = int(v)
            ecut, rho, kg = fixed_ecut, RHO_RATIO_DEFAULT * fixed_ecut, (kv, kv, kv)
            tag = f"conv_scf_k{kv}"
        print(f"[{system}/{scan_type}] running {tag}", flush=True)
        for stage in ("scf", "nscf"):
            t = tag.replace("conv_scf", f"conv_{stage}")
            cf = (str(Path("calculations") / system / tag)
                  if stage == "nscf" else None)
            row = run_stage(stage, system, ecut, rho, kg, t, np_ranks,
                            charge_from=cf)
            if stage == "scf":
                rows.append(row)
            else:
                rows[-1].update({
                    "highest_occupied_ev": row["highest_occupied_ev"],
                    "lowest_unoccupied_ev": row["lowest_unoccupied_ev"],
                    "gap_proxy_ev": row.get("gap_ev"),
                    "wall_s": rows[-1]["wall_s"] + row["wall_s"],
                })
    return rows


def write_csv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def plot(rows: list[dict], system: str, scan_type: str, out_png: Path,
         gap_target_mev: float = 5.0, e_target_mev: float = 1.0) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    if scan_type == "ecut":
        x = [r["ecutwfc_ry"] for r in rows]
        xl = "ecutwfc (Ry)"
    elif scan_type == "ecutrho":
        x = [r["ecutrho_ry"] for r in rows]
        xl = "ecutrho (Ry)"
    else:
        x = [int(r["kgrid"].split("x")[0]) for r in rows]
        xl = "k-grid n (n x n x n)"

    e = np.array([r["total_energy_ry"] for r in rows]) * 13605.693122994  # meV
    nat = 2
    e_per_atom = e / nat
    g = np.array([r["gap_proxy_ev"] if r["gap_proxy_ev"] is not None else np.nan
                  for r in rows]) * 1000.0  # meV

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    axes[0].plot(x, e_per_atom - e_per_atom[-1], "o-", color="#0b5394")
    axes[0].axhline(e_target_mev, ls="--", c="grey", lw=0.8,
                    label=f"+{e_target_mev:g} meV/atom target")
    axes[0].axhline(-e_target_mev, ls="--", c="grey", lw=0.8)
    axes[0].set_xlabel(xl)
    axes[0].set_ylabel("E - E(ref)  [meV/atom]")
    axes[0].set_title(f"{system}: total-energy convergence")
    axes[0].legend(fontsize=8)

    axes[1].plot(x, g - g[-1], "s-", color="#990000")
    axes[1].axhline(gap_target_mev, ls="--", c="grey", lw=0.8,
                    label=f"+{gap_target_mev:g} meV target")
    axes[1].axhline(-gap_target_mev, ls="--", c="grey", lw=0.8)
    axes[1].set_xlabel(xl)
    axes[1].set_ylabel("gap_proxy - gap_proxy(ref)  [meV]")
    axes[1].set_title(f"{system}: gap-proxy convergence (uniform grid HO/LU)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("system", choices=["cdte", "hgte"])
    ap.add_argument("--scan", choices=["ecut", "kgrid", "ecutrho"], required=True)
    ap.add_argument("--values", type=float, nargs="+", required=True)
    ap.add_argument("--fixed-ecut", type=float, default=70.0)
    ap.add_argument("--fixed-kgrid", type=int, nargs=3, default=[8, 8, 8])
    ap.add_argument("--np", type=int, default=6)
    args = ap.parse_args()

    rows = scan(args.system, args.scan, args.values, args.fixed_ecut,
                tuple(args.fixed_kgrid), args.np)
    write_csv(rows, ROOT / "results/tables"
              / f"{args.system}_{args.scan}_convergence.csv")
    plot(rows, args.system, args.scan,
         ROOT / "results/figures" / f"{args.system}_{args.scan}_convergence.png")
    print(json.dumps({"rows": rows}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
