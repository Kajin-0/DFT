"""Command-line entry points for the mct_dft package."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import detector
from .parser import parse_pw_output
from .runner import workspace_root, ensure_within_workspace


def _cmd_empirical(args) -> int:
    for row in detector.empirical_series(args.x, temperatures=args.temps):
        lam = row["lambda_c_um"]
        print(
            f"x={row['x']:.3f}  T={row['T_K']:6.1f} K  "
            f"Eg={row['Eg_eV']:+.8f} eV  "
            f"lambda_c={f'{lam:.5f} um' if lam else 'n/a (Eg<=0)'}"
        )
    return 0


def _cmd_parse(args) -> int:
    path = ensure_within_workspace(args.output)
    text = Path(path).read_text(errors="replace")
    res = parse_pw_output(text)
    out = {
        "ok": res.ok, "converged": res.converged, "job_done": res.job_done,
        "calculation": res.calculation,
        "total_energy_Ry": res.total_energy_ry,
        "total_energy_eV": res.total_energy_ev,
        "fermi_eV": res.fermi_ev,
        "highest_occupied_eV": res.highest_occupied_ev,
        "n_iterations": res.nscf_iterations,
        "n_bands": res.n_bands, "n_kpts": res.n_kpts, "nat": res.nat,
        "noncollinear": res.noncollinear, "spin_orbit": res.spin_orbit,
        "errors": res.errors,
    }
    print(json.dumps(out, indent=2))
    return 0 if res.ok else 2


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="mct-dft")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("empirical", help="empirical Eg(x,T) and cutoff table")
    p.add_argument("--x", type=float, default=0.2)
    p.add_argument("--temps", type=float, nargs="+", default=[0.0, 77.0, 300.0])
    p.set_defaults(func=_cmd_empirical)

    p = sub.add_parser("parse", help="parse a pw.x output file")
    p.add_argument("output")
    p.set_defaults(func=_cmd_parse)

    p = sub.add_parser("root", help="print workspace root")
    p.set_defaults(func=lambda a: (print(workspace_root()) or 0))

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
