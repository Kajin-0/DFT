#!/usr/bin/env python3
"""No-SOC vs SOC band overlay for a system (from two existing runs)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    system, out_png = sys.argv[1], sys.argv[2]
    runs = {}
    for tag, style, label in (("nosoc", dict(color="#0b5394", lw=1.0),
                               "PBE (scalar-rel., no SOC)"),
                              ("soc", dict(color="#990000", lw=1.0),
                               "PBE + SOC")):
        p = ROOT / "calculations" / system / f"production_bands_{tag}"
        z = np.load(p / "bands.npz")
        runs[tag] = (z, style, label, p)

    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    for tag, (z, style, label, p) in runs.items():
        kpts, eigs = z["kpts"], z["eigs_ev"]
        d = np.linalg.norm(np.diff(kpts, axis=0), axis=1)
        x = np.concatenate([[0.0], np.cumsum(d)])
        # align each set's VBM to 0 using its own band counting
        nv = 18 if tag == "soc" else 9
        ref = eigs[:, nv - 1].max()
        for ib in range(eigs.shape[1]):
            ax.plot(x, eigs[:, ib] - ref, **style,
                    label=label if ib == 0 else None)
    ax.axhline(0, color="k", lw=0.5, ls="--")
    ax.set_ylabel("E - VBM (eV)")
    ax.set_xlabel("k (path)")
    ax.set_ylim(-8, 6)
    ax.legend(fontsize=8)
    ax.set_title(f"{system.upper()}: SOC comparison "
                 "(VBM-aligned; gaps differ? see labels)", fontsize=10)
    fig.tight_layout()
    out = ROOT / out_png
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
