#!/usr/bin/env python3
"""Empirical benchmark figure: Eg(x,T) curves with the x=0.20 reference."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mct_dft.detector import hansen_eg, cutoff_wavelength_um  # noqa: E402

x = np.linspace(0.15, 0.35, 400)
fig, ax = plt.subplots(figsize=(5.6, 4.0))
for T, ls, lab in ((0.0, "-", "T = 0 K (static-lattice reference)"),
                   (77.0, "--", "T = 77 K"),
                   (300.0, ":", "T = 300 K")):
    ax.plot(x, [hansen_eg(xx, T) for xx in x], ls, label=lab)
for T, ls, lab in ((0.0, "-", None), (77.0, "--", None), (300.0, ":", None)):
    eg = hansen_eg(0.2, T)
    ax.plot([0.2], [eg], "o", color="k", ms=4)
ax.axvline(0.2, color="0.7", lw=0.5)
ax.set_xlabel("x in Hg$_{1-x}$Cd$_x$Te")
ax.set_ylabel("Eg (eV)")
ax.set_title("Empirical HgCdTe band gap (Hansen)")
ax.legend(fontsize=8)
fig.tight_layout()
out = ROOT / "results/figures/empirical_hansen.png"
fig.savefig(out, dpi=200)
print(f"wrote {out}")
