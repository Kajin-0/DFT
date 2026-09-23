"""Quantum ESPRESSO input-file writers (pw.x, bands.x, dos.x, epsilon.x).

Plain-text generation, no hidden state. Strings are assembled from explicit
python dicts; numbers are formatted with full precision.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from ase import Atoms

from .runner import ensure_within_workspace


def _fmt_block(name: str, params: dict) -> str:
    lines = [f"&{name}"]
    for key, val in params.items():
        if val is None:
            continue
        if isinstance(val, bool):
            val = ".true." if val else ".false."
        elif isinstance(val, str):
            val = f"'{val}'"
        lines.append(f"   {key} = {val}")
    lines.append("/")
    return "\n".join(lines)


def pw_input(
    atoms: Atoms,
    pseudopotentials: dict[str, str],  # symbol -> upf filename
    ecutwfc_ry: float,
    ecutrho_ry: float | None = None,
    calculation: str = "scf",
    kgrid: tuple[int, int, int] | None = (8, 8, 8),
    # shift 0 0 0 keeps Gamma in the mesh — mandatory for tracing band-edge
    # eigenvalues of a direct-gap zincblende like CdTe/HgTe
    kshift: tuple[int, int, int] = (0, 0, 0),
    kpoints_explicit: np.ndarray | None = None,  # (nk, 4) crystal coords+weight
    occupations: str = "fixed",
    smearing: str | None = None,
    degauss_ry: float | None = None,
    nbnd: int | None = None,
    soc: bool = False,
    conv_thr: float = 1.0e-10,
    mixing_beta: float = 0.3,
    electron_maxstep: int = 200,
    diagonalization: str = "david",
    diago_full_acc: bool = False,
    extra_system: dict | None = None,
    extra_control: dict | None = None,
    relax: bool = False,
    vc_relax: bool = False,
    pseudo_dir: str = "./",
    outdir: str = "./out",
    verbosity: str = "high",
) -> str:
    """Build a complete pw.x input deck as a string."""
    control = {
        "calculation": "vc-relax" if vc_relax else ("relax" if relax else calculation),
        "pseudo_dir": pseudo_dir,
        "outdir": outdir,
        "verbosity": verbosity,
        "tprnfor": relax or vc_relax,
        "tstress": vc_relax,
        "wf_collect": True,
    }
    if extra_control:
        control.update(extra_control)

    system = {
        "ibrav": 0,  # explicit CELL_PARAMETERS, supports arbitrary SQS cells
        "nat": len(atoms),
        "ntyp": len(set(atoms.get_chemical_symbols())),
        "ecutwfc": ecutwfc_ry,
        "occupations": occupations,
    }
    if ecutrho_ry is not None:
        system["ecutrho"] = ecutrho_ry
    if occupations not in ("fixed",):
        system["occupations"] = occupations
        system["smearing"] = smearing or "mv"
        system["degauss"] = degauss_ry if degauss_ry is not None else 0.005
    if nbnd is not None:
        system["nbnd"] = nbnd
    if soc:
        # fully-relativistic US/PAW potentials: spin-orbit needs noncollinear
        system["noncolin"] = True
        system["lspinorb"] = True
    if extra_system:
        for k, v in extra_system.items():
            if v is not None:
                system[k] = v

    electrons = {
        "conv_thr": conv_thr,
        "mixing_beta": mixing_beta,
        "electron_maxstep": electron_maxstep,
        "diagonalization": diagonalization,
        # Needed for US+SOC band runs: avoids 'S matrix not positive definite'
        "diago_full_acc": True if diago_full_acc else None,
    }

    parts = [_fmt_block("CONTROL", control), _fmt_block("SYSTEM", system),
             _fmt_block("ELECTRONS", electrons)]
    if relax or vc_relax:
        parts.append(_fmt_block("IONS", {"ion_dynamics": "bfgs"}))
    if vc_relax:
        parts.append(_fmt_block("CELL", {"cell_dynamics": "bfgs",
                                         "cell_factor": 3.0}))

    symbols = list(dict.fromkeys(atoms.get_chemical_symbols()))
    parts.append("ATOMIC_SPECIES\n" + "\n".join(
        f"  {s:3s}  {atomic_mass(s):.4f}  {pseudopotentials[s]}"
        for s in symbols))
    cell = atoms.cell.array
    lines = ["CELL_PARAMETERS angstrom"]
    lines += ["  " + " ".join(f"{v: .12f}" for v in row) for row in cell]
    parts.append("\n".join(lines))

    pos = atoms.get_scaled_positions()
    syms = atoms.get_chemical_symbols()
    lines = ["ATOMIC_POSITIONS crystal"]
    lines += [f"  {s:3s}  {p[0]: .12f} {p[1]: .12f} {p[2]: .12f}"
              for s, p in zip(syms, pos)]
    parts.append("\n".join(lines))

    if kpoints_explicit is not None:
        kp = np.atleast_2d(kpoints_explicit)
        lines = [f"K_POINTS crystal", f"{len(kp)}"]
        lines += [f"  {row[0]: .12f} {row[1]: .12f} {row[2]: .12f}  {row[3]: .6f}"
                  for row in kp]
        parts.append("\n".join(lines))
    else:
        part = "K_POINTS automatic\n"
        part += f"  {kgrid[0]} {kgrid[1]} {kgrid[2]}  {kshift[0]} {kshift[1]} {kshift[2]}"
        parts.append(part)
    return "\n\n".join(parts) + "\n"


# -- atomic masses (adequate for QE; it uses them only for dynamics/labels) --
_ATOMIC_MASS = {"H": 1.008, "Si": 28.085, "Cd": 112.414, "Hg": 200.592,
                "Te": 127.60}


def atomic_mass(symbol: str) -> float:
    return _ATOMIC_MASS.get(symbol, 100.0)


def bands_x_input(filband: str = "bands.dat", lsym: bool = False,
                  outdir: str = "./out", prefix: str = "pwscf") -> str:
    return _fmt_block("BANDS", {
        "outdir": outdir, "prefix": prefix,
        "filband": filband, "lsym": lsym, "no_overlap": True}) + "\n"


def dos_x_input(fildos: str = "dos.dat", e1: float | None = None,
                e2: float | None = None, de: float = 0.01,
                outdir: str = "./out", prefix: str = "pwscf") -> str:
    p: dict = {"fildos": fildos, "DeltaE": de,
               "outdir": outdir, "prefix": prefix}
    if e1 is not None:
        p["Emin"] = e1
    if e2 is not None:
        p["Emax"] = e2
    return _fmt_block("DOS", p) + "\n"


def projwfc_x_input(filpdos: str = "pdos", degauss: float = 0.01,
                    outdir: str = "./out", prefix: str = "pwscf") -> str:
    return _fmt_block("PROJWFC", {"filpdos": filpdos, "DeltaE": 0.01,
                                  "ngauss": 0, "degauss": degauss,
                                  "outdir": outdir, "prefix": prefix}) + "\n"


def epsilon_x_input(nk1: int, nk2: int, nk3: int,
                    intersmear: float = 0.136,  # eV
                    wmax: float = 30.0,         # eV
                    nw: int = 600,
                    outdir: str = "./out", prefix: str = "pwscf",
                    ) -> str:
    """epsilon.x (IP/RPA dielectric tensor) input deck, QE 7.x layout."""
    return (
        _fmt_block("INPUTPP", {
            "calculation": "eps", "outdir": outdir, "prefix": prefix})
        + "\n"
        + _fmt_block("ENERGY_GRID", {
            "smeartype": "gauss",
            "intersmear": intersmear,
            "intrasmear": 0.0,
            "wmin": 0.0, "wmax": wmax, "nw": nw, "shift": 0.0})
        + "\n"
    )


def write_text(path, text: str, root=None) -> Path:
    p = ensure_within_workspace(path, root)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as fh:
        fh.write(text)
    return p
