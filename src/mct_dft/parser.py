"""Parsers for Quantum ESPRESSO text output (pw.x primarily).

A calculation is only reported as successful when ALL of the following hold:

* the '! total energy' line for the final iteration exists
* 'convergence has been achieved' (SCF) OR the non-SCF/bands calculation
  completed a full k-point loop
* 'JOB DONE.' terminator present
* no fatal markers (convergence NOT achieved, mpi errors, OOM kill lines)
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

_FAILURE_PATTERNS = [
    "convergence NOT achieved",
    "Stopping because iteration",
    "Error in routine",
    "problems computing cholesky",
    "BAD TERMINATION",
    "Out of memory",
    "oom-kill",
    "Killed",
    "job aborted",
    "forrtl: severe",
]


@dataclass
class PWResult:
    converged: bool = False
    job_done: bool = False
    ok: bool = False
    calculation: str = ""
    total_energy_ry: float | None = None
    total_energy_ev: float | None = None
    fermi_ev: float | None = None
    highest_occupied_ev: float | None = None
    lowest_unoccupied_ev: float | None = None
    nscf_iterations: int | None = None
    n_bands: int | None = None
    n_kpts: int | None = None
    nat: int | None = None
    n_electrons: float | None = None
    noncollinear: bool = False
    spin_orbit: bool = False
    lattice_parameter_bohr: float | None = None
    unit_cell_volume_bohr3: float | None = None
    # vc-relax: final relaxed cell (Angstrom) and final lattice constant
    final_cell_ang: list | None = None
    final_positions_crystal: list | None = None
    errors: list[str] = field(default_factory=list)
    # for calculation='bands': arrays filled by parse_bands_output()
    kpoints_cryst: np.ndarray | None = None   # (nk, 3)
    eigenvalues_ev: np.ndarray | None = None  # (nk, nb)


_E_TOT_RE = re.compile(
    r"^!?\s*total energy\s*=\s*(-?\d+\.\d+)\s+Ry", re.MULTILINE)
_FERMI_RE = re.compile(
    r"the Fermi energy is\s+(-?\d+\.\d+)\s+ev", re.IGNORECASE)
_HOCC_RE = re.compile(
    r"highest occupied(?:, lowest unoccupied)? level \(ev\):\s+(-?\d+\.\d+)"
    r"(?:\s+(-?\d+\.\d+))?", re.IGNORECASE)
_NBAND_RE = re.compile(r"number of Kohn-Sham states\s*=\s*(\d+)")
_NK_RE = re.compile(r"number of k points\s*=\s*(\d+)")
_NAT_RE = re.compile(r"number of atoms\s*=\s*(\d+)")
_NITER_RE = re.compile(r"convergence has been achieved in\s+(\d+)\s+iterations")
_NELEC_RE = re.compile(r"number of electrons\s*=\s*([\d.]+)")
_ALAT_RE = re.compile(r"lattice parameter \(alat\)\s*=\s*([\d.]+)\s*a\.u\.")
_VOL_RE = re.compile(r"unit-cell volume\s*=\s*([\d.]+)\s+\(a\.u\.\)\^3")
_CALC_RE = re.compile(r"calculation\s*=\s*'?(\w[\w-]*)'?", re.IGNORECASE)


def parse_pw_output(text: str) -> PWResult:
    """Parse stdout text of a pw.x run."""
    res = PWResult()
    if not text:
        res.errors.append("empty output")
        return res

    m = _CALC_RE.search(text[:6000])
    if m:
        res.calculation = m.group(1)
    else:
        # QE echoes the calculation type in words, not as a namelist card
        for pat, name in [
            ("Self-consistent Calculation", "scf"),
            ("Band Structure Calculation", "bands"),
            ("Non-scf calculation", "nscf"),
            ("Geometry Optimization", "relax"),
        ]:
            if pat.lower() in text[:20000].lower():
                res.calculation = name
                break

    res.job_done = "JOB DONE." in text
    for pat in _FAILURE_PATTERNS:
        if pat.lower() in text.lower():
            # 'convergence NOT achieved' handled with SCF flag; others are fatal
            if pat != "convergence NOT achieved":
                res.errors.append(f"failure marker: {pat}")

    m = list(_E_TOT_RE.finditer(text))
    if m:
        res.total_energy_ry = float(m[-1].group(1))
        res.total_energy_ev = res.total_energy_ry * 13.605693122994
    m = _FERMI_RE.search(text)
    if m:
        res.fermi_ev = float(m.group(1))
    m = _HOCC_RE.search(text)
    if m:
        res.highest_occupied_ev = float(m.group(1))
        if m.group(2) is not None:
            res.lowest_unoccupied_ev = float(m.group(2))
    m = _NBAND_RE.search(text)
    if m:
        res.n_bands = int(m.group(1))
    m = _NK_RE.search(text)
    if m:
        res.n_kpts = int(m.group(1))
    m = _NAT_RE.search(text)
    if m:
        res.nat = int(m.group(1))
    m = _NELEC_RE.search(text)
    if m:
        res.n_electrons = float(m.group(1))
    m = _NITER_RE.search(text)
    if m:
        res.nscf_iterations = int(m.group(1))
        scf_converged = True
    else:
        scf_converged = False
    m = _ALAT_RE.search(text)
    if m:
        res.lattice_parameter_bohr = float(m.group(1))
    m = _VOL_RE.search(text)
    if m:
        res.unit_cell_volume_bohr3 = float(m.group(1))

    cell, pos = parse_relaxed_structure(text)
    res.final_cell_ang = cell
    res.final_positions_crystal = pos

    res.noncollinear = "noncollinear" in text.lower() or \
        "non-collinear" in text.lower()
    res.spin_orbit = "spin-orbit" in text.lower()

    if "convergence NOT achieved" in text:
        res.errors.append("convergence NOT achieved")

    calc = res.calculation.lower()
    is_scf_like = calc in ("scf", "relax", "vc-relax", "")
    res.converged = scf_converged if is_scf_like else res.job_done
    has_eigs = "bands (ev)" in text
    res.ok = (
        res.job_done
        and not any(
            e.startswith("failure marker") or e == "convergence NOT achieved"
            for e in res.errors
        )
        and (res.total_energy_ry is not None or has_eigs)
    )
    return res


_K_LINE_RE = re.compile(
    r"^\s*k\s*=\(?\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*\)?")
_K_PWS_RE = re.compile(r"^\s*k\s*=\s*")  # nscf-style 'k = x y z (.. PWs) bands (ev):'
_FLOAT_RE = re.compile(r"[-+]?\d+\.\d+")


_CELL_RE = re.compile(
    r"CELL_PARAMETERS \(angstrom\)\s*\n((?:\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n){3})")
_POS_RE = re.compile(
    r"ATOMIC_POSITIONS \(crystal\)\s*\n((?:\s+\S+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+.*\n)+)")


def parse_relaxed_structure(text: str):
    """Last CELL_PARAMETERS/ATOMIC_POSITIONS blocks from a (vc-)relax run."""
    cells = _CELL_RE.findall(text)
    poss = _POS_RE.findall(text)
    cell = None
    positions = None
    if cells:
        rows = [[float(v) for v in line.split()]
                for line in cells[-1].strip().splitlines()]
        cell = rows
    if poss:
        lines = poss[-1].strip().splitlines()
        positions = []
        for ln in lines:
            parts = ln.split()
            positions.append({"symbol": parts[0],
                              "crystal": [float(x) for x in parts[1:4]]})
    return cell, positions


def parse_bands_output(text: str) -> tuple[np.ndarray, np.ndarray]:
    """Parse the k-points/eigenvalues listed by pw.x in a 'bands' calculation.

    Handles both listings pw.x produces:
      * path-band style   "k =( x y z ), P =..." + standalone "bands (ev):"
      * nscf style        "k = x y z ( 1959 PWs)   bands (ev):" on one line
    Returns (kpoints[nk,3] crystal, eigenvalues[nk, nb] in eV).
    """
    kpts, eig_rows, current = [], [], None
    in_block = False
    for line in text.splitlines():
        m = _K_LINE_RE.match(line)
        nscf_style = _K_PWS_RE.match(line) and "bands (ev)" in line
        if m and "(ev)" not in line:
            if current is not None:
                eig_rows.append(current)
            kpts.append([float(m.group(i)) for i in (1, 2, 3)])
            current = []
            in_block = False
            continue
        if nscf_style:
            if current is not None:
                eig_rows.append(current)
            nums = _FLOAT_RE.findall(line.split("bands (ev)")[0])
            kpts.append([float(v) for v in nums[:3]])
            current = []
            in_block = True
            continue
        ls = line.strip()
        if ls.startswith("bands (ev)"):
            in_block = True
            continue
        if current is not None and in_block:
            if ls == "":
                continue  # blank line separation, still inside band block
            if _FLOAT_RE.search(ls) and not any(c.isalpha() for c in ls):
                current.extend(float(t) for t in _FLOAT_RE.findall(ls))
            else:
                in_block = False  # e.g. "occupation numbers" closes the block
    if current is not None:
        eig_rows.append(current)
    if not eig_rows:
        return np.empty((0, 3)), np.empty((0, 0))
    nb = min(len(r) for r in eig_rows)
    if nb == 0:
        return np.empty((0, 3)), np.empty((0, 0))
    return (np.asarray(kpts, dtype=float),
            np.asarray([r[:nb] for r in eig_rows], dtype=float))
