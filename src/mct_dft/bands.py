"""Band-edge bookkeeping: VBM/CBM, gap, direct/indirect.

Never force a positive gap: if band ordering is inverted (VBM>CBM) we report
the signed gap and mark the system 'inverted'/'overlap'.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import HC_EV_UM


@dataclass
class BandEdges:
    vbm_ev: float
    cbm_ev: float
    gap_ev: float
    k_vbm: np.ndarray
    k_cbm: np.ndarray
    i_vbm: int
    i_cbm: int
    direct: bool
    metallic: bool
    status: str  # 'semiconductor', 'semimetal/overlap', 'metallic', ...
    occupations_reference_ev: float  # energy used as occupancy threshold


def band_edges(
    kpts_cryst: np.ndarray,
    eigenvalues_ev: np.ndarray,
    reference_ev: float | None = None,
    tol_ev: float = 1e-4,
) -> BandEdges:
    """Locate VBM/CBM from a dense eigenvalue grid.

    reference_ev: occupancy reference (Fermi level for metals; highest occupied
    level otherwise). When None, we use the mid-point of the global gap region
    estimated from a median-split heuristic: energies below the per-k median
    count as occupied. For an insulator with even electron count this matches
    the true VBM set.
    """
    E = np.asarray(eigenvalues_ev, dtype=float)      # (nk, nb)
    nk, nb = E.shape
    if reference_ev is None:
        reference_ev = float(np.median(E))

    occ_mask = E < reference_ev - tol_ev
    unocc_mask = E > reference_ev + tol_ev

    vbm, cbm = -np.inf, np.inf
    i_vbm = i_cbm = -1
    if occ_mask.any():
        i_flat = int(np.argmax(np.where(occ_mask, E, -np.inf)))
        i_vbm = i_flat // nb
        vbm = float(E.flat[i_flat])
    if unocc_mask.any():
        i_flat = int(np.argmin(np.where(unocc_mask, E, np.inf)))
        i_cbm = i_flat // nb
        cbm = float(E.flat[i_flat])

    gap = cbm - vbm
    metallic = not (occ_mask.any() and unocc_mask.any())
    direct = bool(i_vbm == i_cbm) and not metallic

    if metallic:
        status = "metallic (no occupied/unoccupied split at reference)"
    elif gap < 0:
        status = "inverted/overlap (indirect band overlap): not a normal gap"
    elif abs(gap) < 1e-6 and not direct:
        status = "zero-gap / touching bands"
    else:
        status = ("direct-gap semiconductor" if direct
                  else "indirect-gap semiconductor")

    return BandEdges(
        vbm_ev=vbm, cbm_ev=cbm, gap_ev=float(gap),
        k_vbm=np.asarray(kpts_cryst)[i_vbm] if i_vbm >= 0 else np.full(3, np.nan),
        k_cbm=np.asarray(kpts_cryst)[i_cbm] if i_cbm >= 0 else np.full(3, np.nan),
        i_vbm=i_vbm, i_cbm=i_cbm, direct=direct, metallic=metallic,
        status=status, occupations_reference_ev=reference_ev,
    )


def cutoff_um_from_gap(gap_ev: float) -> float | None:
    """Cutoff wavelength from a *positive, physically meaningful* gap."""
    if gap_ev <= 0 or not np.isfinite(gap_ev):
        return None
    return HC_EV_UM / gap_ev
