"""UPF pseudopotential metadata extraction and SOC validation.

We never trust filenames: relativistic character and spin-orbit content are
read from the UPF header itself (``<PP_HEADER ...>`` attributes), and
production SOC runs are refused unless every selected potential is fully
relativistic and carries spin-orbit projectors.
"""
from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field


@dataclass
class UPFMeta:
    path: str
    element: str = ""
    functional: str = ""
    pseudo_type: str = ""      # NC / US / PAW
    relativistic: str = ""     # 'scalar' or 'full'
    has_so: bool = False
    z_valence: float = float("nan")
    l_max: int = -1
    is_ultrasoft: bool = False
    is_paw: bool = False
    core_correction: bool = False
    suggested_ecutwfc_ry: float | None = None
    suggested_ecutrho_ry: float | None = None
    extras: dict = field(default_factory=dict)

    def soc_capable(self) -> bool:
        return (
            self.relativistic.lower() == "full"
            and self.has_so
        )


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


_ATTR_RE = re.compile(r'(\w[\w-]*)\s*=\s*"([^"]*)"')


def parse_upf_header(path: str, header_lines: int = 400) -> UPFMeta:
    """Parse the PP_HEADER block of a UPF (v1 pseudo-XML or v2 XML) file."""
    meta = UPFMeta(path=str(path))
    with open(path, "r", errors="replace") as fh:
        text = "".join(fh.readline() for _ in range(header_lines))

    m = re.search(r"<PP_HEADER\b([^>]*)>", text, flags=re.IGNORECASE)
    if m is None:
        # UPF v1 plain header: first line carries element/type fragments.
        lines = text.splitlines()
        meta.extras["v1_first_line"] = lines[0] if lines else ""
        if len(lines) > 1:
            meta.extras["v1_second_line"] = lines[1]
        return meta

    attrs = dict(_ATTR_RE.findall(m.group(1)))
    meta.element = attrs.get("element", "").strip()
    meta.functional = attrs.get("functional", attrs.get("xc", "")).strip()
    meta.pseudo_type = attrs.get("pseudo_type", "").strip()
    meta.relativistic = attrs.get("relativistic", "").strip()
    meta.has_so = attrs.get("has_so", "").strip().lower() in ("t", "true", "1")
    meta.is_ultrasoft = attrs.get("is_ultrasoft", "").strip().lower() in (
        "t", "true", "1")
    meta.is_paw = attrs.get("is_paw", "").strip().lower() in ("t", "true", "1")
    meta.core_correction = attrs.get("core_correction", "").strip().lower() in (
        "t", "true", "1")
    try:
        meta.z_valence = float(attrs.get("z_valence", "nan"))
    except ValueError:
        pass
    try:
        meta.l_max = int(attrs.get("l_max", attrs.get("lmax", "-1")))
    except ValueError:
        pass
    # suggested cutoffs (present in some library releases)
    for k in list(attrs):
        kl = k.lower().replace("-", "_")
        if "wfc" in kl and ("cut" in kl or "e" == kl[-2:-1]):
            try:
                meta.suggested_ecutwfc_ry = float(attrs[k])
            except ValueError:
                pass
        if ("rho" in kl and "cut" in kl) or kl == "ecutrho":
            try:
                meta.suggested_ecutrho_ry = float(attrs[k])
            except ValueError:
                pass
    return meta


def validate_soc_set(metas: list[UPFMeta], xc_expected: str = "PBE") -> list[str]:
    """Reject an SOC production run unless every PP is fully relativistic etc.

    Returns a list of problem strings; empty list == OK.
    """
    problems: list[str] = []
    for m in metas:
        tag = f"{m.element or '?'} ({m.path})"
        if m.relativistic.lower() != "full":
            problems.append(
                f"{tag}: relativistic='{m.relativistic or 'unknown'}', "
                "need 'full' for SOC production runs"
            )
        if not m.has_so:
            problems.append(f"{tag}: has_so is false/missing")
        if m.element == "":
            problems.append(f"{tag}: element not identified in PP_HEADER")
        if (
            xc_expected
            and m.functional
            and xc_expected.upper() not in m.functional.upper()
        ):
            problems.append(
                f"{tag}: functional '{m.functional}' does not match "
                f"required '{xc_expected}'"
            )
    # family consistency: all PPs should share pseudo_type
    ptypes = {m.pseudo_type.upper() for m in metas if m.pseudo_type}
    if len(ptypes) > 1:
        problems.append(f"mixed pseudo_type across set: {sorted(ptypes)}")
    return problems
