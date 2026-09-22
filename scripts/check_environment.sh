#!/usr/bin/env bash
# Read-only environment audit. Writes results/provenance/environment.json and
# reports/environment.md. Never modifies anything outside the repository.
set -uo pipefail
cd "$(dirname "$0")/.."
REPO_ROOT="$PWD"
source scripts/env.sh

PY=.venv/bin/python
QE_BIN="$REPO_ROOT/.local/qe-env/bin"

qe_banner() { # $1 = executable name -> first banner line with version
  local exe="$QE_BIN/$1"
  if [ -x "$exe" ]; then
    (echo "" | timeout 10 "$exe" 2>&1 | grep -m1 -oE "Program [A-Z0-9._]+ v\.[0-9.]+" ) || echo "present (no banner)"
  else
    echo "MISSING"
  fi
}

PW_VER=$(qe_banner pw.x)
BANDS_VER=$(qe_banner bands.x)
DOS_VER=$(qe_banner dos.x)
PROJWFC_VER=$(qe_banner projwfc.x)
EPSILON_VER=$(qe_banner epsilon.x)

echo "=== environment audit $(date -u +%FT%TZ) ==="
echo "python : $($PY --version 2>&1)"
for pkg in numpy scipy matplotlib ase icet spglib seekpath pyyaml pytest; do
  v=$($PY -c "import importlib.metadata as m;print(m.version('$pkg'))" 2>/dev/null || echo MISSING)
  echo "$pkg : $v"
done
echo "pw.x      : $PW_VER"
echo "bands.x   : $BANDS_VER"
echo "dos.x     : $DOS_VER"
echo "projwfc.x : $PROJWFC_VER"
echo "epsilon.x : $EPSILON_VER"
echo "mpirun    : $(command -v mpirun >/dev/null && mpirun --version 2>&1 | head -1 || echo MISSING)"
echo "cpus      : $(nproc) logical ($(lscpu 2>/dev/null | awk -F: '/^Core/{gsub(/ /,"",$2);print $2" cores/socket"}'))"
echo "model     : $(lscpu 2>/dev/null | awk -F: '/Model name/{sub(/^ +/,"",$2);print $2}')"
echo "mem       : $(free -h | awk '/^Mem/{print $2" total, "$7" available"}')"
echo "disk      : $(df -h "$REPO_ROOT" | awk 'NR==2{print $4" free of "$2}')"

$PY - <<'PYEOF'
import json, os, shutil, subprocess, sys
from pathlib import Path
from importlib import metadata as md

root = Path(os.environ["WORKSPACE_ROOT"])
qe = root / ".local/qe-env/bin"

def banner(name):
    exe = qe / name
    if not exe.exists():
        return None
    try:
        p = subprocess.run([str(exe)], input="", capture_output=True,
                           text=True, timeout=15)
        for line in (p.stdout + p.stderr).splitlines():
            if "v." in line and "Program" in line:
                return line.strip()
    except Exception:
        pass
    return "present"

pkgs = {}
for pkg in ["numpy","scipy","matplotlib","ase","icet","spglib","seekpath","pyyaml","pytest"]:
    try: pkgs[pkg] = md.version(pkg)
    except Exception: pkgs[pkg] = None

def free_gb(path):
    s = os.statvfs(path)
    return round(s.f_bavail * s.f_frsize / 2**30, 1)

mem = {}
with open("/proc/meminfo") as fh:
    for line in fh:
        if line.startswith(("MemTotal", "MemAvailable")):
            mem[line.split(":")[0]] = int(line.split()[1]) // 1024  # MB

rec = {
    "git_commit": subprocess.run(["git","rev-parse","HEAD"], cwd=root,
                                 capture_output=True, text=True).stdout.strip() or None,
    "python": sys.version.split()[0],
    "python_packages": pkgs,
    "qe": {e: banner(e + ".x") for e in ["pw","bands","dos","projwfc","epsilon"]},
    "qe_bin_dir": str(qe),
    "mpi": shutil.which("mpirun"),
    "cpu_logical": os.cpu_count(),
    "mem_mb": mem,
    "disk_free_gb": free_gb(root),
}
out = root / "results/provenance/environment.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(rec, indent=2))
print(f"[audit] wrote {out}")
PYEOF
echo "[audit] done"
