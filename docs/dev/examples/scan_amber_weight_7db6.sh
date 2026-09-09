#!/usr/bin/env bash
# Scan --amber_weight on 7db6 (MT1-Gi1 + ramelteon) and summarise FSC / E_AMBER / geometry per weight.
# Data are fetched by run_7db6_amber_openff.sh conventions (tests/7db6, MD5 verified there).
#
# Usage: bash docs/dev/examples/scan_amber_weight_7db6.sh [outdir] [ncycle] [parallel] [weights...]
#   EXTRA_ARGS="--cross_validation" to refine against half map 1 and report FSC against half map 2 (overfitting check)
set -euo pipefail
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "$0")/../../.." && pwd)}"
OUTDIR="${1:-$PROJECT_ROOT/work/7db6_weight_scan}"
NCYCLE="${2:-5}"
PAR="${3:-3}"
EXTRA_ARGS="${EXTRA_ARGS:-}"
shift 3 2>/dev/null || true
WEIGHTS=("$@"); [ ${#WEIGHTS[@]} -gt 0 ] || WEIGHTS=(0.02 0.05 0.1 0.2 0.5 1.0 2.0)
PY="$PROJECT_ROOT/.venv-openff/bin/python"
DATA="$PROJECT_ROOT/tests/7db6"
export CLIBD_MON="${CLIBD_MON:-$PROJECT_ROOT/third_party/monomers}"
export PATH="$PROJECT_ROOT/.venv-openff/bin:$PATH"
for f in emd_30627_half_map_1.map.gz emd_30627_half_map_2.map.gz 7db6.cif.gz; do
    [ -f "$DATA/$f" ] || { echo "missing $DATA/$f: run run_7db6_amber_openff.sh first"; exit 1; }
done
PLATFORM=$("$PY" -c "import openmm; print('CPU' if 'CPU' in [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())] else 'Reference')")
mkdir -p "$OUTDIR"

run_one() {  # weight
    local w="$1" d="$OUTDIR/w$1"
    mkdir -p "$d" && cd "$d"
    if [ "$w" = "0" ]; then
        "$PY" -m servalcat refine_spa_norefmac --model "$DATA/7db6.cif.gz" \
            --halfmaps "$DATA/emd_30627_half_map_1.map.gz" "$DATA/emd_30627_half_map_2.map.gz" \
            -d 3.3 --ncycle "$NCYCLE" --hout $EXTRA_ARGS -o refined > refine.log 2>&1
    else
        "$PY" -m servalcat refine_spa_norefmac --model "$DATA/7db6.cif.gz" \
            --halfmaps "$DATA/emd_30627_half_map_1.map.gz" "$DATA/emd_30627_half_map_2.map.gz" \
            -d 3.3 --ncycle "$NCYCLE" --amber_enable --amber_weight "$w" --amber_his_state HIE \
            --amber_ligand_ff openff-2.2.1 --amber_ligand_charge nagl --amber_platform "$PLATFORM" \
            --hout $EXTRA_ARGS -o refined > refine.log 2>&1
    fi
    echo "done w=$w"
}
export -f run_one; export PY DATA NCYCLE PLATFORM OUTDIR EXTRA_ARGS
printf "%s\n" 0 "${WEIGHTS[@]}" | xargs -P "$PAR" -I{} bash -c 'run_one {}'

# ---- summary
"$PY" - "$OUTDIR" 0 "${WEIGHTS[@]}" <<'PYEOF'
import sys, os, json, re, servalcat, gemmi, numpy
outdir, weights = sys.argv[1], sys.argv[2:]
print("%6s %9s %11s %11s %9s %9s %9s %s" % ("w_ff", "FSC(full)", "E_ff cyc1", "E_ff last", "bondZ", "angleZ", "vdwZ", "fval decreased"))
for w in weights:
    d = os.path.join(outdir, "w" + w)
    log = open(os.path.join(d, "refined.log")).read()
    fsc = [float(l.split("=")[1]) for l in log.splitlines() if l.strip().startswith("FSCaverage(full)")]
    other = {}
    for l in log.splitlines():
        m = re.match(r"\s*FSCaverage\((\w+)\)\s*=\s*([-\d.]+)", l)
        if m and m.group(1) != "full":
            other[m.group(1)] = float(m.group(2))
    ff = [float(l.split("=")[1]) for l in log.splitlines() if l.startswith(" ff= ")]
    stats = json.load(open(os.path.join(d, "refined_stats.json")))
    g = stats[-1]["geom"]["summary"]
    bz = g["r.m.s.Z"].get("Bond distances, non H"); az = g["r.m.s.Z"].get("Bond angles, non H")
    vz = next((v for k, v in g["r.m.s.Z"].items() if k.lower().startswith("vdw")), float("nan"))
    dec = [c.get("fval_decreased") for c in stats[1:]]
    e1 = ff[1] if len(ff) > 1 else float("nan"); el = ff[-1] if ff else float("nan")
    extra = " ".join("%s=%.4f" % kv for kv in sorted(other.items()))
    print("%6s %9.4f %11.0f %11.0f %9.3f %9.3f %9.3f %s %s" % (w, fsc[-1], e1, el, bz, az, vz, "".join("T" if x else "F" for x in dec), extra))
PYEOF
