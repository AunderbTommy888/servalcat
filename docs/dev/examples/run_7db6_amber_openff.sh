#!/usr/bin/env bash
# Example: AMBER + OpenFF prior refinement of 7db6 (melatonin receptor MT1-Gi1 with ramelteon, EMD-30627, 3.3 A)
#
# Requirements (see docs/dev/build_uv_memo_ja.md):
#   - OpenFF environment at $PROJECT_ROOT/.venv-openff (openmm, openmmforcefields, openff-toolkit, rdkit, openff-nagl)
#   - monomer library at $PROJECT_ROOT/third_party/monomers
#   - network access for the first run (about 18 MB of half maps are cached in tests/7db6)
#
# Usage:
#   cd "$PROJECT_ROOT"
#   bash docs/dev/examples/run_7db6_amber_openff.sh [outdir] [ncycle]
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "$0")/../../.." && pwd)}"
OUTDIR="${1:-$PROJECT_ROOT/work/7db6_amber_openff}"
NCYCLE="${2:-5}"
PY="$PROJECT_ROOT/.venv-openff/bin/python"
DATA="$PROJECT_ROOT/tests/7db6"
export CLIBD_MON="${CLIBD_MON:-$PROJECT_ROOT/third_party/monomers}"
export PATH="$PROJECT_ROOT/.venv-openff/bin:$PATH"   # sqm/antechamber for --amber_ligand_charge am1bcc

# ---- 1) data (same files and MD5 sums as tests/test_amber_phase1.py::TestAmberLigandRealData)
mkdir -p "$DATA"
fetch() {  # url md5
    local f="$DATA/$(basename "$1")"
    [ -f "$f" ] || curl -sSL -o "$f" "$1"
    echo "$2  $f" | md5sum -c --quiet
}
fetch https://files.wwpdb.org/pub/emdb/structures/EMD-30627/other/emd_30627_half_map_1.map.gz 319dc338b6a44f9be02f5c4b2c07c1b4
fetch https://files.wwpdb.org/pub/emdb/structures/EMD-30627/other/emd_30627_half_map_2.map.gz 322e28cb5fd70bc893c01c02fbec6a2b
fetch https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/db/7db6.cif.gz 86fce1c3112c758ff4acc6c0f4bc0e54

# ---- 2) platform: CPU if OpenMM provides it, otherwise Reference
PLATFORM=$("$PY" -c "import openmm; print('CPU' if 'CPU' in [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())] else 'Reference')")

# ---- 3) refinement with AMBER prior; ligand JEV is parameterised automatically with OpenFF from the monomer library
#         (replace '--amber_weight 0.1' by '--amber_weight_auto' to set the weight from the gradient norm ratio; see build memo)
mkdir -p "$OUTDIR" && cd "$OUTDIR"
"$PY" -m servalcat refine_spa_norefmac \
    --model "$DATA/7db6.cif.gz" \
    --halfmaps "$DATA/emd_30627_half_map_1.map.gz" "$DATA/emd_30627_half_map_2.map.gz" \
    -d 3.3 \
    --ncycle "$NCYCLE" \
    --amber_enable \
    --amber_weight 0.1 \
    --amber_his_state HIE \
    --amber_ligand_ff openff-2.2.1 \
    --amber_ligand_charge nagl \
    --amber_platform "$PLATFORM" \
    --hout \
    -o refined_amber > refine_amber.log 2>&1

# ---- 4) reference run without the AMBER prior (same data, same cycles)
"$PY" -m servalcat refine_spa_norefmac \
    --model "$DATA/7db6.cif.gz" \
    --halfmaps "$DATA/emd_30627_half_map_1.map.gz" "$DATA/emd_30627_half_map_2.map.gz" \
    -d 3.3 \
    --ncycle "$NCYCLE" \
    --hout \
    -o refined_noamber > refine_noamber.log 2>&1

# ---- 5) summary
echo "== AMBER prior setup (refined_amber.log)"
grep -E "AMBER prior:| ligand " refined_amber.log
echo
echo "== E_AMBER per evaluation (kJ/mol)"
grep -E "^ ff= " refined_amber.log | awk '{printf "%s ", $2} END {print ""}'
echo
echo "== FSCaverage(full)  amber / no-amber"
echo "$(grep 'FSCaverage(full)' refined_amber.log | tail -1 | awk '{print $NF}')  /  $(grep 'FSCaverage(full)' refined_noamber.log | tail -1 | awk '{print $NF}')"
echo
echo "== geometry: bond rmsZ per cycle (amber / no-amber) and ligand contacts"
"$PY" - <<'PYEOF'
import json, servalcat, gemmi
for tag in ("refined_amber", "refined_noamber"):
    st = json.load(open(f"{tag}_stats.json"))
    z = [c["geom"]["summary"]["r.m.s.Z"]["Bond distances, non H"] if isinstance(c["geom"]["summary"], dict) else None for c in st]
    fv = [c.get("fval") for c in st]
    print(f"{tag:16s} bond rmsZ: {['%.3f' % v for v in z if v is not None]}  fval decreased: {[c.get('fval_decreased') for c in st[1:]]}")
    s = gemmi.read_structure(f"{tag}.mmcif")
    jev = [r for ch in s[0] for r in ch if r.name == "JEV"][0]
    ns = gemmi.NeighborSearch(s[0], s.cell, 5).populate()
    close = set()
    for a in jev:
        for m in ns.find_atoms(a.pos, '\0', radius=4.0):
            cra = m.to_cra(s[0])
            if cra.residue.name != "JEV" and cra.atom.element.name != "H":
                close.add(f"{cra.chain.name}/{cra.residue.name}{cra.residue.seqid.num}")
    print(f"{'':16s} JEV heavy-atom contacts within 4 A: {len(close)} residues")
PYEOF
echo
echo "outputs in $OUTDIR: refined_amber.mmcif refined_amber_maps.mtz refined_amber_normalized_fofc.mrc (and refined_noamber.*)"
