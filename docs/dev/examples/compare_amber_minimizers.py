"""Summarise refinement runs of the three AMBER minimiser variants.

Usage: .venv-openff/bin/python docs/dev/examples/compare_amber_minimizers.py <run_dir> [<run_dir> ...]
Each run_dir must contain refine.log and refined_stats.json.
"""
import sys, os, re, json
import numpy


def accepted_ff(log):
    txt = open(log).read()
    blocks = re.split(r"====== CYCLE\s+\d+ ======", txt)
    if len(blocks) < 2:
        return []
    first = [float(l.split("=")[1]) for l in blocks[1].splitlines() if l.startswith(" ff= ")]
    out = [first[0]] if first else []
    for b in blocks[1:]:
        ffs = [float(l.split("=")[1]) for l in b.splitlines() if l.startswith(" ff= ")]
        if ffs:
            out.append(ffs[-1])
    return out


def summarise(d):
    log = os.path.join(d, "refine.log")
    txt = open(log).read()
    st = json.load(open(os.path.join(d, "refined_stats.json")))
    g = st[-1]["geom"]["summary"]["r.m.s.Z"]
    fsc = {m.group(1): float(m.group(2)) for m in re.finditer(r"FSCaverage\((\w+)\)\s*=\s*([-\d.]+)", txt)}
    ff = accepted_ff(log)
    wall = None
    tf = os.path.join(d, "time.log")
    if os.path.exists(tf):
        last = [l.strip() for l in open(tf) if l.strip()]
        if last:
            try: wall = float(last[-1])
            except ValueError: pass
    mini = "lbfgs" if "Minimizer: L-BFGS-B" in txt else "gn"
    offdiag = "hessian off-diagonal: True" in txt
    nit = [int(m.group(1)) for m in re.finditer(r"L-BFGS-B: nit=(\d+)", txt)]
    nfev = [int(m.group(1)) for m in re.finditer(r"nfev=(\d+)", txt)]
    return {"dir": os.path.basename(d), "minimizer": mini, "offdiag": offdiag,
            "ff": ff, "fsc": fsc, "bond": g.get("Bond distances, non H"),
            "angle": g.get("Bond angles, non H"),
            "vdw": next((v for k, v in g.items() if k.lower().startswith("vdw")), float("nan")),
            "not_minimised": txt.count("WARNING: function not minimised"),
            "fval": [c.get("fval") for c in st[1:]],
            "wall": wall, "nit": nit, "nfev": nfev,
            "ff_weight": next((float(m.group(1)) for m in
                               re.finditer(r"ff_weight determined automatically: ([\d.eE+-]+)", txt)), None)}


rows = [summarise(d) for d in sys.argv[1:]]
print("%-16s %-7s %-8s %8s %11s %11s %7s %7s %7s %6s %8s" % (
    "run", "minim.", "offdiag", "w_ff", "E_ff cyc1", "E_ff last", "bondZ", "angleZ", "vdwZ", "notmin", "wall(s)"))
for r in rows:
    print("%-16s %-7s %-8s %8s %11s %11s %7.3f %7.3f %7.3f %6d %8s" % (
        r["dir"], r["minimizer"], r["offdiag"], "%.3f" % r["ff_weight"] if r["ff_weight"] else "-",
        "%.0f" % r["ff"][1] if len(r["ff"]) > 1 else "-", "%.0f" % r["ff"][-1] if r["ff"] else "-",
        r["bond"], r["angle"], r["vdw"], r["not_minimised"],
        "%.0f" % r["wall"] if r["wall"] else "-"))
print()
for r in rows:
    print("%-16s FSC %s" % (r["dir"], " ".join("%s=%.4f" % kv for kv in sorted(r["fsc"].items()))))
    print("%-16s E_ff per cycle: %s" % ("", " ".join("%.0f" % x for x in r["ff"])))
    print("%-16s f per cycle   : %s" % ("", " ".join("%.5e" % x for x in r["fval"] if x is not None)))
    if r["nit"]:
        print("%-16s L-BFGS nit %s  nfev %s" % ("", r["nit"], r["nfev"]))
