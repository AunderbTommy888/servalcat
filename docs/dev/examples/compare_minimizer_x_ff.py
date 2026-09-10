"""Summarise the 2x2 comparison of minimiser (Gauss-Newton / L-BFGS) x force field (off / AMBER).

Usage: .venv-openff/bin/python docs/dev/examples/compare_minimizer_x_ff.py <dir containing the run dirs>
Run directories are expected as {dy0,db6}_{gn,lbfgs}_{noff,ff} (see work/minimizer_x_ff/).
"""
import sys, os, re, json, numpy
def rows(base, tags):
    out = []
    for tag, lab in tags:
        d = os.path.join(base, tag)
        if not os.path.exists(os.path.join(d, "refined_stats.json")):
            out.append((lab, None)); continue
        txt = open(os.path.join(d, "refine.log")).read()
        st = json.load(open(os.path.join(d, "refined_stats.json")))
        z = st[-1]["geom"]["summary"]["r.m.s.Z"]
        g = lambda k: next((v for kk, v in z.items() if kk.lower().startswith(k)), float("nan"))
        f = {m.group(1): float(m.group(2)) for m in re.finditer(r"FSCaverage\((\w+)\)\s*=\s*([-\d.]+)", txt)}
        ff = [c.get("ff") for c in st[1:] if c.get("ff") is not None]
        wall = None
        tf = os.path.join(d, "time.log")
        if os.path.exists(tf):
            v = [l.strip() for l in open(tf) if l.strip()]
            try: wall = float(v[-1])
            except Exception: pass
        nit = [int(m.group(1)) for m in re.finditer(r"L-BFGS-B: nit=(\d+)", txt)]
        out.append((lab, dict(full=f.get("full"), work=f.get("half1"), free=f.get("half2"),
                              bond=z.get("Bond distances, non H"), angle=z.get("Bond angles, non H"),
                              vdw=g("vdw"), chir=g("chir"), ff=ff[-1] if ff else None,
                              notmin=txt.count("WARNING: function not minimised"), wall=wall, nit=nit)))
    return out
S = sys.argv[1]
for sysname, tags in (("7dy0 (3.1 A, apo)", [("dy0_gn_noff", "GN,    no FF"), ("dy0_gn_ff", "GN,    AMBER"),
                                             ("dy0_lbfgs_noff", "LBFGS, no FF"), ("dy0_lbfgs_ff", "LBFGS, AMBER")]),
                      ("7db6 (3.3 A, +ligand)", [("db6_gn_noff", "GN,    no FF"), ("db6_gn_ff", "GN,    AMBER"),
                                                 ("db6_lbfgs_noff", "LBFGS, no FF"), ("db6_lbfgs_ff", "LBFGS, AMBER")])):
    print("=== %s ===" % sysname)
    print("  %-13s %9s %8s %8s %8s %8s %8s %8s %10s %6s %7s" % (
        "setting", "FSC(full)", "work", "free", "w-f", "bondZ", "angleZ", "chirZ", "E_ff", "notmin", "wall"))
    for lab, r in rows(S, tags):
        if r is None:
            print("  %-13s (missing)" % lab); continue
        wf = (r["work"] - r["free"]) if (r["work"] and r["free"]) else float("nan")
        print("  %-13s %9.4f %8.4f %8.4f %8.4f %8.3f %8.3f %8.3f %10s %6d %6.0fs%s" % (
            lab, r["full"] or float("nan"), r["work"] or float("nan"), r["free"] or float("nan"), wf,
            r["bond"] if r["bond"] is not None else float("nan"),
            r["angle"] if r["angle"] is not None else float("nan"), r["chir"],
            "%.0f" % r["ff"] if r["ff"] else "-", r["notmin"], r["wall"] or 0,
            "  nit=%s" % r["nit"] if r["nit"] else ""))
