"""Plot the AMBER weight / Hessian study on 7db6 and 7dy0 from the run directories in work/.

Usage: .venv-openff/bin/python docs/dev/examples/plot_amber_weight_study.py [work_dir] [fig_dir]
Produces PNG figures and a JSON table of the plotted numbers.
"""
import sys, os, re, json, glob
import numpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORK = sys.argv[1] if len(sys.argv) > 1 else "work"
FIG = sys.argv[2] if len(sys.argv) > 2 else "docs/dev/figures"
os.makedirs(FIG, exist_ok=True)

# ---- palette (docs/dev dataviz reference, light mode)
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
SEQ = ["#86b6ef", "#6da7ec", "#3987e5", "#2a78d6", "#1c5cab", "#184f95", "#0d366b"]  # ordinal blue ramp 250..700
plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": INK, "text.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURFACE,
                     "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "grid.color": GRID, "grid.linewidth": 0.8,
                     "axes.grid": True, "axes.grid.axis": "y", "axes.axisbelow": True, "lines.linewidth": 2,
                     "lines.markersize": 6, "legend.frameon": False})

# ---- parsers
def accepted_ff(log):
    """E_AMBER at start and accepted value after each cycle."""
    txt = open(log).read()
    blocks = re.split(r"====== CYCLE\s+\d+ ======", txt)
    ffs0 = [float(l.split("=")[1]) for l in blocks[1].splitlines() if l.startswith(" ff= ")] if len(blocks) > 1 else []
    out = [ffs0[0]] if ffs0 else []
    for b in blocks[1:]:
        ffs = [float(l.split("=")[1]) for l in b.splitlines() if l.startswith(" ff= ")]
        if ffs:
            out.append(ffs[-1])
    return out

def fsc(log):
    d = {}
    for l in open(log):
        m = re.match(r"\s*FSCaverage\((\w+)\)\s*=\s*([-\d.]+)", l)
        if m:
            d[m.group(1)] = float(m.group(2))
    return d

def rmsz(stats_json):
    st = json.load(open(stats_json))
    g = st[-1]["geom"]["summary"]["r.m.s.Z"]
    vdw = next((v for k, v in g.items() if k.lower().startswith("vdw")), float("nan"))
    return {"bond": g["Bond distances, non H"], "angle": g["Bond angles, non H"], "vdw": vdw}

def ff_weight(log):
    for l in open(log):
        m = re.search(r"ff_weight determined automatically: ([\d.eE+-]+)", l)
        if m:
            return float(m.group(1))
    return None

def scan(dirname, logname="refined.log"):
    rows = []
    for d in sorted(glob.glob(os.path.join(WORK, dirname, "w*")), key=lambda p: float(os.path.basename(p)[1:])):
        w = float(os.path.basename(d)[1:])
        log = os.path.join(d, logname)
        rows.append({"w": w, "fsc": fsc(log), "ff": accepted_ff(log), "rmsz": rmsz(os.path.join(d, "refined_stats.json"))})
    return rows

data = {"scan": scan("7db6_weight_scan"), "scan_cv": scan("7db6_weight_scan_cv")}
auto = {}
for tag in ("full", "cv"):
    log = os.path.join(WORK, "7db6_weight_auto", tag, "refine.log")
    auto[tag] = {"w": ff_weight(log), "fsc": fsc(log), "ff": accepted_ff(log),
                 "rmsz": rmsz(os.path.join(WORK, "7db6_weight_auto", tag, "refined_stats.json"))}
data["auto"] = auto
dy0 = {}
for label, path in (("const 10, H at electron pos. (old default)", "7dy0_hessian_study/amber5_hd10/run.log"),
                    ("const 1000, H at electron pos.", "7dy0_hessian_study/amber5_hd1000/run.log"),
                    ("const 1000, H at nucleus pos.", "7dy0_hessian_study/patched5_hd1000/run.log"),
                    ("bonded, HIP, w=0.05 (default)", "7dy0_hessian_study/new5_HIP/run.log"),
                    ("bonded, HIE, w=0.05", "7dy0_hessian_study/new5_HIE/run.log"),
                    ("bonded, HIE, w=auto (0.47)", "7dy0_weight_auto/refine.log")):
    p = os.path.join(WORK, path)
    if os.path.exists(p):
        dy0[label] = {"ff": accepted_ff(p), "fsc": fsc(p)}
data["dy0"] = dy0
json.dump(data, open(os.path.join(FIG, "amber_weight_study_data.json"), "w"), indent=1, ensure_ascii=False)

# ---- helpers
def style(ax, title, xlabel=None, ylabel=None):
    ax.set_title(title, loc="left", fontsize=10, pad=8)
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    ax.tick_params(length=0)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(AXIS)

def weight_axis(ax, ws):
    """ordinal positions for weights (0 included) with tick labels"""
    ax.set_xticks(range(len(ws)))
    ax.set_xticklabels(["0\n(off)" if w == 0 else ("%g" % w) for w in ws])
    ax.set_xlim(-0.4, len(ws) - 0.6)

def end_labels(ax, x, items, min_gap, fontsize=8):
    """items: list of (y, text, color). Place labels at x, pushed apart vertically by >= min_gap (data units)."""
    items = sorted(items, key=lambda t: t[0])
    ys = [t[0] for t in items]
    for i in range(1, len(ys)):
        if ys[i] - ys[i-1] < min_gap:
            ys[i] = ys[i-1] + min_gap
    # re-centre so the block stays around the original mean
    shift = (numpy.mean(ys) - numpy.mean([t[0] for t in items]))
    ys = [y - shift for y in ys]
    for (y0, text, color), y in zip(items, ys):
        ax.text(x, y, text, color=color, fontsize=fontsize, va="center")


def interp_pos(w, ws):
    """position of weight w on the ordinal axis by log interpolation between neighbours"""
    pos = [i for i, x in enumerate(ws) if x > 0]; vals = [ws[i] for i in pos]
    return numpy.interp(numpy.log(w), numpy.log(vals), pos)

rows = data["scan"]; ws = [r["w"] for r in rows]
# ================= Fig 1: metrics vs weight =================
fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
ax = axes[0]
y = [r["fsc"]["full"] for r in rows]
ax.plot(range(len(ws)), y, color=CAT[0], marker="o", markeredgecolor=SURFACE, markeredgewidth=1.5)
xa = interp_pos(auto["full"]["w"], ws)
ax.plot([xa], [auto["full"]["fsc"]["full"]], color=CAT[1], marker="D", linestyle="none", markeredgecolor=SURFACE, markeredgewidth=1.5)
ax.annotate("auto (R=0.3)\nw=%.2f" % auto["full"]["w"], (xa, auto["full"]["fsc"]["full"]), textcoords="offset points", xytext=(8, 8), fontsize=8, color=INK2)
ax.set_ylim(0.75, 0.84)
style(ax, "FSC(full) vs AMBER weight", "w_ff", "FSCaverage(full)"); weight_axis(ax, ws)

ax = axes[1]
e1 = [r["ff"][1] if len(r["ff"]) > 1 else numpy.nan for r in rows]
e5 = [r["ff"][-1] if len(r["ff"]) > 1 else numpy.nan for r in rows]
ax.axhline(0, color=AXIS, linewidth=1)
ax.plot(range(len(ws)), numpy.array(e1) / 1000, color=CAT[0], marker="o", markeredgecolor=SURFACE, markeredgewidth=1.5, label="after cycle 1")
ax.plot(range(len(ws)), numpy.array(e5) / 1000, color=CAT[1], marker="o", markeredgecolor=SURFACE, markeredgewidth=1.5, label="after cycle 5")
ax.plot([xa, xa], [auto["full"]["ff"][1] / 1000, auto["full"]["ff"][-1] / 1000], color=CAT[2], marker="D", linestyle="none", markeredgecolor=SURFACE, markeredgewidth=1.5, label="auto (R=0.3)")
ax.legend(loc="upper right", fontsize=8)
style(ax, "E_AMBER vs weight", "w_ff", "E_AMBER (10^3 kJ/mol)"); weight_axis(ax, ws)

ax = axes[2]
for i, (k, lab) in enumerate((("bond", "bond"), ("angle", "angle"), ("vdw", "VDW"))):
    yy = [r["rmsz"][k] for r in rows]
    ax.plot(range(len(ws)), yy, color=CAT[i], marker="o", markeredgecolor=SURFACE, markeredgewidth=1.5, label=lab)
    ax.text(len(ws) - 0.85, yy[-1], lab, color=INK2, fontsize=8, va="center")
ax.set_ylim(0.4, 1.5)
style(ax, "Geometry rmsZ vs weight", "w_ff", "r.m.s. Z after cycle 5"); weight_axis(ax, ws)
fig.tight_layout(w_pad=2)
fig.savefig(os.path.join(FIG, "fig_7db6_weight_scan.png"), dpi=160); plt.close(fig)

# ================= Fig 2: E_AMBER trajectories per weight =================
fig, ax = plt.subplots(figsize=(7, 4))
ax.axhline(0, color=AXIS, linewidth=1)
nz = [r for r in rows if r["w"] > 0]
labels = []
for i, r in enumerate(nz):
    y = numpy.array(r["ff"]) / 1000
    ax.plot(range(len(y)), y, color=SEQ[i], marker="o", markeredgecolor=SURFACE, markeredgewidth=1.5, label="w=%g" % r["w"])
    labels.append((y[-1], "w=%g" % r["w"], INK2))
y = numpy.array(auto["full"]["ff"]) / 1000
ax.plot(range(len(y)), y, color=CAT[1], marker="D", markeredgecolor=SURFACE, markeredgewidth=1.5, label="auto (R=0.3), w=%.2f" % auto["full"]["w"])
labels.append((y[-1], "auto w=%.2f" % auto["full"]["w"], CAT[1]))
end_labels(ax, 5.12, labels, min_gap=4.0)
ax.legend(loc="upper left", fontsize=8, ncol=2)
ax.set_xticks(range(6)); ax.set_xlim(-0.2, 6.4)
style(ax, "7db6: accepted E_AMBER per cycle", "cycle", "E_AMBER (10^3 kJ/mol)")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_7db6_eamber_trajectory.png"), dpi=160); plt.close(fig)

# ================= Fig 3: cross-validation =================
cv = data["scan_cv"]; wcv = [r["w"] for r in cv]
fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.6))
ax = axes[0]
h1 = [r["fsc"]["half1"] for r in cv]; h2 = [r["fsc"]["half2"] for r in cv]
ax.plot(range(len(wcv)), h1, color=CAT[0], marker="o", markeredgecolor=SURFACE, markeredgewidth=1.5, label="work (half1)")
ax.plot(range(len(wcv)), h2, color=CAT[1], marker="o", markeredgecolor=SURFACE, markeredgewidth=1.5, label="free (half2)")
xa = interp_pos(auto["cv"]["w"], wcv)
ax.plot([xa, xa], [auto["cv"]["fsc"]["half1"], auto["cv"]["fsc"]["half2"]], color=CAT[2], marker="D", linestyle="none", markeredgecolor=SURFACE, markeredgewidth=1.5, label="auto (R=0.3)")
ax.text(len(wcv) - 0.85, h1[-1], "work", color=INK2, fontsize=8, va="center"); ax.text(len(wcv) - 0.85, h2[-1], "free", color=INK2, fontsize=8, va="center")
ax.legend(loc="lower left", fontsize=8); ax.set_ylim(0.69, 0.80)
style(ax, "Cross-validation FSC (refined against half1)", "w_ff", "FSCaverage"); weight_axis(ax, wcv)
ax = axes[1]
gap = numpy.array(h1) - numpy.array(h2)
ax.plot(range(len(wcv)), gap, color=CAT[0], marker="o", markeredgecolor=SURFACE, markeredgewidth=1.5)
ax.plot([xa], [auto["cv"]["fsc"]["half1"] - auto["cv"]["fsc"]["half2"]], color=CAT[2], marker="D", linestyle="none", markeredgecolor=SURFACE, markeredgewidth=1.5)
ax.set_ylim(0, 0.08)
style(ax, "work - free (overfitting indicator)", "w_ff", "delta FSC"); weight_axis(ax, wcv)
fig.tight_layout(w_pad=2); fig.savefig(os.path.join(FIG, "fig_7db6_cv.png"), dpi=160); plt.close(fig)

# ================= Fig 4: 7dy0 Hessian / hydrogen / weight study =================
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1, 1.25]})
ax = axes[0]
ax.axhline(0, color=AXIS, linewidth=1)
for i, (label, d) in enumerate(dy0.items()):
    y = numpy.array(d["ff"]) / 1000
    ax.plot(range(len(y)), y, color=CAT[i], marker="o", markeredgecolor=SURFACE, markeredgewidth=1.5, label=label)
ax.set_xticks(range(6)); ax.set_xlim(-0.2, 5.4)
ax.legend(loc="upper right", fontsize=7.5)
style(ax, "7dy0: all variants (overview)", "cycle", "E_AMBER (10^3 kJ/mol)")
ax = axes[1]
labels = []
for i, (label, d) in enumerate(dy0.items()):
    y = numpy.array(d["ff"]) / 1000
    if y.max() > 0:
        continue  # the diverging old default is shown in the overview only
    ax.plot(range(len(y)), y, color=CAT[i], marker="o", markeredgecolor=SURFACE, markeredgewidth=1.5)
    labels.append((y[-1], label, INK2))
end_labels(ax, 5.12, labels, min_gap=0.32, fontsize=7.5)
ax.set_xticks(range(6)); ax.set_xlim(-0.2, 8.8); ax.set_ylim(-9.0, -4.0)
style(ax, "7dy0: stable variants (zoom)", "cycle", "E_AMBER (10^3 kJ/mol)")
fig.tight_layout(w_pad=2); fig.savefig(os.path.join(FIG, "fig_7dy0_hessian_study.png"), dpi=160); plt.close(fig)
# ================= Fig 5: replacement-version weight scan =================
def scan_rep(dirname):
    rows = []
    for d in sorted(glob.glob(os.path.join(WORK, dirname, "w*")), key=lambda p: float(os.path.basename(p)[1:])):
        w = float(os.path.basename(d)[1:])
        log = os.path.join(d, "refined.log")
        st = json.load(open(os.path.join(d, "refined_stats.json")))
        z = st[-1]["geom"]["summary"]["r.m.s.Z"]
        rows.append({"w": w, "fsc": fsc(log), "ff": accepted_ff(log),
                     "bond": z.get("Bond distances, non H", numpy.nan),
                     "angle": z.get("Bond angles, non H", numpy.nan),
                     "chir": next((v for k, v in z.items() if k.lower().startswith("chir")), numpy.nan)})
    return rows

if os.path.isdir(os.path.join(WORK, "7db6_weight_scan_replace")):
    rep = scan_rep("7db6_weight_scan_replace")
    repcv = scan_rep("7db6_weight_scan_replace_cv")
    data["scan_replace"] = rep
    data["scan_replace_cv"] = repcv
    json.dump(data, open(os.path.join(FIG, "amber_weight_study_data.json"), "w"), indent=1, ensure_ascii=False)
    ws = [r["w"] for r in rep]
    cvw = [r["w"] for r in repcv]
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.7))
    ax = axes[0]
    ax.plot(range(len(ws)), [r["fsc"]["full"] for r in rep], color=CAT[0], marker="o",
            markeredgecolor=SURFACE, markeredgewidth=1.5, label="FSC(full)")
    pos = {w: i for i, w in enumerate(ws)}
    xs = [pos[r["w"]] for r in repcv if r["w"] in pos]
    ax.plot(xs, [r["fsc"]["half2"] for r in repcv if r["w"] in pos], color=CAT[1], marker="D",
            markeredgecolor=SURFACE, markeredgewidth=1.5, label="free FSC(half2)")
    ax.legend(loc="upper right", fontsize=8)
    style(ax, "Full replacement: FSC vs w_ff", "w_ff", "FSCaverage"); weight_axis(ax, ws)
    ax = axes[1]
    ax.plot(range(len(ws)), [r["ff"][1]/1000 if len(r["ff"]) > 1 else numpy.nan for r in rep], color=CAT[0],
            marker="o", markeredgecolor=SURFACE, markeredgewidth=1.5, label="after cycle 1")
    ax.plot(range(len(ws)), [r["ff"][-1]/1000 if r["ff"] else numpy.nan for r in rep], color=CAT[1],
            marker="o", markeredgecolor=SURFACE, markeredgewidth=1.5, label="after cycle 5")
    ax.legend(loc="lower left", fontsize=8)
    style(ax, "Full replacement: E_AMBER vs w_ff", "w_ff", "E_AMBER (10^3 kJ/mol)"); weight_axis(ax, ws)
    ax = axes[2]
    for i, (k, lab) in enumerate((("bond", "bond"), ("angle", "angle"), ("chir", "chirality"))):
        yy = [r[k] for r in rep]
        ax.plot(range(len(ws)), yy, color=CAT[i], marker="o", markeredgecolor=SURFACE, markeredgewidth=1.5, label=lab)
    ax.axhspan(0.4, 1.4, color=CAT[0], alpha=0.08, zorder=0)
    ax.text(len(ws)-1.05, 1.45, "AMBER vs dictionary\nintrinsic offset", fontsize=7, color=MUTED, ha="right", va="bottom")
    ax.legend(loc="upper right", fontsize=8)
    style(ax, "Full replacement: rmsZ vs dictionary", "w_ff", "r.m.s. Z after cycle 5"); weight_axis(ax, ws)
    fig.tight_layout(w_pad=2)
    fig.savefig(os.path.join(FIG, "fig_7db6_weight_scan_replace.png"), dpi=160); plt.close(fig)
    print("  replacement scan figure written")

print("figures written to", FIG)
for k, v in dy0.items(): print("  7dy0", k, "->", ["%.0f" % x for x in v["ff"]], v["fsc"].get("full"))
print("  auto full w=%.4f fsc=%s ff=%s" % (auto["full"]["w"], auto["full"]["fsc"], ["%.0f" % x for x in auto["full"]["ff"]]))
