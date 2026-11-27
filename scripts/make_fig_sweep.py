"""Fig. 3: listing time and peak memory against cell size, binary fcc, all compositions.

    python3 scripts/make_fig_sweep.py scripts/sweep_binary_fcc_M2.txt \
        scripts/sweep_binary_fcc_enumlibjl050_M2.txt scripts/sweep_binary_fcc_2628_M2.txt \
        scripts/sweep_binary_fcc_3032_M2.txt -o fig_sweep.pdf

Rows from the first log whose code is enumlibjl are skipped: that run used Enumlib.jl 0.4.2.
"""
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("main"); ap.add_argument("more", nargs="+", help="further sweep logs, merged in order")
ap.add_argument("-o", default="fig_sweep.pdf")
a = ap.parse_args()

def read(path, skip=()):
    d = {}
    for line in open(path):
        if line.startswith("#") or line.startswith("code"):
            continue
        p = line.split()
        if len(p) < 6 or p[5] != "ok" or p[0] in skip:
            continue
        d.setdefault(p[0], {})[int(p[1])] = (float(p[3]), float(p[4]))
    return d

D = read(a.main, skip=("enumlibjl",))
for src in (read(f) for f in a.more):
    for code, pts in src.items():
        D.setdefault(code, {}).update(pts)

# fixed categorical order; each code also has its own marker, so identity never rests on colour alone
STY = {"enumlib": ("#eb6834", "s", "enumlib", 0.001),
       "enumlibjl": ("#1baf7a", "^", "Enumlib.jl", 0.30),
       "dsenum": ("#eda100", "D", "dsenum", 0.16),
       "necklace": ("#2a78d6", "o", "This work", 0.085)}
plt.rcParams.update({"font.size": 8, "font.family": "serif", "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "legend.frameon": False})
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.9))
for code in ("enumlib", "enumlibjl", "dsenum", "necklace"):
    c, mk, name, start = STY[code]
    ns = sorted(D[code]); t = [D[code][n][0] for n in ns]; m = [D[code][n][1] for n in ns]
    kw = dict(color=c, marker=mk, ms=4.2, lw=1.3, mec="white", mew=0.5, zorder=3 if code == "necklace" else 2)
    ax1.plot(ns, t, label=name, **kw)
    ax2.plot(ns, m, label=f"{name} (start-up {start:g} GB)", **kw)
    ax2.axhline(start, color=c, lw=0.6, ls=(0, (3, 3)), zorder=1)
ax2.axhline(24, color="#888", lw=0.7, ls="-", zorder=1)
ax2.text(10.2, 31, "machine memory, 24 GB", fontsize=6.5, color="#555")
for ax, lab, yl in ((ax1, "(a)", "Listing time (s)"), (ax2, "(b)", "Peak memory (GB)")):
    ax.set_yscale("log")
    ax.set_xlabel("Sites per cell, $n$")
    ax.set_ylabel(yl)
    ax.set_xticks(range(10, max(max(v) for v in D.values()) + 1, 2))
    ax.grid(True, which="major", color="#e8e8e8", lw=0.5, zorder=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.text(-0.16, 1.0, lab, transform=ax.transAxes, va="top", fontsize=9)
ax1.legend(loc="upper left", fontsize=7, handlelength=2.0)
ax2.set_ylim(5e-4, 200)
ax2.legend(loc="upper left", bbox_to_anchor=(0.02, 0.86), fontsize=6.3, handlelength=2.0)
fig.tight_layout(w_pad=2.2)
fig.savefig(a.o)
fig.savefig(a.o.rsplit(".", 1)[0] + ".png", dpi=220)
for code in D:
    print(code, sorted(D[code]))
