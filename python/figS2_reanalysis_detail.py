"""Supplementary Figure 2 | Published nitrogen trajectories behind Fig. 5g and 5i–k (panels a–c).
Run from python/:  python figS2_reanalysis_detail.py -> ../figures/Supplementary_Figure_2_reanalysis.{pdf,png,svg}"""
import os, numpy as np, pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter
from style import plt, COL, HIST, FS, FS_LETTER, save

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "..", "data"); OUT = os.path.join(HERE, "..", "figures")
gas = pd.read_csv(os.path.join(DATA, "sennett_final_incubation_gas.csv"))
M = pd.read_csv(os.path.join(DATA, "figure3_derived_metrics.csv"))
co2 = pd.read_csv(os.path.join(DATA, "liu_2025_panel_h_co2_efflux.csv"))

FW, FH = 7.2, 0.52 + 1.45 + 0.6 + 1.45 + 0.55
fig = plt.figure(figsize=(FW, FH))
def ax_at(x, y_top, w, h): return fig.add_axes([x / FW, (FH - y_top - h) / FH, w / FW, h / FH])
def strip(y_top, col, label, claim, x0=0.08, x1=FW - 0.06, cx=0.36):
    a = fig.add_axes([x0 / FW, (FH - y_top - 0.17) / FH, (x1 - x0) / FW, 0.17 / FH]); a.set_axis_off()
    a.add_patch(Rectangle((0, 0), 1, 1, color=col, transform=a.transAxes))
    a.text(0.008 * (FW - 0.14) / (x1 - x0), 0.5, label, color="white", fontweight="bold", va="center", fontsize=FS)
    a.text(cx, 0.5, claim, color="white", va="center", fontsize=FS)
def head(ax, L, title):
    bb = ax.get_position()
    fig.text(bb.x0 - 0.42 / FW, bb.y1 + 0.08 / FH, L, fontsize=FS_LETTER, fontweight="bold", va="baseline")
    fig.text(bb.x0 - 0.26 / FW, bb.y1 + 0.08 / FH, title, fontsize=FS, fontweight="bold", va="baseline")

AXH = 1.45; ws = 2.75; xs = [0.62, 4.15, 0.62]
strip(0.05, COL["Memory"], "OXYGEN-CONDITIONING EXPERIMENT", "")
yt = 0.05 + 0.17 + 0.30
yt2 = yt + AXH + 0.6
ax = ax_at(xs[0], yt, ws, AXH)
for (h, r), q in gas.groupby(["history", "replicate"]):
    q = q.sort_values("time_h"); ax.plot(q.time_h, q.N2_cum, color=HIST[h], lw=0.9)
ax.set_xlim(0, 112); ax.set_ylim(0, 11.5); ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel("Time in common challenge (h)"); ax.set_ylabel("Cumulative N$_2$ (µmol N)")
for h, yy in zip(HIST, (0.30, 0.19, 0.08)):
    ax.text(0.97, yy, h, color=HIST[h], fontweight="bold", transform=ax.transAxes, ha="right")
head(ax, "a", "")

ax = ax_at(xs[1], yt, ws, AXH)
for (h, r), q in gas.groupby(["history", "replicate"]):
    q = q.sort_values("time_h"); ax.plot(q.time_h, q.NO, color=HIST[h], lw=0.9)
ax.set_xlim(0, 112); ax.set_xticks([0, 25, 50, 75, 100]); ax.set_ylim(bottom=0)
ax.set_xlabel("Time in common challenge (h)"); ax.set_ylabel("NO in vial (µmol N)")
head(ax, "b", "")

ax = ax_at(xs[2], yt2, ws, AXH)
metrics = [("N2_final", "Terminal\nN$_2$", "s"), ("N2O_integral_1_24h", "1–24 h\n∫N$_2$O dt", "^"), ("N2O_time_integral", "Complete-record\n∫N$_2$O dt", "o")]
ax.axhline(1, color="0.6", lw=0.6)
for j, (m, lab, mk) in enumerate(metrics):
    ox = M[M.history == "Ox"][m].mean()
    for i, h in enumerate(HIST):
        vals = M[M.history == h][m].values / ox
        x = j + (i - 1) * 0.24
        ax.scatter(np.full(3, x), vals, s=9, marker=mk, facecolor="white", edgecolor=HIST[h], lw=0.6, zorder=3)
        ax.scatter([x], [vals.mean()], s=22, marker=mk, color=HIST[h], zorder=4)
ax.set_yscale("log"); ax.set_ylim(0.5, 20)
ax.yaxis.set_major_locator(FixedLocator([0.5, 1, 2, 5, 10, 20])); ax.yaxis.set_major_formatter(FixedFormatter(["0.5", "1", "2", "5", "10", "20"]))
ax.yaxis.set_minor_formatter(NullFormatter())
ax.set_xticks(range(3), [m[1] for m in metrics]); ax.set_xlim(-0.5, 2.5); ax.set_ylabel("Ratio to Ox mean")
head(ax, "c", "")

save(fig, OUT, "Supplementary_Figure_2_reanalysis")
