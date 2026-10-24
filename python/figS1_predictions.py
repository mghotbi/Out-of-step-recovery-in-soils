"""Supplementary Figure 1 | Restricted-model consequences P1-P4 (P5 is evaluated by the experiment in main Fig. 6). All curves analytic (theory.py).
Run from python/:  python figS1_predictions.py -> ../figures/Supplementary_Figure_1_predictions.{pdf,png,svg}"""
import os, numpy as np
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from style import plt, COL, FS, FS_LETTER, save
import theory as th
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
FW, FH = 7.2, 9.5 - 1.88
fig = plt.figure(figsize=(FW, FH))
def ax_at(x, y_top, w, h): return fig.add_axes([x / FW, (FH - y_top - h) / FH, w / FW, h / FH])
def band(y_top, code, title, col):
    a = fig.add_axes([0.08 / FW, (FH - y_top - 0.17) / FH, (FW - 0.14) / FW, 0.17 / FH]); a.set_axis_off()
    a.add_patch(Rectangle((0, 0), 0.04, 1, color=col, transform=a.transAxes))
    a.text(0.02, 0.5, code, color="white", fontweight="bold", ha="center", va="center")
    a.text(0.05, 0.5, title, fontweight="bold", va="center")
def letter(ax, L):
    bb = ax.get_position(); fig.text(bb.x0 - 0.45 / FW, bb.y1 + 0.02 / FH, L, fontsize=FS_LETTER, fontweight="bold", va="bottom")
ROWH = 1.88; AXH = 1.10; XL = [0.62, 4.15]; AW = 2.85
def row_top(i): return 0.05 + i * ROWH
logfmt = lambda ax, axis="x": None

# ---------- P1 ----------
y0 = row_top(0); band(y0, "P1", "Access changes amplitudes; kinetics moves rate constants", COL["Connectivity"])
Q = np.array([50, 50.]); kref = np.array([1, 0.01])
cases = [("reference", [1, 1], kref, "0.15", 1.2, "-"), ("all access × 0.5", [0.5, 0.5], kref, (0.18, 0.56, 0.58, 0.45), 3.5, "-"),
         ("slow access × 0.1", [1, 0.1], kref, COL["Connectivity"], 1.2, "--"), ("kinetics ÷ 3", [1, 1], kref / 3, COL["Kinetics"], 1.2, "-")]
t = np.exp(np.linspace(np.log(0.01), np.log(2000), 600))
ax = ax_at(XL[0], y0 + 0.4, AW, AXH)
for nm, a_, k_, c, lw, ls in [cases[1], cases[0], cases[2], cases[3]]:
    a_ = np.array(a_); ax.plot(t, th.c_acc(t, Q, a_, k_) / np.sum(Q * a_), color=c, lw=lw, ls=ls, label=nm)
ax.set_xscale("log"); ax.set_ylim(0, 1.6); ax.set_yticks([0, 0.5, 1]); ax.set_xlabel("Event duration τ (h)"); ax.set_ylabel("Fraction of plateau")
h_, l_ = ax.get_legend_handles_labels(); order = [1, 0, 2, 3]
ax.legend([h_[i] for i in order], [l_[i] for i in order], loc="upper left", ncol=2, columnspacing=0.8)
letter(ax, "a")
ax = ax_at(XL[1], y0 + 0.4, AW, AXH)
for j, off in ((0, 1), (2, 10 ** 0.06), (3, 1)):
    nm, a_, k_, c, lw, ls = cases[j]; amp = Q * np.array(a_); xs = k_ * off
    ax.vlines(xs, 0, amp, color=c if j else "0.15", lw=1.8, ls=ls); ax.plot(xs, amp, "o", color=c if j else "0.15", ms=3.5)
ax.set_xscale("log"); ax.set_xlim(1e-3, 3); ax.set_ylim(0, 62)
ax.text(0.018, 57, "access", color=COL["Connectivity"], va="top"); ax.text(0.018, 33, "kinetics", color=COL["Kinetics"], va="top")
ax.set_xlabel("Rate constant k$_i$ (h$^{-1}$)"); ax.set_ylabel("Amplitude Q$_i$α$_i$")
letter(ax, "b")

# ---------- P2 ----------
y0 = row_top(1); band(y0, "P2", "Event duration can reverse exchange rankings", COL["Capacity"])
tb = np.linspace(0, 120, 1201); tx = th.crossover_time(20, 1, 60, 0.02)
ax = ax_at(XL[0], y0 + 0.4, AW, AXH)
ax.plot(tb, th.c_acc(tb, 20, 1, 1), color=COL["Kinetics"], label="small, fast"); ax.plot(tb, th.c_acc(tb, 60, 1, 0.02), color=COL["Capacity"], label="large, slow")
ax.axvline(tx, color="0.3", ls="--", lw=0.6); ax.text(tx + 3, 4, f"τ* = {tx:.1f} h")
ax.set_xlim(0, 120); ax.set_ylim(0, 80); ax.set_xlabel("Event duration τ (h)"); ax.set_ylabel("Exchanged (mmol e$^{-}$ kg$^{-1}$)")
ax.legend(loc="upper center", ncol=2); letter(ax, "c")
ax = ax_at(XL[1], y0 + 0.4, AW - 0.45, AXH)
xr = np.exp(np.linspace(np.log(0.05), np.log(20), 121)); yr = np.exp(np.linspace(np.log(0.002), np.log(20), 121))
TT = np.full((yr.size, xr.size), np.nan)
for i, yk in enumerate(yr):
    for j, xp in enumerate(xr):
        if (1 * 1 - xp * yk) * (1 - xp) < 0: TT[i, j] = th.crossover_time(1, 1, xp, yk)
cmap = LinearSegmentedColormap.from_list("cr", ["#FFF7BC", "#FEC44F", "#D95F0E", "#662506"])
pc = ax.pcolormesh(xr, yr, np.clip(TT, 0.1, 1000), cmap=cmap, norm=LogNorm(0.1, 1000), shading="auto", rasterized=True)
ax.set_xscale("log"); ax.set_yscale("log"); ax.plot(3, 0.02, "o", mfc="white", mec="k", ms=5)
ax.text(0.07, 0.012, "no crossing", color="0.3"); ax.text(2.5, 3, "no crossing", color="0.3"); ax.text(2.2, 0.0035, "crossing", color="white")
ax.set_xlabel("Plateau ratio (Qα)$_2$/(Qα)$_1$"); ax.set_ylabel("Rate ratio k$_2$/k$_1$")
cb = fig.colorbar(pc, ax=ax, fraction=0.06, pad=0.03); cb.set_label("τ* k$_1$"); cb.outline.set_linewidth(0.5)
letter(ax, "d")

# ---------- P3 ----------
y0 = row_top(2); band(y0, "P3", "A late terminal sink creates a redox debt", COL["Debt"])
ts = np.linspace(0, 30, 1201); f0 = 0.1 * th.x_traj(ts, 0); f4 = 0.1 * th.x_traj(ts, 4)
ax = ax_at(XL[0], y0 + 0.4, AW, AXH)
ax.axvspan(0, 12, color="0.93", lw=0); ax.text(6, 0.93, "N$_2$O formation", ha="center", color="0.35")
ax.fill_between(ts, f0, np.maximum(f4, f0), color=COL["Debt"], alpha=0.3, lw=0)
for L, c in ((0, "0.25"), (4, COL["Kinetics"]), (6, "#7b1f1f")):
    ax.plot(ts, 0.1 * th.x_traj(ts, L), color=c, lw=1.2, label=f"sink lag {L} h: E = {float(th.escape(L)):.2f} µmol N")
ax.text(0.6, 0.62, "debt (L = 4 h)", fontweight="bold", color=COL["Debt"])
ax.set_xlim(0, 26); ax.set_ylim(0, 1); ax.set_xlabel("Time (h)"); ax.set_ylabel("N$_2$O escape flux (µmol N h$^{-1}$)")
ax.legend(loc="upper right", bbox_to_anchor=(1, 0.86)); letter(ax, "e")
ax = ax_at(XL[1], y0 + 0.4, AW, AXH)
Ls = np.linspace(0, 30, 301); ax.plot(Ls, th.excess_escape(Ls), color=COL["Debt"], lw=1.4, label="exact, equation (2)")
Lq = Ls[Ls <= 9]; ax.plot(Lq, (1 / 1.1) * 0.1 * Lq ** 2 / 2, "--", color="0.3", lw=0.8, label="short-lag approximation")
ceil = 1 * 12 / 1.1; ax.axhline(ceil, color="0.6", ls=":", lw=0.8); ax.text(0.5, ceil + 0.35, "sink never returns", color="0.4")
ax.plot(4, float(th.excess_escape(4)), "o", color=COL["Kinetics"], ms=3.5)
ax.set_xlim(0, 30); ax.set_ylim(0, 12.5); ax.set_xlabel("Sink lag L (h)"); ax.set_ylabel("Excess escape (µmol N)")
ax.legend(loc="lower right"); letter(ax, "f")

# ---------- P4 ----------
y0 = row_top(3); band(y0, "P4", "Recurrence faster than relaxation accumulates legacy", COL["Memory"])
ax = ax_at(XL[0], y0 + 0.4, AW, AXH)
n = np.arange(1, 13)
for R, c in zip((0.25, 1, 4), ("#c9bce0", "#9a7cc4", COL["Memory"])):
    ax.plot(n, th.legacy_pre(12, 1, R), "-o", color=c, ms=2.5, label=f"$R$ = {R}")
    ax.axhline(float(th.carried_legacy(R)), color=c, ls=":", lw=0.7)
ax.set_xlim(1, 12); ax.set_ylim(0, 3.7); ax.set_xlabel("Event number n"); ax.set_ylabel("Pre-event legacy (increments)")
ax.legend(loc="upper left"); letter(ax, "g")
ax = ax_at(XL[1], y0 + 0.4, AW, AXH)
R = np.exp(np.linspace(np.log(0.1), np.log(31.6), 300)); r = np.exp(-1 / R)
ax.plot(R, th.carried_legacy(R), color=COL["Memory"], lw=1.4, label="equation (3)")
RR = R[R > 1.5]; ax.plot(RR, RR - 0.5, "--", color="0.3", lw=0.8, label="large-$R$ limit")
ax.plot(R, r / (1 - r * (1 - 1 / 5)), "-.", color="0.55", lw=0.9, label="bounded (max = 5)")
ax.axhline(1, color="0.7", ls=":", lw=0.6); ax.axvline(1 / np.log(2), color="0.7", ls=":", lw=0.6)
ax.text(1 / np.log(2) * 1.08, 0.013, "1/ln 2", color="0.4")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(0.01, 40)
ax.set_xlabel("Persistence ratio $R$ = T$_m$/Δ"); ax.set_ylabel("Legacy B (increments)")
ax.legend(loc="lower right"); letter(ax, "h")

save(fig, OUT, "Supplementary_Figure_1_predictions")
