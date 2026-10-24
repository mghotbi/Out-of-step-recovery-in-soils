"""Figure 1 | Out-of-step recovery and its cumulative consequences (schematic).
Run from python/:  python fig1_concept.py -> ../figures/Figure_1_out_of_step.{pdf,png,svg}
In a, function is the product of the three normalized components (an assumption for illustration);
in b, the printed debts are integrals of the plotted curves."""
import os, numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from style import plt, COL, FS, FS_LETTER, save
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")

FW, FH = 7.2, 4.05
fig = plt.figure(figsize=(FW, FH))
def ax_at(x, y_top, w, h): return fig.add_axes([x / FW, (FH - y_top - h) / FH, w / FW, h / FH])
def letter(ax, L, dx=0.45):
    bb = ax.get_position(); fig.text(bb.x0 - dx / FW, bb.y1 + 0.06 / FH, L, fontsize=FS_LETTER, fontweight="bold", va="baseline")

# ---------------- a: components recover at different times ----------------
t = np.linspace(0, 12, 2401)
bulk = 1 - np.exp(-t / 0.3); conn = 1 - np.exp(-t / 1.2); cap = 1 - 0.6 * np.exp(-t / 3)
kin = np.where(t < 1.5, 0, 1 - np.exp(-np.maximum(t - 1.5, 0) / 0.8)); F = conn * cap * kin
tb = t[np.argmax(bulk >= 0.95)]; tF = t[np.argmax(F >= 0.95)]
ax = ax_at(0.55, 0.25, 2.85, 1.75)
ax.axvspan(tb, tF, color="#F3EBD9", lw=0, zorder=0)
for x_ in (tb, tF): ax.axvline(x_, color="0.45", ls=":", lw=0.7, zorder=1)
ax.text(tb, 1.24, "$t_{bulk}$", ha="center", va="bottom"); ax.text(tF, 1.24, "$t_{F}$", ha="center", va="bottom")
for y, c, lw, lab in ((bulk, COL["Bulk"], 1.1, "bulk (water, E$_h$)"), (conn, COL["Connectivity"], 1.2, "connectivity, α"),
                      (cap, COL["Capacity"], 1.2, "capacity, Q"), (kin, COL["Kinetics"], 1.2, "kinetics, k"),
                      (F, COL["ink"], 2.0, "function")):
    ax.plot(t, y, color=c, lw=lw, label=lab, zorder=3)
ax.set_xlim(0, 12); ax.set_ylim(0, 1.22); ax.set_yticks([0, 0.5, 1]); ax.set_xticks([])
ax.set_xlabel("Time after the event (arbitrary units)"); ax.set_ylabel("Fraction recovered")
ax.text((tb + tF) / 2, 0.05, "out-of-step\nwindow", ha="center", va="bottom", fontweight="bold")
ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.02), handlelength=1.4, labelspacing=0.22, borderaxespad=0.1)
letter(ax, "a")

# ---------------- b: same endpoint, different cumulative outcome ----------------
tc = np.linspace(0, 10, 2001); ref = np.ones_like(tc)
def traj(lag, tau, over):
    y = np.ones_like(tc); d = (tc >= 1) & (tc < 2.2); y[d] = 1 - 0.75 * (tc[d] - 1) / 1.2
    a = tc >= 2.2; s = tc[a] - 2.2; sl = np.maximum(s - lag, 0)
    y[a] = np.where(s < lag, 0.25, 1 - 0.75 * np.exp(-sl / tau) + over * sl * np.exp(-sl / 0.9)); return y
fast = traj(0, 0.45, 0); slow = traj(1.6, 0.9, 1.1); dt = tc[1] - tc[0]
D1 = np.sum(ref - fast) * dt; D2 = np.sum(ref - slow) * dt
ax = ax_at(4.2, 0.25, 2.85, 1.75)
ax.axvspan(1, 2.2, color="0.9", lw=0); ax.text(1.6, 1.4, "event", ha="center", va="center", color="0.35")
ax.fill_between(tc, np.minimum(fast, 1), 1, color="0.5", alpha=0.25, lw=0)
ax.fill_between(tc, np.minimum(slow, 1), 1, color=COL["Debt"], alpha=0.3, lw=0)
ax.fill_between(tc, 1, np.maximum(slow, 1), facecolor="none", edgecolor=COL["Debt"], hatch="/////", lw=0)
ax.plot(tc, ref, color="0.5", ls="--", lw=0.7); ax.plot(tc, fast, color="0.25", lw=1.4); ax.plot(tc, slow, color=COL["Debt"], lw=1.4)
ax.annotate("fast path", xy=(2.62, 0.78), xytext=(3.3, 1.34), color="0.25", ha="left", va="center", arrowprops=dict(arrowstyle="-", lw=0.5, color="0.25"))
ax.text(4.55, 0.42, "slow path", color=COL["Debt"], ha="left", va="center")
ax.annotate("lag", xy=(3.1, 0.25), xytext=(3.1, 0.1), ha="center", va="center", arrowprops=dict(arrowstyle="-", lw=0.5))
ax.annotate("overshoot", xy=(5.3, 1.12), xytext=(6.6, 1.36), va="center", arrowprops=dict(arrowstyle="-", lw=0.5))
ax.text(9.9, 0.66, "$D_W = s\\int_W (F - F_{ref})\\,dt$", ha="right")
ax.text(9.9, 0.48, f"fast: $|D_W|$ = {D1:.2f}", ha="right", color="0.25")
ax.text(9.9, 0.33, f"slow: $|D_W|$ = {D2:.2f}", ha="right", color=COL["Debt"])
ax.set_xlim(0, 10); ax.set_ylim(0, 1.5); ax.set_xticks([]); ax.set_yticks([1], ["$F_{ref}$"])
ax.set_xlabel("Time within window W (arbitrary units)"); ax.set_ylabel("Function F")
letter(ax, "b")

# ---------------- c: the present state and its update between events ----------------
ax = ax_at(0.1, 2.62, 7.0, 1.35); ax.set_axis_off(); ax.set_xlim(0, 100); ax.set_ylim(0, 34)
bb = ax.get_position()
fig.text(bb.x0 + 0.03 / FW, bb.y1 + 0.02 / FH, "c", fontsize=FS_LETTER, fontweight="bold", va="baseline")
def box(x, y0, w, h, col, fc=None):
    ax.add_patch(FancyBboxPatch((x, y0), w, h, boxstyle="round,pad=0,rounding_size=1.2", fc=fc or col, ec="none", alpha=0.13 if fc is None else 1))
    ax.add_patch(FancyBboxPatch((x, y0), w, h, boxstyle="round,pad=0,rounding_size=1.2", fc="none", ec=col, lw=1.0))
def tile(x, w, col, name, sym):
    box(x, 19.5, w, 9.5, col)
    ax.text(x + w / 2, 25.6, name, ha="center", va="center", fontweight="bold", color=col)
    ax.text(x + w / 2, 21.9, sym, ha="center", va="center", style="italic", color=col)
ax.add_patch(Rectangle((1.3, 18.2), 66.4, 12.1, fc="none", ec="0.6", ls=(0, (2, 2)), lw=0.6))
ax.text(34.5, 31.2, "present state $z_n$", ha="center", va="bottom", color="0.3")
tile(2.5, 20.5, COL["Capacity"], "CAPACITY", "accepting inventory, Q")
tile(24.5, 20.5, COL["Connectivity"], "CONNECTIVITY", "connected fraction, α")
tile(46.5, 20.5, COL["Kinetics"], "KINETICS", "exchange coefficient, k")
# equation (1) fed by the three requirements
box(14, 1.5, 41, 7.5, "#5b7a99", fc="#EEF3F7")
ax.text(34.5, 5.25, "$C_{acc}(\\tau) = Σ_i\\,Q_i\\,\\alpha_i\\,(1 - e^{-k_i \\tau})$   (equation 1)", ha="center", va="center")
for x in (12.75, 34.75, 56.75):
    ax.add_patch(FancyArrowPatch((x, 18.0), (34.5 + (x - 34.75) * 0.55, 9.3), arrowstyle="-|>", mutation_scale=6, color="0.35", lw=0.7))
# memory and the update loop
box(73, 4.0, 25, 26.3, COL["Memory"])
ax.text(85.5, 26.6, "MEMORY", ha="center", va="center", fontweight="bold", color=COL["Memory"])
ax.text(85.5, 22.9, "update rule, M", ha="center", va="center", style="italic", color=COL["Memory"])
ax.text(85.5, 16.2, "$z_{n+1} = M(z_n, H_n)$", ha="center", va="center")
ax.text(85.5, 9.0, "persistence ratio $R = T_m/\\Delta$", ha="center", va="center", color=COL["Memory"])
ax.add_patch(FancyArrowPatch((67.9, 26.0), (72.8, 26.0), arrowstyle="-|>", mutation_scale=7, color=COL["Memory"], lw=0.9))
ax.text(70.4, 27.4, "$H_n$", ha="center", va="bottom", color=COL["Memory"])
ax.add_patch(FancyArrowPatch((72.8, 9.0), (61.5, 18.0), connectionstyle="arc3,rad=-0.5", arrowstyle="-|>", mutation_scale=7, color=COL["Memory"], lw=0.9))
ax.text(64.5, 3.2, "update to $z_{n+1}$", ha="center", va="center", color=COL["Memory"])
save(fig, OUT, "Figure_1_out_of_step")
print(f"D fast {D1:.2f}, slow {D2:.2f}; window {tb:.2f}-{tF:.2f}")
