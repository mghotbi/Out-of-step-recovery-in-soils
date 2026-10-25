"""Figure 6 | A test of predictive sufficiency (schematic; no data). Minimal-text version; details are in the legend.
Run from python/:  python fig6_test.py -> ../figures/Figure_6_test.{pdf,png,svg}"""
import os
from matplotlib.patches import FancyBboxPatch, Rectangle
from style import plt, FS, FS_LETTER, save
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
FW, FH = 7.2, 3.55
fig = plt.figure(figsize=(FW, FH))
INK, MUTE = "#1E252E", "#444D58"
def L(x, y, s): fig.text(x / FW, (FH - y) / FH, s, fontsize=FS_LETTER, fontweight="bold", va="top", color=INK)

# ---------------- a ----------------
L(0.08, 0.05, "a")
steps = [("1", "Randomize", "histories", ("#DDE2E8", "#4A5563")),
         ("2", "Recover", "vary interval Δ", ("#E3D6F0", "#6A4C9C")),
         ("3", "Measure", "state z", ("#F8E0B4", "#A8661A")),
         ("4", "Challenge", "common event", ("#F4C9C3", "#A3312A")),
         ("5", "Predict", "held-out units", ("#CFE0F1", "#2F5F8F"))]
ax = fig.add_axes([0.3 / FW, (FH - 0.82) / FH, (FW - 0.42) / FW, 0.72 / FH]); ax.set_axis_off(); ax.set_xlim(0, 100); ax.set_ylim(0, 10)
w, g = 18.0, 2.5
for i, (n, ttl, sub, (fc, ec)) in enumerate(steps):
    x = i * (w + g)
    ax.add_patch(FancyBboxPatch((x, 0.3), w, 9.2, boxstyle="round,pad=0,rounding_size=1.0", fc=fc, ec=ec, lw=0.9))
    ax.text(x + w / 2, 6.4, f"{n}  {ttl}", ha="center", va="center", fontweight="bold", color=ec)
    ax.text(x + w / 2, 3.2, sub, ha="center", va="center", color=INK)
    if i < 4: ax.plot(x + w + g / 2, 4.9, marker=">", ms=3.6, color=INK, mec="none")
ax2 = fig.add_axes([0.3 / FW, (FH - 1.12) / FH, (FW - 0.42) / FW, 0.22 / FH]); ax2.set_axis_off(); ax2.set_xlim(0, 100); ax2.set_ylim(0, 1)
for x, m, dsc in ((0, "M0", "bulk state + forcing"), (35, "M1", "M0 + measured state"), (70, "M2", "M1 + history")):
    ax2.text(x, 0.5, m, fontweight="bold", va="center", color=INK); ax2.text(x + 3.6, 0.5, dsc, va="center", color=INK)

# ---------------- b ----------------
L(0.08, 1.35, "b")
ax = fig.add_axes([1.05 / FW, (FH - 3.05) / FH, 2.45 / FW, 1.55 / FH])
ax.axvspan(-0.05, 0.05, color="#9CCBA6", alpha=0.55, lw=0)
ax.axvline(-0.1, color=INK, ls=(0, (3, 2)), lw=0.8); ax.axvline(0, color="#4A5563", lw=0.7)
for yy, (lo, mid, hi), c, lab in ((2, (-0.30, -0.22, -0.14), "#B22A1F", "M1 − M0"), (1, (-0.03, 0.005, 0.04), "#5B2E9E", "M2 − M1")):
    ax.errorbar(mid, yy, xerr=[[mid - lo], [hi - mid]], fmt="o", ms=4.2, color=c, ecolor=c, elinewidth=1.8, capsize=2.2, capthick=1.0)
ax.text(-0.103, 2.72, "−δ", ha="right", va="center", color=INK); ax.text(0.0, 2.72, "±ε", ha="center", va="center", color="#24642F", fontweight="bold", bbox=dict(boxstyle="square,pad=0.08", fc="#C9E3CF", ec="none"))
ax.set_xlim(-0.42, 0.2); ax.set_ylim(0.45, 2.95); ax.set_yticks([1, 2]); ax.set_yticklabels(["M2 − M1\n(+ history)", "M1 − M0\n(+ state)"])
ax.tick_params(axis="y", length=0, pad=4); ax.spines["left"].set_visible(False)
ax.set_xticks([-0.4, -0.2, 0, 0.2]); ax.set_xticklabels(["−0.4", "−0.2", "0", "0.2"])
ax.set_xlabel("Change in held-out error")

# ---------------- c ----------------
L(3.95, 1.35, "c")
ax = fig.add_axes([3.95 / FW, (FH - 3.35) / FH, 3.15 / FW, 1.95 / FH]); ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
x0, yb, cw, ch, g = 0.30, 0.05, 0.34, 0.33, 0.02
cells = [(0, 1, "Measured state\nadequate for the\ntested prediction", ("#BFE0C6", "#2E7D3E")), (1, 1, "State or model\nincomplete", ("#DCCDEF", "#6A4C9C")),
         (0, 0, "Bulk model\nadequate", ("#D9DEE4", "#4A5563")), (1, 0, "Carrier\nunresolved", ("#F4C9C3", "#A3312A"))]
for cx, cy, txt, (fc, ec) in cells:
    X = x0 + cx * (cw + g); Y = yb + cy * (ch + g)
    ax.add_patch(Rectangle((X, Y), cw, ch, fc=fc, ec=ec, lw=0.9)); ax.text(X + cw / 2, Y + ch / 2, txt, ha="center", va="center", linespacing=1.15, color=INK)
ax.text(x0 + cw + g / 2, 0.93, "History improves prediction?", ha="center", va="center", fontweight="bold", color=INK)
ax.text(x0 + cw / 2, 0.81, "No", ha="center", va="center", color=INK); ax.text(x0 + 1.5 * cw + g, 0.81, "Yes", ha="center", va="center", color=INK)
ax.text(0.0, yb + ch + g / 2, "State\nuseful?", ha="left", va="center", fontweight="bold", color=INK, linespacing=1.1)
ax.text(x0 - 0.03, yb + 1.5 * ch + g, "Yes", ha="right", va="center", color=INK)
ax.text(x0 - 0.03, yb + ch / 2, "No", ha="right", va="center", color=INK)
save(fig, OUT, "Figure_6_test")
