"""Supplementary Fig. 1: illustrative example (file name kept for stability). Second-challenge N2O escape versus recovery interval.
Assumptions (illustrative): identical formation, escape and post-activation consumption in both cores
(p = 1 umol N h-1, T = 12 h, l = 0.1 h-1, k = 1 h-1, equation (2)); core A's sink-activation lag relaxes
exponentially between events, L(D) = L0 exp(-D/Tm), with Tm = 10 d and L0 set so that L = 4 h at D = 4 d."""
import numpy as np, matplotlib.pyplot as plt
from style import COL, save, letter, FS
from theory import escape
OUT = "../figures"
Tm = 10.0; L0 = 4.0 * np.exp(4.0 / Tm)   # lag of exactly 4 h at the short interval (4 d)
D = np.linspace(0, 40, 400)
EA = escape(0.0)
EB = np.array([escape(L0 * np.exp(-d / Tm)) for d in D])
short, long_ = 4.0, 30.0
eps = 0.05 * EA   # illustrative equivalence precision (+/-5%)

fig = plt.figure(figsize=(89 / 25.4, 58 / 25.4))
ax = fig.add_axes([0.16, 0.19, 0.80, 0.71])
ax.axhspan(EA - eps, EA + eps, color="#D9DEE4", lw=0)
ax.plot(D, np.full_like(D, EA), color=COL["Connectivity"], lw=1.4)
ax.plot(D, EB, color=COL["Kinetics"], lw=1.4)
for d in (short, long_):
    ax.axvline(d, color="#4A5563", lw=0.6, ls=(0, (3, 2)))
    eb = escape(L0 * np.exp(-d / Tm))
    ax.plot([d], [eb], "o", ms=3.6, color=COL["Kinetics"]); ax.plot([d], [EA], "o", ms=3.6, color=COL["Connectivity"])
ax.annotate("", xy=(short + 0.8, escape(L0*np.exp(-short/Tm))), xytext=(short + 0.8, EA),
            arrowprops=dict(arrowstyle="<->", lw=0.7, color="#222222", shrinkA=0, shrinkB=0))
ax.text(long_, 2.55, "long", ha="center", va="bottom"); ax.text(short, 2.55, "short", ha="center", va="bottom")
ax.text(13.0, 1.36, "Core A · delayed sink", ha="left", va="bottom", color=COL["Kinetics"])
ax.text(13.0, EA - 0.08, "Core B · active sink", ha="left", va="top", color=COL["Connectivity"])
EBs = escape(L0 * np.exp(-short / Tm))
ax.text(short + 0.9, EBs + 0.02, "Δ = 4 d: %.2f µmol N" % EBs, ha="left", va="bottom", color=COL["Kinetics"])
ax.text(short - 0.5, (EA + EBs) / 2 - 0.08, "excess\n%.2f" % (EBs - EA), ha="right", va="center", color="#222222", linespacing=1.1)
ax.set_title("Hypothetical model response", loc="right", color="#5A6270", pad=2)
ax.text(30.8, EA - eps - 0.02, "illustrative\ntolerance", ha="left", va="top", color="#5A6270", linespacing=1.1)
ax.set_xlim(0, 40); ax.set_ylim(0.8, 2.75)
ax.set_xlabel("Recovery interval before second challenge, Δ (d)")
ax.set_ylabel("Cumulative N$_2$O escape (µmol N)")
import sys
PANEL = sys.argv[1] if len(sys.argv) > 1 else ""   # pass "b" once the artwork panel a is added
if PANEL:
    letter(ax, PANEL, x=-0.10)
save(fig, OUT, "Supplementary_Figure_1_worked_example" + ("_panel_" + PANEL if PANEL else ""))
print("A", round(EA, 3), "B short", round(escape(L0*np.exp(-short/Tm)), 3), "lag", round(L0*np.exp(-short/Tm), 2),
      "B long", round(escape(L0*np.exp(-long_/Tm)), 3), "lag", round(L0*np.exp(-long_/Tm), 2))
