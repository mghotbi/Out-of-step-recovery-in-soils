"""Figure 5 | Published observations organized by the four framework components.
Every panel shows source data from ONE study; no panel combines studies. Studies: Li et al. 2025 (a-c, e-f),
Patzner et al. 2020 (d), Sennett et al. 2024 (g, i-k), Liu et al. 2025 (h).
Run from python/:  python fig5_evidence.py -> ../figures/Figure_5_evidence.{pdf,png,svg}"""
import os, numpy as np, pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter
from scipy.stats import spearmanr
from style import plt, COL, HIST, FS, FS_LETTER, save
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "..", "data"); OUT = os.path.join(HERE, "..", "figures")
rd = lambda f: pd.read_csv(os.path.join(DATA, f))
logticks = lambda axis, v: (axis.set_major_locator(FixedLocator(v)), axis.set_major_formatter(FixedFormatter([f"{x:g}" for x in v])), axis.set_minor_formatter(NullFormatter()))
CMP = ["Iron plaque", "Rhizosphere", "Bulk soil"]; CLAB = ["Plaque", "Rhizo-\nsphere", "Bulk\nsoil"]
FW = 7.2; BAND, GAP, AXH, XLAB = 0.17, 0.26, 1.25, 0.50
ROW = BAND + GAP + AXH + XLAB; FH = 4 * ROW + 0.05
fig = plt.figure(figsize=(FW, FH))
def band(y_top, comp, mnem, src):
    a = fig.add_axes([0.08 / FW, (FH - y_top - BAND) / FH, (FW - 0.14) / FW, BAND / FH]); a.set_axis_off()
    a.add_patch(Rectangle((0, 0), 1, 1, color=COL[comp], transform=a.transAxes))
    a.text(0.008, 0.5, comp.upper(), color="white", fontweight="bold", va="center")
def row(y_top, xs, w):
    return [fig.add_axes([x / FW, (FH - y_top - BAND - GAP - AXH) / FH, w / FW, AXH / FH]) for x in xs]
def letter(ax, L, dx=0.45):
    bb = ax.get_position(); fig.text(bb.x0 - dx / FW, bb.y1 + 0.07 / FH, L, fontsize=FS_LETTER, fontweight="bold", va="bottom")
X4, W4 = (0.52, 2.30, 4.08, 5.86), 1.20
X2, W2 = (0.62, 4.15), 2.75
X3, W3 = (0.55, 2.95, 5.35), 1.62

# ---------------- CAPACITY: a–c Li et al. 2025; d Patzner et al. 2020 ----------------
y = 0.05; band(y, "Capacity", "reserve", "a–c, Li et al. 2025;  d, Patzner et al. 2020")
ax1, ax2, ax3, ax4 = row(y, X4, W4)
eec = rd("li_2025_eec_primary_source.csv").set_index("compartment")
frac = rd("li_2025_fe_fractions_primary_source.csv"); pcf = frac[frac.fraction_index == 2].set_index("compartment").percent_of_total_Fe
goe = rd("li_2025_goethite_eec_primary_source.csv").set_index("goethite")
v = eec.loc[CMP, "eec"].values; sd = eec.loc[CMP, "sd"].values
ax1.bar(range(3), v, 0.62, color=COL["Capacity"], alpha=0.4, edgecolor=COL["Capacity"], lw=0.8)
ax1.errorbar(range(3), v, sd, fmt="none", ecolor="k", elinewidth=0.6, capsize=2)
for i in range(3): ax1.text(i, v[i] + sd[i] + 0.12, f"×{v[i] / v[2]:.2f}", ha="center", va="bottom")
ax1.set_xticks(range(3), CLAB); ax1.set_ylim(0, 4.4); ax1.set_yticks([0, 1, 2, 3, 4]); ax1.set_xlim(-0.6, 2.6)
ax1.set_ylabel("EAC + EDC (mmol e$^{-}$ g$^{-1}$)"); letter(ax1, "a")
ax2.bar(range(3), pcf.loc[CMP].values, 0.62, color="#7A4E1D", alpha=0.75, edgecolor="#7A4E1D", lw=0.8)
for i, c in enumerate(CMP): ax2.text(i, pcf[c] + 2, f"{pcf[c]:.0f}%", ha="center", va="bottom")
ax2.set_xticks(range(3), CLAB); ax2.set_ylim(0, 100); ax2.set_xlim(-0.6, 2.6); ax2.set_ylabel("Poorly crystalline Fe (%)"); letter(ax2, "b")
base = goe.loc["without redox cycling"].EEC; post = goe.loc["after 30 d diel redox cycling"].EEC
ax3.bar([0, 1], [1.0, post / base], 0.6, color=["0.72", COL["Capacity"]], edgecolor="none")
ax3.text(1, post / base + 0.04, f"×{post / base:.2f}", ha="center", va="bottom")
ax3.set_xticks([0, 1], ["Goethite", "After 30 d\ncycling"]); ax3.set_ylim(0, 2.0); ax3.set_yticks([0, 0.5, 1, 1.5, 2]); ax3.set_xlim(-0.6, 1.6)
ax3.set_ylabel("EEC relative to goethite"); letter(ax3, "c")
fe = rd("patzner_panel_p_mineral_fe_oc_measured.csv"); fe = fe[fe.metric == "Reactive Fe"]
hz = ["Organic horizon", "Transition zone", "Mineral horizon"]; st = ["Palsa", "Bog", "Fen"]; shades = ["#EAD7B7", "#C99A55", "#7A4E1D"]
bottom = np.zeros(3)
for j, h in enumerate(hz):
    vals = np.array([100 * fe[(fe.stage == s) & (fe.horizon == h)].value.iloc[0] / fe[fe.stage == s].value.sum() for s in st])
    ax4.bar(range(3), vals, 0.62, bottom=bottom, color=shades[j], edgecolor="white", lw=0.5, label=h.split()[0]); bottom += vals
ax4.set_xticks(range(3), st); ax4.set_ylim(0, 100); ax4.set_yticks([0, 25, 50, 75, 100]); ax4.set_xlim(-0.55, 3.75); ax4.set_ylabel("Share of reactive Fe (%)")
fen = [100 * fe[(fe.stage == "Fen") & (fe.horizon == h)].value.iloc[0] / fe[fe.stage == "Fen"].value.sum() for h in hz]
ax4.text(2.4, fen[0] + fen[1] + fen[2] / 2, "Mineral", va="center"); ax4.text(2.4, fen[0] + fen[1] / 2, "Transition", va="center")
ax4.annotate("Organic", xy=(2.31, fen[0] / 2), xytext=(2.4, 6), va="center", arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4"))
ax4.spines["bottom"].set_bounds(-0.55, 2.55)
letter(ax4, "d")

# ---------------- CONNECTIVITY: e, f Li et al. 2025 ----------------
y += ROW; band(y, "Connectivity", "reach", "Li et al. 2025")
ax1, ax2 = row(y, X2, W2)
diel = rd("li_2025_diel_feII_P_primary_source.csv")
cc = {"Iron plaque": COL["Connectivity"], "Rhizosphere": "#7FBFC3", "Bulk soil": "0.55"}
for c, g in diel.groupby("compartment", sort=False):
    m = g.FeII.mean()
    ax1.errorbar(g.time_h, g.FeII / m, g.FeII_sd / m, color=cc[c], marker="o", ms=3, lw=1.1, elinewidth=0.5, capsize=1.5)
    ax1.text(62.5, {"Iron plaque": 1.40, "Rhizosphere": 1.12, "Bulk soil": 0.86}[c], {"Iron plaque": "Plaque"}.get(c, c), color=cc[c], va="center")
ax1.axhline(1, color="0.7", lw=0.5, ls="--"); ax1.set_xlim(-3, 80); ax1.set_xticks([0, 12, 24, 36, 48, 60]); ax1.set_ylim(0.4, 1.6); ax1.set_yticks([0.5, 1, 1.5])
ax1.set_xlabel("Time (h)"); ax1.set_ylabel("Fe(II) / compartment mean"); letter(ax1, "e")
# f: diel variation against distance from the root (improved from bar chart): one line per solute, plus retained capacity (same study)
dist = [0, 1, 2]; dlab = ["Plaque\n(root surface)", "Rhizosphere\n(< 4 mm)", "Bulk soil\n(> 2 cm)"]
sw = {v_: [np.mean(np.abs(np.diff(diel[diel.compartment == c][v_].values))) for c in CMP] for v_ in ("FeII", "P")}
for v_, lab, col, mk in (("FeII", "Fe(II) variation", COL["Connectivity"], "o"), ("P", "Available-P variation", "#7FBFC3", "s")):
    r = np.array(sw[v_]) / sw[v_][0]
    ax2.plot(dist, r, marker=mk, ms=4, lw=1.3, color=col, label=lab)
    for x, yv in zip(dist[1:], r[1:]):
        if v_ == "P": ax2.text(x + 0.09, yv + 0.05, f"{yv:.2f}", color=col, va="center", ha="left")
        elif x == 1: ax2.text(x - 0.09, yv - 0.04, f"{yv:.2f}", color=col, va="top", ha="right")
        else: ax2.text(x - 0.12, 0.012, f"{yv:.2f}", color=col, va="bottom", ha="right")
capr = eec.loc[CMP, "eec"].values / eec.loc["Iron plaque", "eec"]
ax2.plot(dist, capr, marker="D", ms=3.5, lw=1.0, ls="--", color=COL["Capacity"], mfc="white", label="Assay capacity (EAC + EDC)")
ax2.text(2.09, capr[2] + 0.03, f"{capr[2]:.2f}", color=COL["Capacity"], va="bottom")
ax2.set_xticks(dist, dlab); ax2.set_xlim(-0.3, 2.45); ax2.set_ylim(0, 1.15); ax2.set_yticks([0, 0.25, 0.5, 0.75, 1])
ax2.set_ylabel("Relative to plaque"); ax2.legend(loc="upper right", handlelength=1.8, labelspacing=0.3); letter(ax2, "f")

# ---------------- KINETICS: g Sennett et al. 2024; h Liu et al. 2025 ----------------
y += ROW; band(y, "Kinetics", "rate", "g, Sennett et al. 2024;  h, Liu et al. 2025")
ax1, ax2 = row(y, X2, W2)
M = rd("figure3_derived_metrics.csv")
rho_c = spearmanr(M.N2_half_time_h, M.N2O_integral_1_24h).statistic; rho_w = spearmanr(M.N2_half_time_h, M.N2O_time_integral).statistic
for h in HIST:
    q = M[M.history == h]
    ax1.scatter(q.N2_half_time_h, q.N2O_integral_1_24h, s=20, color=HIST[h], edgecolor="white", lw=0.4, zorder=4, label=h)
ax1.set_xscale("log"); ax1.set_xlim(10, 110); ax1.set_ylim(35, 95)
logticks(ax1.xaxis, [10, 20, 50, 100]); ax1.set_yticks([40, 60, 80])
ax1.set_xlabel("Terminal-record N$_2$ half-time (h)"); ax1.set_ylabel("N$_2$O amount × time, 1–24 h (µmol N h)")
ax1.text(0.97, 0.95, f"ρ = {rho_c:.2f}", transform=ax1.transAxes, ha="right", va="top")
ax1.legend(loc="upper left", handletextpad=0.1, labelspacing=0.25); letter(ax1, "g")
co2 = rd("figure3_derived_co2_windows.csv")
soils = {"Paddy soil": "#5b3a1a", "Sandy soil": "#c2a36b"}
for s_, c in soils.items():
    q = co2[co2.soil == s_].set_index("window")
    ster = q.loc[["0.5-6 h", "0.5-48 h"], "sterilized_pct"].values; tba = q.loc[["0.5-6 h", "0.5-48 h"], "tba_pct"].values
    ax2.plot([0, 1], ster, "-o", color=c, ms=3.5, lw=1.3); ax2.plot([0, 1], tba, "--s", color=c, ms=3.2, lw=1.0, mfc="white")
    for yv in ster: pass
ax2.set_xlim(-0.35, 1.35); ax2.set_ylim(0, 125); ax2.set_yticks([0, 25, 50, 75, 100]); ax2.set_xticks([0, 1], ["0.5–6 h", "0.5–48 h"])
ax2.set_xlabel("Integration window after rewetting"); ax2.set_ylabel("Treated / untreated CO$_2$ (%)")
hd = [Line2D([], [], color="0.3", marker="o", ms=3.5, lw=1.3), Line2D([], [], color="0.3", ls="--", marker="s", ms=3.2, mfc="white", lw=1.0),
      Line2D([], [], color=soils["Paddy soil"], lw=2.5), Line2D([], [], color=soils["Sandy soil"], lw=2.5)]
ax2.legend(hd, ["γ-sterilized", "TBA-treated", "Paddy soil", "Sandy soil"], loc="upper center", ncol=4, columnspacing=0.8, handlelength=1.5, borderaxespad=0.1)
letter(ax2, "h")

# ---------------- MEMORY: i–k Sennett et al. 2024 ----------------
y += ROW; band(y, "Memory", "recall", "Sennett et al. 2024")
ax1, ax2, ax3 = row(y, X3, W3)
gas = rd("sennett_final_incubation_gas.csv"); plfa = rd("sennett_plfa_total.csv")
for (h, r), q in gas.groupby(["history", "replicate"]):
    q = q.sort_values("time_h"); ax1.plot(q.time_h, q.N2O, color=HIST[h], lw=0.9)
ax1.set_xlim(0, 112); ax1.set_ylim(0, 14); ax1.set_xticks([0, 25, 50, 75, 100]); ax1.set_yticks([0, 4, 8, 12])
ax1.set_xlabel("Time in common challenge (h)"); ax1.set_ylabel("N$_2$O (µmol N per vial)")
for h, yy in zip(HIST, (0.90, 0.79, 0.68)): ax1.text(0.88, yy, h, color=HIST[h], transform=ax1.transAxes, fontweight="bold")
letter(ax1, "i")
pl = plfa.assign(v=plfa.total_plfa_pmol_g / 1000)
metrics = [("PLFA", None), ("Terminal\nN$_2$", "N2_final"), ("N$_2$O, 1–24 h", "N2O_integral_1_24h")]
ax2.axvspan(1.5, 2.5, color="0.94", lw=0); ax2.axhline(1, color="0.6", lw=0.6, ls="--")
for j, (lab, col) in enumerate(metrics):
    for k, h in enumerate(("SA", "LA")):
        vv = (pl[pl.treatment == h].v.values / pl[pl.treatment == "Ox"].v.mean()) if col is None else (M[M.history == h][col].values / M[M.history == "Ox"][col].mean())
        x = j + (k - 0.5) * 0.32
        ax2.scatter(np.full(vv.size, x), vv, s=8, facecolor="white", edgecolor=HIST[h], lw=0.6, zorder=3)
        ax2.scatter([x], [vv.mean()], s=20, color=HIST[h], zorder=4, label=f"{h}/Ox" if j == 0 else None)
ax2.set_ylim(0.4, 1.75); ax2.set_yticks([0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6]); ax2.set_xticks(range(3), [m[0] for m in metrics]); ax2.set_xlim(-0.5, 2.5)
ax2.text(0.5, 1.66, "biomass and\nterminal yield", ha="center", va="center", linespacing=1.0); ax2.text(2.0, 1.66, "accumulation", ha="center", va="center")
ax2.set_ylabel("Ratio to Ox mean"); ax2.legend(loc="lower left", handletextpad=0.1, labelspacing=0.25); letter(ax2, "j")
cw = rd("running_integrals_common_window_v16.csv")
cw["history"] = cw.series.str.split("_").str[0]
piv = cw.pivot_table(index="end_h", columns="series", values="amount_time_integral")
mean_h = {h: piv[[c for c in piv.columns if c.startswith(h + "_")]].mean(axis=1) for h in HIST}
ax3.axhline(1, color="0.6", lw=0.6, ls="--")
for c in piv.columns:
    h = c.split("_")[0]
    if h != "Ox": ax3.plot(piv.index, piv[c] / mean_h["Ox"], color=HIST[h], lw=0.5, alpha=0.45)
for h in ("SA", "LA"): ax3.plot(piv.index, mean_h[h] / mean_h["Ox"], color=HIST[h], lw=1.6, label=f"{h}/Ox")
for W in (12, 24):
    r = float(np.interp(W, piv.index, mean_h["SA"] / mean_h["Ox"])); ax3.plot([W], [r], "o", ms=3, color="k", zorder=4)
    ax3.text(W + (0 if W == 12 else -0.8), r * (0.72 if W == 12 else 1.05), f"{r:.2f}", ha="center" if W == 12 else "right", va="bottom" if W == 24 else "center")
ax3.set_yscale("log"); ax3.set_ylim(0.2, 3); logticks(ax3.yaxis, [0.2, 0.5, 1, 2]); ax3.set_xlim(0, 25); ax3.set_xticks([0, 6, 12, 18, 24])
ax3.set_xlabel("Window end W (h), from 1 h"); ax3.set_ylabel("N$_2$O amount × time / Ox")
ax3.legend(loc="upper left", handlelength=1.2, labelspacing=0.25); letter(ax3, "k")
save(fig, OUT, "Figure_5_evidence")
print("Fig5 written; SA/Ox 12 h", float(np.interp(12, piv.index, mean_h["SA"] / mean_h["Ox"])), "24 h", float(np.interp(24, piv.index, mean_h["SA"] / mean_h["Ox"])))
