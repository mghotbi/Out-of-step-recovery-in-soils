"""Supplementary Figure 3 | Assay capacity tracks reactive Fe(III), not total Fe inventory.
Reanalysis of Yu et al. 2026 (Commun. Earth Environ. 7, 290) source data, workbook MOESM4, sheet 'Figure 2'
(13 floodplain sediments). Derived table: ../data/yu_2026_fe_fractions_eac.csv.
Run from python/:  python figS3_capacity_inventory.py -> ../figures/Supplementary_Figure_3_capacity_inventory.{pdf,png,svg}"""
import os, numpy as np, pandas as pd
from scipy import stats
from style import plt, COL, FS_LETTER, save
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "..", "data"); OUT = os.path.join(HERE, "..", "figures")
d = pd.read_csv(os.path.join(DATA, "yu_2026_fe_fractions_eac.csv"))
FW, FH = 183 / 25.4, 2.25
fig, axs = plt.subplots(1, 3, figsize=(FW, FH)); fig.subplots_adjust(left=0.075, right=0.985, bottom=0.21, top=0.9, wspace=0.42)
panels = [("Fe_oxyhydroxide_fractions", "Fe in adsorbed and (oxyhydr)oxide\nfractions (mmol g$^{-1}$)", "a"),
          ("Fe_AAO_phyllosilicate_fractions", "Fe in later, phyllosilicate-dominated\nfractions (mmol g$^{-1}$)", "b"),
          ("HCl_FeIII", "6 M HCl-extractable Fe(III) (mmol g$^{-1}$)", "c")]
for ax, (col, lab, L) in zip(axs, panels):
    x, y = d[col].values, d.EAC.values
    if col == "HCl_FeIII":
        ax.errorbar(x, y, d.EAC_sd, d.HCl_FeIII_sd, fmt="none", ecolor="0.6", elinewidth=0.5, capsize=0)
    ax.scatter(x, y, s=16, color=COL["Capacity"], edgecolor="white", lw=0.4, zorder=3)
    lr = stats.linregress(x, y); xx = np.linspace(x.min(), x.max(), 50)
    ax.plot(xx, lr.intercept + lr.slope * xx, color="0.25", lw=0.9, ls="-" if lr.pvalue < 0.05 else "--")
    ax.text(0.04, 0.96, f"R$^2$ = {lr.rvalue**2:.2f}\nP = {lr.pvalue:.3f}", transform=ax.transAxes, va="top")
    if col == "HCl_FeIII":
        lim = [0.1, 0.46]; ax.plot(lim, lim, color="0.6", lw=0.7, ls=":"); ax.text(0.17, 0.145, "1:1", color="0.5", ha="left", va="top")
        ax.set_xlim(0.1, 0.46)
    ax.set_ylim(0.1, 0.42); ax.set_xlabel(lab); ax.set_ylabel("EAC (mmol e$^{-}$ g$^{-1}$)")
    bb = ax.get_position(); fig.text(bb.x0 - 0.055, bb.y1 + 0.02, L, fontsize=FS_LETTER, fontweight="bold", va="bottom")
save(fig, OUT, "Supplementary_Figure_3_capacity_inventory")
print("FigS3 capacity written")
