"""Recompute every number quoted in the manuscript text and figure legends from the source data and the analytical models.
Run from python/:  python verify_reported_values.py  -> ../data/verification_reported_values.csv (must print ALL PASS)"""
import os, re, json, numpy as np, pandas as pd
from scipy.stats import spearmanr
import theory as th
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, ".."); D = os.path.join(ROOT, "data")
rd = lambda f: pd.read_csv(os.path.join(D, f))
rows = []
def check(name, value, quoted, tol): rows.append(dict(check=name, computed=value, quoted=quoted, tolerance=tol, pass_=abs(value - quoted) <= tol))
def flag(name, ok, detail=""): rows.append(dict(check=name, computed=detail, quoted="", tolerance="", pass_=bool(ok)))
eec = rd("li_2025_eec_primary_source.csv").set_index("compartment").eec
check("EEC rhizo/bulk", eec["Rhizosphere"] / eec["Bulk soil"], 1.25, .005); check("EEC plaque/bulk", eec["Iron plaque"] / eec["Bulk soil"], 1.92, .005)
check("EEC bulk/plaque", eec["Bulk soil"] / eec["Iron plaque"], .52, .005)
fr = rd("li_2025_fe_fractions_primary_source.csv"); pc = fr[fr.fraction_index == 2].set_index("compartment").percent_of_total_Fe
check("poorly crystalline plaque %", pc["Iron plaque"], 75, .5); check("poorly crystalline rhizo %", pc["Rhizosphere"], 30, 3); check("poorly crystalline bulk %", pc["Bulk soil"], 30, 3)
g = rd("li_2025_goethite_eec_primary_source.csv").set_index("goethite").EEC
check("goethite ratio", g["after 30 d diel redox cycling"] / g["without redox cycling"], 1.62, .005)
fe = rd("patzner_panel_p_mineral_fe_oc_measured.csv"); fe = fe[fe.metric == "Reactive Fe"]
sh = lambda st, hz: 100 * fe[(fe.stage == st) & (fe.horizon == hz)].value.iloc[0] / fe[fe.stage == st].value.sum()
check("Patzner palsa mineral %", sh("Palsa", "Mineral horizon"), 91, .5); check("Patzner bog transition %", sh("Bog", "Transition zone"), 66, .5)
diel = rd("li_2025_diel_feII_P_primary_source.csv")
sw = lambda v, c: np.mean(np.abs(np.diff(diel[diel.compartment == c][v].values)))
for v, c, q in (("FeII", "Rhizosphere", .40), ("FeII", "Bulk soil", .08), ("P", "Rhizosphere", .52), ("P", "Bulk soil", .14)):
    check(f"{v} swing {c}/plaque", sw(v, c) / sw(v, "Iron plaque"), q, .005)
M = rd("figure3_derived_metrics.csv")
check("half-time min", M.N2_half_time_h.min(), 14, .5); check("half-time max", M.N2_half_time_h.max(), 78, .5)
check("Spearman rho complete record", spearmanr(M.N2_half_time_h, M.N2O_time_integral).statistic, .98, .005)
r124 = spearmanr(M.N2_half_time_h, M.N2O_integral_1_24h); check("Spearman rho 1-24 h", r124.statistic, .37, .005); check("P 1-24 h", r124.pvalue, .33, .005)
check("1-24 h integral min", M.N2O_integral_1_24h.min(), 45, .5); check("1-24 h integral max", M.N2O_integral_1_24h.max(), 79, .5)
mean = M.groupby("history").mean(numeric_only=True)
for col, h, q, t in (("N2_final", "SA", .91, .005), ("N2_final", "LA", .86, .005), ("N2O_integral_1_24h", "SA", 1.19, .005), ("N2O_integral_1_24h", "LA", 1.24, .005)):
    check(f"{col} {h}/Ox", mean[col][h] / mean[col]["Ox"], q, t)
pl = rd("sennett_plfa_total.csv").groupby("treatment").total_plfa_pmol_g.mean()
check("PLFA SA/Ox", pl["SA"] / pl["Ox"], .93, .005); check("PLFA LA/Ox", pl["LA"] / pl["Ox"], .94, .005)
w = rd("window_comparisons_v16.csv").set_index(["end_h", "history"]).ratio_to_Ox
check("SA/Ox 1–12 h", w[(12, "SA")], .48, .005); check("SA/Ox 1–24 h", w[(24, "SA")], 1.19, .005)
co2 = rd("figure3_derived_co2_windows.csv").set_index(["soil", "window"]).sterilized_pct
for k, q in ((("Paddy soil", "0.5-6 h"), 43), (("Paddy soil", "0.5-48 h"), 20), (("Sandy soil", "0.5-6 h"), 91), (("Sandy soil", "0.5-48 h"), 49)): check(f"Liu sterilized {k}", co2[k], q, .5)
check("tau*", th.crossover_time(20, 1, 60, .02), 20.3, .05)
for L, q in ((0, 1.09), (4, 1.73), (6, 2.44)): check(f"E({L})", float(th.escape(L)), q, .005)
check("B(R=1)", float(th.carried_legacy(1)), .58, .005)
out = pd.DataFrame(rows); out.to_csv(os.path.join(D, "verification_reported_values.csv"), index=False)
print(out.to_string(index=False)); print("ALL PASS" if out.pass_.all() else "FAILURES")
