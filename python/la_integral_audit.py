"""Audit of the long-anoxic (LA) complete-record N2O integral against Table 1 of Sennett et al. (2024).
Run from python/:  python la_integral_audit.py
Prints the complete-record integrals for each history and the integration end time at which the LA mean equals the published 134."""
import os, numpy as np, pandas as pd
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
g = pd.read_csv(os.path.join(D, "sennett_final_incubation_gas.csv"))
def integral(q, t_end=None):
    q = q.sort_values("time_h"); t, y = q.time_h.values, q.N2O.values
    if t_end is None: return np.trapezoid(y, t)
    keep = t <= t_end
    return np.trapezoid(np.r_[y[keep], np.interp(t_end, t, y)], np.r_[t[keep], t_end])
published = {"Ox": 49.9, "SA": 592.0, "LA": 134.0}
for h in ("Ox", "SA", "LA"):
    reps = g[g.history == h].groupby("replicate")
    print(f"{h}: complete-record mean {np.mean([integral(q) for _, q in reps]):.3f} (published {published[h]}), "
          f"record ends {[round(q.time_h.max(), 1) for _, q in reps]} h")
la = [q for _, q in g[g.history == "LA"].groupby("replicate")]
lo, hi = 24.0, 70.0
for _ in range(80):
    mid = (lo + hi) / 2
    lo, hi = (mid, hi) if np.mean([integral(q, mid) for q in la]) < 134 else (lo, mid)
print(f"LA mean equals 134 when integrated to {mid:.1f} h; a difference confined to the endpoint after 24 h would not alter the separately calculated 1-24 h integrals.")
