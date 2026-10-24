"""Closed-form results of the out-of-step recovery theory (Python port of R/theory_4R.R)."""
import numpy as np
from scipy.optimize import brentq

def c_acc(tau, Q, alpha, k):
    """Eq. (1): cumulative accessible exchange."""
    tau = np.atleast_1d(np.asarray(tau, float))
    Q, alpha, k = (np.atleast_1d(np.asarray(v, float)) for v in (Q, alpha, k))
    return np.array([np.sum(Q * alpha * (-np.expm1(-k * t))) for t in tau])

def crossover_time(Pa, ka, Pb, kb, tmax=1e6):
    f = lambda t: Pa * (1 - np.exp(-ka * t)) - Pb * (1 - np.exp(-kb * t))
    if np.sign(f(1e-9)) == np.sign(f(tmax)):
        return np.nan
    return brentq(f, 1e-9, tmax, xtol=1e-12)

def x_traj(t, L, p0=1.0, l=0.1, k=1.0, T=12.0):
    """Exact piecewise solution of dX/dt = p(t) - [l + k s(t)] X."""
    out = []
    cuts = sorted(set([0.0, L, T]))
    for tj in np.atleast_1d(t):
        X, t0 = 0.0, 0.0
        for cc in cuts[1:] + [np.inf]:
            seg_end = min(cc, tj)
            if seg_end > t0:
                p = p0 if t0 < T else 0.0
                r = l + (k if t0 >= L else 0.0)
                X = p / r + (X - p / r) * np.exp(-r * (seg_end - t0))
                t0 = seg_end
            if cc >= tj:
                break
        out.append(X)
    return np.array(out)

def lag_integral(L, p0=1.0, l=0.1, T=12.0):
    L = np.asarray(L, float); Lp = np.minimum(L, T)
    I1 = (p0 / l) * (Lp - (1 - np.exp(-l * Lp)) / l)
    XT = (p0 / l) * (1 - np.exp(-l * Lp))
    return I1 + np.where(L > T, XT * (1 - np.exp(-l * (L - T))) / l, 0.0)

def escape(L, p0=1.0, l=0.1, k=1.0, T=12.0):
    """Eq. (2): cumulative escape after the pool clears."""
    return l * p0 * T / (l + k) + (k / (l + k)) * l * lag_integral(L, p0, l, T)

def excess_escape(L, **kw):
    return escape(L, **kw) - escape(0.0, **kw)

def carried_legacy(R):
    """Eq. (3): B = 1/(exp(1/R) - 1)."""
    return 1.0 / np.expm1(1.0 / np.asarray(R, float))

def legacy_pre(n, delta, Tm, h=1.0):
    r = np.exp(-delta / Tm); n = np.arange(1, n + 1)
    return h * r * (1 - r ** (n - 1)) / (1 - r)

if __name__ == "__main__":
    print("escape L=0,4,6:", [round(float(escape(L)), 2) for L in (0, 4, 6)])
    print("excess 4,6:", [round(float(excess_escape(L)), 2) for L in (4, 6)])
    print("B(R=1):", round(float(carried_legacy(1.0)), 2))
    print("crossover:", round(crossover_time(20, 1.0, 60, 0.02), 1))
