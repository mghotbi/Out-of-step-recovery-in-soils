# Out-of-step recovery: closed-form results used in the manuscript.
# Capacity Q, connectivity alpha, kinetics k; memory is the update between events.
# Base R only. Every function is analytic; this file is the authoritative R implementation.

# ---------------------------------------------------------------- Eq. (1)
c_acc <- function(tau, Q, alpha, k) sapply(tau, function(t) sum(Q * alpha * (-expm1(-k * t))))
plateau <- function(Q, alpha) sum(Q * alpha)                 # capacity and connectivity set it jointly
half_time <- function(k) log(2) / k                          # single pool, independent of Q and alpha
aggregate_half_time <- function(Q, alpha, k, tmax = 1e6) {   # multipool: time to half the plateau
  target <- 0.5 * plateau(Q, alpha)
  uniroot(function(t) c_acc(t, Q, alpha, k) - target, c(0, tmax), tol = 1e-12)$root
}
amplitudes <- function(Q, alpha) Q * alpha                   # P1: access changes amplitudes, not rate constants

# ---------------------------------------------------------------- P2: crossover of two single pools
crossover_exists <- function(Pa, ka, Pb, kb) (Pa * ka - Pb * kb) * (Pa - Pb) < 0
crossover_time <- function(Pa, ka, Pb, kb, tmax = 1e6) {
  f <- function(t) Pa * (1 - exp(-ka * t)) - Pb * (1 - exp(-kb * t))
  if (sign(f(1e-9)) == sign(f(tmax))) return(NA_real_)
  uniroot(f, c(1e-9, tmax), tol = 1e-12)$root
}

# ---------------------------------------------------------------- Eq. (2): late terminal sink
# dX/dt = p(t) - [l + k s(t)] X ; p = p0 for t < T, else 0 ; s = 0 for t < L, else 1 ; escape flux = l X
x_traj <- function(t, L, p0 = 1, l = 0.1, k = 1, T = 12) {  # exact piecewise-analytic solution
  cuts <- sort(unique(c(0, L, T)))
  sapply(t, function(tj) {
    X <- 0; t0 <- 0
    for (cc in c(cuts[-1], Inf)) {
      seg_end <- min(cc, tj)
      if (seg_end > t0) {
        p <- if (t0 < T) p0 else 0
        r <- l + (if (t0 >= L) k else 0)
        X <- if(r==0) X + p*(seg_end-t0) else p / r + (X - p / r) * exp(-r * (seg_end - t0))
        t0 <- seg_end
      }
      if (cc >= tj) break
    }
    X
  })
}
lag_integral <- function(L, p0 = 1, l = 0.1, T = 12) {
  Lp <- pmin(L, T)
  if(l==0) return(p0*Lp^2/2 + pmax(L-T,0)*p0*T)
  I1 <- (p0 / l) * (Lp - (1 - exp(-l * Lp)) / l)
  X_T <- (p0 / l) * (1 - exp(-l * Lp))
  I1 + ifelse(L > T, X_T * (1 - exp(-l * (L - T))) / l, 0)
}
escape <- function(L, p0 = 1, l = 0.1, k = 1, T = 12)        # exact Eq. (2)
  l * p0 * T / (l + k) + (k / (l + k)) * l * lag_integral(L, p0, l, T)
redox_debt <- function(L, p0 = 1, l = 0.1, k = 1, T = 12) escape(L, p0, l, k, T) - escape(0, p0, l, k, T)
redox_debt_small_lag <- function(L, p0 = 1, l = 0.1, k = 1) (k / (k + l)) * l * p0 * L^2 / 2  # overestimates
debt_ceiling <- function(p0 = 1, l = 0.1, k = 1, T = 12) k * p0 * T / (k + l)

# ---------------------------------------------------------------- Eq. (3): recurrence (recall)
# carrier relaxes as exp(-t/Tm) between events; each event adds a constant increment h
legacy_pre <- function(n_events, delta, Tm, h = 1) {         # sampled immediately before each event
  r <- exp(-delta / Tm); n <- seq_len(n_events)
  h * r * (1 - r^(n - 1)) / (1 - r)
}
carried_legacy <- function(R) 1 / expm1(1 / R)               # B = 1/(exp(1/R) - 1), R = Tm/delta
post_event_legacy <- function(R) 1 + carried_legacy(R)       # A = 1/(1 - exp(-1/R))
events_to_95pct <- function(R) R * log(20)
