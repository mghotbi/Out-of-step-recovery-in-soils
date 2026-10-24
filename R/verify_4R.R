# Verify every number quoted in the manuscript and Supplementary Note 1.
# Stops with an error if any check fails.  Run from the R/ folder:  source("verify_4R.R")
source("theory_4R.R")
failed <- character(0); rows <- list()
check <- function(name, value, expected, tol) {
  ok <- abs(value - expected) <= tol
  cat(sprintf("%-8s %-48s %12.5f  (expected %s)\n", if (ok) "OK" else "MISMATCH", name, value, format(expected)))
  rows[[length(rows) + 1]] <<- data.frame(check = name, value = value, expected = expected, ok = ok)
  if (!ok) failed <<- c(failed, name)
}

# Box 1, Eq. (1)
A <- c_acc(48, 120 * c(0.2, 0.8), c(0.25, 0.25), c(1, 0.01))
B <- c_acc(48, 120 * c(0.7, 0.3), c(0.70, 0.70), c(1, 0.01))
check("SI example: soil A, 48 h (mmol e-/kg)", A, 15.15, 0.005)
check("SI example: soil B, 48 h (mmol e-/kg)", B, 68.41, 0.005)
check("SI example: ratio B/A", B / A, 4.52, 0.005)

# P1
Q <- c(50, 50); k <- c(1, 0.01)
t_ref <- aggregate_half_time(Q, c(1, 1), k)
check("P1: two-pool half-time, alpha=(1,1) (h)", t_ref, 3.40, 0.005)
check("P1: proportional access, half-time unchanged", aggregate_half_time(Q, c(0.5, 0.5), k) - t_ref, 0, 1e-6)
check("P1: selective access alpha=(1,0.1) (h)", aggregate_half_time(Q, c(1, 0.1), k), 0.80, 0.005)
tt <- seq(0, 500, length.out = 2001)
nrm <- function(a) c_acc(tt, Q, a, k) / plateau(Q, a)
check("P1: proportional change, max |normalized diff|", max(abs(nrm(c(1, 1)) - nrm(c(0.5, 0.5)))), 0, 1e-12)

# P2
check("P2: crossing criterion satisfied (1=yes)", as.numeric(crossover_exists(20, 1, 60, 0.02)), 1, 0)
check("P2: crossing time tau* (h)", crossover_time(20, 1, 60, 0.02), 20.27, 0.005)
set.seed(1); grid <- c(seq(1e-6, 1, length.out = 2000), exp(seq(0, log(1e5), length.out = 4000))); mism <- 0
for (i in 1:2000) {
  ab <- runif(2, 1, 100); kk <- 10^runif(2, -3, 1)
  d <- ab[1] * (1 - exp(-kk[1] * grid)) - ab[2] * (1 - exp(-kk[2] * grid))
  mism <- mism + (any(diff(sign(d)) != 0) != crossover_exists(ab[1], kk[1], ab[2], kk[2]))
}
check("P2: criterion vs brute force (2000 random pairs)", mism, 0, 0)

# P3, Eq. (2), against an independent RK4 integration
ode_escape <- function(L, p0 = 1, l = 0.1, kk = 1, T = 12, dt = 5e-3, tend = 150) {
  f <- function(t, X) (if (t < T) p0 else 0) - (l + (if (t >= L) kk else 0)) * X
  X <- 0; E <- 0; t <- 0
  while (t < tend - 1e-12) {
    k1 <- f(t, X); k2 <- f(t + dt/2, X + dt*k1/2); k3 <- f(t + dt/2, X + dt*k2/2); k4 <- f(t + dt, X + dt*k3)
    Xn <- X + dt * (k1 + 2*k2 + 2*k3 + k4) / 6
    E <- E + l * 0.5 * (X + Xn) * dt; X <- Xn; t <- t + dt
  }
  E
}
check("P3: escape, on-time sink (umol N)", escape(0), 1.09, 0.005)
check("P3: escape, 4-h lag (umol N)", escape(4), 1.73, 0.005)
check("Fig 3d: escape, 6-h lag / dormant sink (umol N)", escape(6), 2.44, 0.005)
check("Fig 3d: cumulative ratio dormant/active", escape(6) / escape(0), 2.24, 0.005)
check("P3: escape L=4, analytic minus ODE", escape(4) - ode_escape(4), 0, 5e-4)
check("P3: exact debt, 4-h lag (umol N)", redox_debt(4), 0.639, 0.0005)
check("P3: short-lag approximation, 4-h lag", redox_debt_small_lag(4), 0.727, 0.0005)
check("P3: approximation overestimate (%)", 100 * (redox_debt_small_lag(4) / redox_debt(4) - 1), 13.8, 0.05)
check("P3: debt ceiling, sink never returns", debt_ceiling(), 10.909, 0.0005)

# P4, Eq. (3)
check("P4: carried legacy B at R=1", carried_legacy(1), 0.582, 0.0005)
check("P4: post-event A at R=1", post_event_legacy(1), 1.582, 0.0005)
check("P4: B at R=10", carried_legacy(10), 9.508, 0.0005)
check("P4: B = 1 at R = 1/ln2", carried_legacy(1 / log(2)), 1, 1e-12)
s <- legacy_pre(400, 0.5, 1)
check("P4: iterated series converges to B (R=2)", s[400] - carried_legacy(2), 0, 1e-9)

out <- do.call(rbind, rows)
if (dir.exists(file.path("..", "data"))) write.csv(out, file.path("..", "data", "verification_4R_R.csv"), row.names = FALSE)
if (length(failed)) stop(length(failed), " check(s) failed: ", paste(failed, collapse = "; "))
cat("\nALL CHECKS PASSED\n")
