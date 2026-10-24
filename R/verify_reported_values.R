# Independent base-R check of every number quoted in the manuscript text and figure legends.
# Mirrors python/verify_reported_values.py. Run from R/:  Rscript verify_reported_values.R
# Writes ../data/verification_reported_values_R.csv and stops with an error if any check fails.
source("theory_4R.R")
D <- file.path("..", "data"); rd <- function(f) read.csv(file.path(D, f), check.names = FALSE)
rows <- list()
check <- function(name, value, quoted, tol) {
  ok <- abs(value - quoted) <= tol
  cat(sprintf("%-6s %-40s %10.4f  quoted %s\n", if (ok) "OK" else "FAIL", name, value, format(quoted)))
  rows[[length(rows) + 1]] <<- data.frame(check = name, computed = value, quoted = quoted, tolerance = tol, pass = ok)
}
# Li et al. 2025: assay capacity, Fe fractions, goethite, diel changes
eec <- rd("li_2025_eec_primary_source.csv"); e <- setNames(eec$eec, eec$compartment)
check("EEC rhizosphere/bulk", e[["Rhizosphere"]] / e[["Bulk soil"]], 1.25, .005)
check("EEC plaque/bulk", e[["Iron plaque"]] / e[["Bulk soil"]], 1.92, .005)
check("EEC bulk/plaque", e[["Bulk soil"]] / e[["Iron plaque"]], .52, .005)
fr <- rd("li_2025_fe_fractions_primary_source.csv"); fr <- fr[fr$fraction_index == 2, ]; pc <- setNames(fr$percent_of_total_Fe, fr$compartment)
check("poorly crystalline Fe, plaque %", pc[["Iron plaque"]], 75, .5)
check("poorly crystalline Fe, rhizosphere %", pc[["Rhizosphere"]], 30, 3)
check("poorly crystalline Fe, bulk %", pc[["Bulk soil"]], 30, 3)
g <- rd("li_2025_goethite_eec_primary_source.csv"); gg <- setNames(g$EEC, g$goethite)
check("goethite cycled/untreated", gg[["after 30 d diel redox cycling"]] / gg[["without redox cycling"]], 1.62, .005)
diel <- rd("li_2025_diel_feII_P_primary_source.csv")
sw <- function(v, comp) mean(abs(diff(diel[diel$compartment == comp, v])))
for (a in list(c("FeII", "Rhizosphere", .40), c("FeII", "Bulk soil", .08), c("P", "Rhizosphere", .52), c("P", "Bulk soil", .14)))
  check(paste(a[1], "change", a[2], "/ plaque"), sw(a[1], a[2]) / sw(a[1], "Iron plaque"), as.numeric(a[3]), .005)
# Patzner et al. 2020: reactive-Fe horizon shares
fe <- rd("patzner_panel_p_mineral_fe_oc_measured.csv"); fe <- fe[fe$metric == "Reactive Fe", ]
sh <- function(st, hz) 100 * fe$value[fe$stage == st & fe$horizon == hz] / sum(fe$value[fe$stage == st])
check("palsa mineral horizon %", sh("Palsa", "Mineral horizon"), 91, .5)
check("bog transition zone %", sh("Bog", "Transition zone"), 66, .5)
# Sennett et al. 2024: half-times, correlations, ratios
M <- rd("figure3_derived_metrics.csv")
check("N2 half-time min (h)", min(M$N2_half_time_h), 14, .5); check("N2 half-time max (h)", max(M$N2_half_time_h), 78, .5)
rho_full <- cor(M$N2_half_time_h, M$N2O_time_integral, method = "spearman")
rho_124 <- cor(M$N2_half_time_h, M$N2O_integral_1_24h, method = "spearman")
n <- nrow(M); tstat <- rho_124 * sqrt((n - 2) / (1 - rho_124^2)); p_124 <- 2 * pt(-abs(tstat), n - 2)
check("Spearman rho, complete record", rho_full, .98, .005)
check("Spearman rho, 1-24 h", rho_124, .37, .005); check("P, 1-24 h (t approximation)", p_124, .33, .005)
check("1-24 h integral min", min(M$N2O_integral_1_24h), 45, .5); check("1-24 h integral max", max(M$N2O_integral_1_24h), 79, .5)
mh <- aggregate(M[, c("N2_final", "N2O_integral_1_24h")], list(history = M$history), mean); rownames(mh) <- mh$history
check("terminal N2 SA/Ox", mh["SA", "N2_final"] / mh["Ox", "N2_final"], .91, .005)
check("terminal N2 LA/Ox", mh["LA", "N2_final"] / mh["Ox", "N2_final"], .86, .005)
check("N2O 1-24 h SA/Ox", mh["SA", "N2O_integral_1_24h"] / mh["Ox", "N2O_integral_1_24h"], 1.19, .005)
check("N2O 1-24 h LA/Ox", mh["LA", "N2O_integral_1_24h"] / mh["Ox", "N2O_integral_1_24h"], 1.24, .005)
pl <- tapply(rd("sennett_plfa_total.csv")$total_plfa_pmol_g, rd("sennett_plfa_total.csv")$treatment, mean)
check("PLFA SA/Ox", pl[["SA"]] / pl[["Ox"]], .93, .005); check("PLFA LA/Ox", pl[["LA"]] / pl[["Ox"]], .94, .005)
w <- rd("window_comparisons_v16.csv")
check("SA/Ox, 1-12 h", w$ratio_to_Ox[w$end_h == 12 & w$history == "SA"], .48, .005)
check("SA/Ox, 1-24 h", w$ratio_to_Ox[w$end_h == 24 & w$history == "SA"], 1.19, .005)
# Liu et al. 2025: sterilized CO2 as % of untreated
co2 <- rd("figure3_derived_co2_windows.csv")
for (a in list(c("Paddy soil", "0.5-6 h", 43), c("Paddy soil", "0.5-48 h", 20), c("Sandy soil", "0.5-6 h", 91), c("Sandy soil", "0.5-48 h", 49)))
  check(paste("sterilized %", a[1], a[2]), co2$sterilized_pct[co2$soil == a[1] & co2$window == a[2]], as.numeric(a[3]), .5)
# Analytical results, equations (1)-(3)
check("crossing time tau* (h)", crossover_time(20, 1, 60, 0.02), 20.3, .05)
for (a in list(c(0, 1.09), c(4, 1.73), c(6, 2.44))) check(paste0("escape E(", a[1], " h) (umol N)"), escape(a[1]), a[2], .005)
check("carry-over B at R = 1", carried_legacy(1), .58, .005)
check("access confound: 40*0.5 vs 100*0.2", 40 * 0.5 - 100 * 0.2, 0, 1e-12)
out <- do.call(rbind, rows); write.csv(out, file.path(D, "verification_reported_values_R.csv"), row.names = FALSE)
if (!all(out$pass)) stop("Some checks failed: ", paste(out$check[!out$pass], collapse = "; "))
cat("\nALL PASS (", nrow(out), " checks)\n", sep = "")
