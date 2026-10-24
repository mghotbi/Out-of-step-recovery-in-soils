# Verify Supplementary Figure 2 source data in two stages; stops with an error on any failure. Base R only.
# Stage 1  primary workbook -> extracted tables (Supplementary Figure 2a-f; stored cell values, no recalculation)
# Stage 2  extracted tables -> ../data/source_data_supplementary_figure2.csv (panel-level source data)
# The small OOXML reader below handles the cells used here; it is not a general Excel parser.
# Run from the R/ folder:  source("verify_primary_tables.R")
options(encoding = "UTF-8")
DATA <- file.path("..", "data"); WB <- file.path(DATA, "primary_source_workbooks")
xml_text <- function(file, member) {
  info <- unzip(file, list = TRUE); n <- info$Length[match(member, info$Name)]
  con <- unz(file, member, "rb"); on.exit(close(con)); rawToChar(readBin(con, "raw", n = n))
}
unescape <- function(x) { for (p in list(c("&amp;", "&"), c("&lt;", "<"), c("&gt;", ">"), c("&quot;", '"'), c("&apos;", "'"))) x <- gsub(p[1], p[2], x, fixed = TRUE); x }
get_sheet <- function(file, sheet) {
  file <- file.path(WB, file); members <- unzip(file, list = TRUE)$Name; ss <- character()
  if ("xl/sharedStrings.xml" %in% members) {
    z <- xml_text(file, "xl/sharedStrings.xml"); v <- regmatches(z, gregexpr("<si[ >][\\s\\S]*?</si>", z, perl = TRUE))[[1]]
    ss <- vapply(v, function(s) { t <- regmatches(s, gregexpr("<t(?: [^>]*)?>[\\s\\S]*?</t>", s, perl = TRUE))[[1]]
      unescape(paste(gsub("<[^>]+>", "", t), collapse = "")) }, character(1))
  }
  z <- xml_text(file, paste0("xl/worksheets/sheet", sheet, ".xml")); z <- gsub("<c\\b[^>]*/>", "", z, perl = TRUE)
  tags <- regmatches(z, gregexpr('<c\\b[^>]*r="[A-Z]+[0-9]+"[^>]*>[\\s\\S]*?</c>', z, perl = TRUE))[[1]]
  id <- sub('^[\\s\\S]*?\\br="([A-Z]+[0-9]+)"[\\s\\S]*$', "\\1", tags, perl = TRUE)
  val <- vapply(tags, function(t) { if (!grepl("<v>", t, fixed = TRUE)) return(NA_character_)
    v <- sub("^[\\s\\S]*?<v>(.*?)</v>[\\s\\S]*$", "\\1", t, perl = TRUE); if (grepl('t="s"', t, fixed = TRUE)) ss[as.integer(v) + 1L] else v }, character(1))
  setNames(val, id)
}
colname <- function(n) { v <- ""; while (n > 0) { r <- (n - 1) %% 26; v <- paste0(LETTERS[r + 1], v); n <- (n - 1) %/% 26 }; v }
num <- function(s, cells) suppressWarnings(as.numeric(s[cells]))
rd <- function(f, ...) read.csv(file.path(DATA, f), stringsAsFactors = FALSE, ...)
checks <- data.frame(stage = integer(), check = character(), n = integer(), max_abs_difference = numeric(), ok = logical())
record <- function(stage, label, a, b, tol = 1e-8) {
  a <- as.numeric(a); b <- as.numeric(b)
  ok <- length(a) == length(b) && all(is.finite(a)) && all(is.finite(b))
  err <- if (ok && length(a)) max(abs(a - b)) else Inf; ok <- ok && err < tol
  checks[nrow(checks) + 1, ] <<- list(stage, label, length(a), err, ok)
}

# ---------------- Stage 1 ----------------
s <- get_sheet("Li_2025_source.xlsx", 2); d <- rd("li_2025_eec_primary_source.csv")
record(1, "a Li EEC means", d$eec, num(s, paste0("AA", 5:7))); record(1, "a Li EEC SD", d$sd, num(s, paste0("AB", 5:7)))
s <- get_sheet("Li_2025_source.xlsx", 1); d <- rd("li_2025_o2_eh_ph_primary_source.csv")
for (j in 1:3) {
  a <- c("A", "F", "J")[j]; b <- c("B", "G", "K")[j]; v <- c("O2", "Eh", "pH")[j]; rn <- 5:2169
  xx <- num(s, paste0(a, rn)); yy <- num(s, paste0(b, rn)); k <- is.finite(xx) & is.finite(yy); o <- order(xx[k])
  q <- d[d$variable == v, ]; q <- q[order(q$time), ]
  record(1, paste("i Li", v, "times"), q$time, xx[k][o]); record(1, paste("i Li", v, "values"), q$value, yy[k][o])
}
s <- get_sheet("Patzner_2020_stock_source.xlsx", 1); d <- rd("patzner_panel_p_mineral_fe_oc_measured.csv"); d <- d[d$metric == "Reactive Fe", ]
hz <- c("Organic horizon", "Transition zone", "Mineral horizon")
for (j in 1:3) {
  st <- c("Palsa", "Bog", "Fen")[j]; rn <- list(8:10, 19:21, 27:29)[[j]]; q <- d[d$stage == st, ]; q <- q[match(hz, q$horizon), ]
  record(1, paste("b Patzner", st, "stocks"), q$value, num(s, paste0("E", rn))); record(1, paste("b Patzner", st, "errors"), q$error, num(s, paste0("F", rn)))
}
s <- get_sheet("Liu_2025_source.xlsx", 3); d <- rd("liu_2025_panel_h_co2_efflux.csv")
for (j in 1:2) for (k in 1:3) {
  soil <- c("Sandy soil", "Paddy soil")[j]; tr <- c("Original soil", "Sterilized soil", "Soil with TBA")[k]; rn <- list(5:12, 19:26)[[j]]
  q <- d[d$soil == soil & d$treatment == tr, ]; q <- q[order(q$time), ]
  record(1, paste("e Liu", soil, tr, "times"), q$time, num(s, paste0("B", rn)))
  record(1, paste("e Liu", soil, tr, "means"), q$value, num(s, paste0(c("C", "E", "G")[k], rn)))
  record(1, paste("e Liu", soil, tr, "SD"), q$error, num(s, paste0(c("D", "F", "H")[k], rn)))
}
gs <- rd("sennett_final_incubation_gas.csv")
for (k in 1:3) {
  h <- c("Ox", "SA", "LA")[k]; s <- get_sheet("Sennett_2024_source.xlsx", 3 + k); q <- gs[gs$history == h, ]; rn <- 4 + seq_len(nrow(q))
  record(1, paste("d Sennett", h, "final-incubation times"), q$time_h, num(s, paste0("A", rn)))
  record(1, paste("d Sennett", h, "N2O"), q$N2O, num(s, paste0("D", rn)))
  record(1, paste("d Sennett", h, "cumulative N2"), q$N2_cum, num(s, paste0("E", rn)))
  record(1, paste("d Sennett", h, "replicate labels"), as.numeric(identical(q$replicate, unname(s[paste0("B", rn)]))), 1)
}
record(1, "d Sennett row counts (Ox 51, SA 303, LA 164)", as.numeric(table(gs$history)[c("Ox", "SA", "LA")]), c(51, 303, 164), tol = 0.5)
s <- get_sheet("Sennett_2024_source.xlsx", 7); pl <- rd("sennett_plfa_total.csv"); rn <- 3 + seq_len(nrow(pl))
record(1, "e Sennett total PLFA", pl$total_plfa_pmol_g, num(s, paste0("B", rn)))
record(1, "e Sennett PLFA treatment labels", as.numeric(identical(pl$treatment, unname(s[paste0("A", rn)]))), 1)

# ---------------- Stage 2 ----------------
# Independent recomputation of descriptive Supplementary Figure 2 quantities from the extracted tables.
S <- read.csv(file.path(DATA, "source_data_supplementary_figure2.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
P <- function(k, f = NULL) { g <- S[S$panel == k, ]; if (!is.null(f)) g <- g[g$facet == f, ]; g }
tz <- function(x, y) sum(diff(x) * (y[-1] + y[-length(y)]) / 2)
key <- paste(gs$history, gs$replicate); bb <- P("b")
for (k in unique(key)) { q <- gs[key == k, ]; q <- q[order(q$time_h), ]; f <- max(q$N2_cum[length(q$N2_cum)])
  t50 <- approx(q$N2_cum[!duplicated(q$N2_cum)], q$time_h[!duplicated(q$N2_cum)], 0.5 * f, ties = "ordered")$y
  record(2, paste("b", k, "time-integrated N2O"), bb$y[bb$series == k], tz(q$time_h, q$N2O), tol = 1e-6)
  record(2, paste("b", k, "N2 half-time"), bb$x[bb$series == k], t50, tol = 1e-3) }
record(2, "c N2O series", P("c")$y, gs$N2O)
pl <- rd("sennett_plfa_total.csv"); record(2, "e PLFA (nmol g-1)", P("e")$y, pl$total_plfa_pmol_g / 1000)
l <- rd("liu_2025_panel_h_co2_efflux.csv"); ff <- P("f")
for (soil in c("Paddy soil", "Sandy soil")) for (w in list(c(0.5, 6), c(0.5, 48))) {
  cum <- function(tr) { q <- l[l$soil == soil & l$treatment == tr & l$time >= w[1] & l$time <= w[2], ]; q <- q[order(q$time), ]; tz(q$time, q$value) }
  lab <- sprintf("%g-%g h", w[1], w[2])
  record(2, paste("f", soil, lab, "sterilized %"), ff$y[ff$facet == soil & ff$series == paste("Sterilized", lab)], 100 * cum("Sterilized soil") / cum("Original soil"), tol = 1e-6)
  record(2, paste("f", soil, lab, "TBA %"), ff$y[ff$facet == soil & ff$series == paste("TBA", lab)], 100 * cum("Soil with TBA") / cum("Original soil"), tol = 1e-6) }
e <- rd("li_2025_eec_primary_source.csv"); gg <- P("g","EEC mean"); gg <- gg[match(e$compartment, gg$series), ]
record(2, "g EEC absolute means", gg$y, e$eec)
gsd<-P("g","EEC SD");record(2,"g EEC SD",gsd$y[match(e$compartment,gsd$series)],e$sd)
p <- rd("patzner_panel_p_mineral_fe_oc_measured.csv"); p <- p[p$metric == "Reactive Fe", ]; hh <- P("h")
for (st in unique(p$stage)) { q <- p[p$stage == st, ]; h2 <- hh[hh$facet == st, ]; h2 <- h2[match(q$horizon, h2$series), ]
  record(2, paste("h", st, "horizon shares"), h2$y, 100 * q$value / sum(q$value)) }

# Segment-wise clipping integrates each original linear segment, independently of the production helper.
clip_area<-function(x,y,a,b) {
 ans<-0
 for(i in seq_len(length(x)-1)) {
  lo<-max(a,x[i]);hi<-min(b,x[i+1]);if(hi<=lo)next
  slope<-(y[i+1]-y[i])/(x[i+1]-x[i]);zlo<-y[i]+slope*(lo-x[i]);zhi<-y[i]+slope*(hi-x[i])
  ans<-ans+(hi-lo)*(zlo+zhi)/2
 }
 ans
}
DM<-rd("figure3_derived_metrics.csv")
for(k in unique(key)) {
 q<-gs[key==k,];q<-q[order(q$time_h),];d<-DM[paste(DM$history,DM$replicate)==k,]
 record(2,paste("d",k,"common 1-24 h integral"),d$N2O_integral_1_24h,clip_area(q$time_h,q$N2O,1,24),tol=1e-6)
 record(2,paste("d",k,"sampling bounds and n"),c(d$first_time_h,d$last_time_h,d$n_observations),c(min(q$time_h),max(q$time_h),nrow(q)))
}
for(pair in list(c("Final N2 / Ox","N2_final"),c("Time-integrated N2O / Ox","N2O_time_integral"),c("1-24 h N2O integral / Ox","N2O_integral_1_24h"))) {
 d<-P("d",pair[1]);means<-tapply(DM[[pair[2]]],DM$history,mean)
 record(2,paste("d",pair[1]),d$y,unname(means[d$series]/means["Ox"]))
}
write.csv(checks, file.path(DATA, "primary_workbook_verification_R.csv"), row.names = FALSE)
for (i in seq_len(nrow(checks))) cat(sprintf("%-8s stage %d  %-42s n=%5d  max|diff|=%.2e\n", if (checks$ok[i]) "OK" else "MISMATCH",
                                             checks$stage[i], checks$check[i], checks$n[i], checks$max_abs_difference[i]))
cat(sprintf("\n%d comparisons, %d numeric entries; %d failed.\n", nrow(checks), sum(checks$n), sum(!checks$ok)))
cat("Workbook inputs, descriptive panel metrics, integration windows and EEC uncertainty verified; conditional theory is checked separately.\n")
if (any(!checks$ok)) stop("primary-table verification failed: ", paste(checks$check[!checks$ok], collapse = "; "))
