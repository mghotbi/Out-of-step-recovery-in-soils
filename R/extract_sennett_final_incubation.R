# Extract final common anoxic incubation gas kinetics (Sennett et al. 2024, source data Fig. 3a-c) and
# total PLFA biomass (source data Fig. 4a) for Figure 3d,e. Base R only; reads stored workbook values.
# Run from the R/ folder:  source("extract_sennett_final_incubation.R")
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
num <- function(s, cells) suppressWarnings(as.numeric(s[cells]))
gas <- list()
for (k in 1:3) {
  hist <- c("Ox", "SA", "LA")[k]; s <- get_sheet("Sennett_2024_source.xlsx", 3 + k)
  stopifnot(identical(unname(s["A3"]), paste0(hist, "_gas")))
  rn <- 5:400; tt <- num(s, paste0("A", rn)); keep <- is.finite(tt); rn <- rn[keep]
  gas[[k]] <- data.frame(history = hist, replicate = unname(s[paste0("B", rn)]), time_h = num(s, paste0("A", rn)),
                         NO = num(s, paste0("C", rn)), N2O = num(s, paste0("D", rn)), N2_cum = num(s, paste0("E", rn)))
}
gas <- do.call(rbind, gas)
s <- get_sheet("Sennett_2024_source.xlsx", 7); rn <- 4:40; tr <- unname(s[paste0("A", rn)]); keep <- !is.na(tr)
plfa <- data.frame(treatment = tr[keep], total_plfa_pmol_g = num(s, paste0("B", rn[keep])))
write.csv(gas, file.path(DATA, "sennett_final_incubation_gas.csv"), row.names = FALSE)
write.csv(plfa, file.path(DATA, "sennett_plfa_total.csv"), row.names = FALSE)
cat(sprintf("Extracted %d gas observations (%s) and %d PLFA values.\n", nrow(gas),
            paste(names(table(gas$history)), table(gas$history), sep = "=", collapse = ", "), nrow(plfa)))
