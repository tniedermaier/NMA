###############################################################################
# Niedermaier et al. - COPD mortality by GOLD grade
# MASTER ANALYSIS SCRIPT
#
# Purpose:
#   One reproducible script containing only analyses reported in the manuscript
#   and supplement: primary network meta-analysis, subgroup/sensitivity analyses,
#   leave-one-out analysis, network diagnostics/figures, and conversion of the
#   pooled HRs to age- and sex-specific mortality probabilities.
#
# IMPORTANT methodological correction relative to the current manuscript files:
#   Ogata 2025 reported Bonferroni-adjusted 99.17% CIs for all three included
#   pairwise GOLD estimates. Their SEs are therefore reconstructed from the
#   99.17% CI width before pooling. This changes some currently reported values.
#
# The primary inferential method follows the current manuscript:
#   random-effects weighted least squares + study-clustered CR2 standard errors
#   with Satterthwaite degrees of freedom. netmeta is used for tau^2, I^2,
#   network structure, funnel plots and inconsistency diagnostics.
###############################################################################

rm(list = ls(all.names = TRUE))
options(width = 140)

## ============================================================================
## 0. USER SETTINGS
## ============================================================================

# Put this script in the same directory as the input files, or change paths here.
setwd("path\\to\\files")

nma_file <- "Data_Masterthesis_final_v1_smoking_adj_neu.xlsx"
nma_sheet <- "neu"

# German abridged period life tables
life_2022_24_file <- "mortality_rates_germany.xlsx"
life_2023_25_file <- "mortality_rates_germany2.xlsx"
life_2017_19_file <- "mortality_rates_germany_2017_19.xlsx"

# Austrian smoothed period life table
life_austria_file <- "Austria_Geglaettete_Sterbetafeln_2016_2024.ods"
life_austria_sheet <- "2022-2024"

output_dir <- "analysis_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

required_packages <- c("readxl", "netmeta", "clubSandwich", "meta", "readODS")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required packages first: ", paste(missing_packages, collapse = ", "))
}

library(readxl)
library(netmeta)
library(clubSandwich)
library(meta)
library(readODS)

## ============================================================================
## 1. LOAD AND CHECK NMA DATA
## ============================================================================

if (!file.exists(nma_file)) stop("NMA input file not found: ", nma_file)

neu0 <- readxl::read_excel(nma_file, sheet = nma_sheet)
neu0 <- as.data.frame(neu0)
neu1 <- neu0[!grepl("relevant", neu0[, 1]), , drop = FALSE]
neu <- neu1[neu1$include == 1, , drop = FALSE]

neu$log_hr    <- as.numeric(neu$log_hr)
neu$se_log_hr <- as.numeric(neu$se_log_hr)
neu$treat1    <- gsub("\\s+", "", as.character(neu$treat1))
neu$treat2    <- gsub("\\s+", "", as.character(neu$treat2))
neu$PMID      <- as.character(neu$PMID)
neu$author    <- as.character(neu$author)
neu$year      <- as.numeric(neu$year)
neu$setting   <- as.character(neu$setting)
neu$follow_up <- as.numeric(neu$follow_up)

stopifnot(all(is.finite(neu$log_hr)), all(is.finite(neu$se_log_hr)))
stopifnot(!anyNA(neu$treat1), !anyNA(neu$treat2), !anyNA(neu$PMID))

# Natural-log audit against raw HR column, if present.
if ("HR" %in% names(neu)) {
  hr_num <- suppressWarnings(as.numeric(neu$HR))
  ok <- is.finite(hr_num) & hr_num > 0
  max_log_difference <- max(abs(neu$log_hr[ok] - log(hr_num[ok])))
  if (max_log_difference > 1e-10) {
    stop("log_hr is not equal to ln(HR); max absolute difference = ", max_log_difference)
  }
}

cat("Included rows:", nrow(neu), "\n")
cat("Unique studies (PMID):", length(unique(neu$PMID)), "\n")
cat("Treatments:", paste(sort(unique(c(neu$treat1, neu$treat2))), collapse = ", "), "\n")

## ============================================================================
## 2. OGATA 2025: RECONSTRUCT SEs FROM BONFERRONI-ADJUSTED 99.17% CIs
## ============================================================================

# Included Ogata estimates (all vs GOLD 4) as published:
# HR 0.13 (0.02-0.77), HR 0.12 (0.03-0.41), HR 0.34 (0.14-0.83).
# The published intervals are 99.17% CIs after Bonferroni adjustment for six
# pairwise comparisons, not ordinary 95% CIs.

ogata_source <- data.frame(
  HR = c(0.13, 0.12, 0.34),
  lower_99_17 = c(0.02, 0.03, 0.14),
  upper_99_17 = c(0.77, 0.41, 0.83)
)

CL_ogata <- 0.9917
z_ogata <- qnorm(1 - (1 - CL_ogata) / 2)
ogata_source$SE <- (
  log(ogata_source$upper_99_17) - log(ogata_source$lower_99_17)
) / (2 * z_ogata)

ogata_idx <- which(tolower(neu$author) == "ogata")
if (length(ogata_idx) != 3) {
  stop("Expected exactly 3 included Ogata estimates; found ", length(ogata_idx))
}

for (i in ogata_idx) {
  hr_i <- exp(neu$log_hr[i])
  j <- which.min(abs(log(hr_i) - log(ogata_source$HR)))
  if (abs(log(hr_i) - log(ogata_source$HR[j])) > 0.03) {
    stop("Could not safely match an Ogata row to HR 0.13/0.12/0.34")
  }
  neu$se_log_hr[i] <- ogata_source$SE[j]
}

cat("Ogata SEs after correction:\n")
print(neu[ogata_idx, intersect(c("author", "year", "PMID", "treat1", "treat2", "log_hr", "se_log_hr"), names(neu))])

## ============================================================================
## 3. CORE FUNCTIONS: RANDOM-EFFECTS WLS + STUDY-CLUSTERED CR2 INFERENCE
## ============================================================================

make_X <- function(t1, t2, ref) {
  trts <- sort(unique(c(t1, t2)))
  trts <- trts[trts != ref]
  X <- matrix(0, nrow = length(t1), ncol = length(trts))
  colnames(X) <- trts

  for (i in seq_len(nrow(X))) {
    a <- t1[i]
    b <- t2[i]
    if (!is.na(a) && a != ref) X[i, a] <- X[i, a] - 1
    if (!is.na(b) && b != ref) X[i, b] <- X[i, b] + 1
  }
  X
}

make_netmeta <- function(dat, ref = "normal") {
  netmeta::netmeta(
    TE = dat$log_hr,
    seTE = dat$se_log_hr,
    treat1 = dat$treat1,
    treat2 = dat$treat2,
    common = FALSE,
    random = TRUE,
    reference.group = ref,
    sm = "HR",
    data = dat
  )
}

tau2_from_netmeta <- function(dat, ref) {
  nm <- make_netmeta(dat, ref)
  if (!is.null(nm$tau)) return(as.numeric(nm$tau)^2)
  if (!is.null(nm$tau.random)) return(as.numeric(nm$tau.random)^2)
  stop("tau not found in netmeta object")
}

I2_percent <- function(dat, ref) {
  x <- as.numeric(summary(make_netmeta(dat, ref))$I2)
  if (!is.finite(x)) return(NA_real_)
  if (x <= 1) x <- 100 * x
  x
}

fit_CR2 <- function(dat, ref) {
  dat <- as.data.frame(dat)
  X <- make_X(dat$treat1, dat$treat2, ref)
  tau2 <- tau2_from_netmeta(dat, ref)

  y <- dat$log_hr
  w <- 1 / (dat$se_log_hr^2 + tau2)
  fit <- lm(y ~ 0 + ., data = data.frame(y = y, X), weights = w)

  Vcr <- clubSandwich::vcovCR(fit, cluster = dat$PMID, type = "CR2")
  ct <- clubSandwich::coef_test(fit, vcov = Vcr, test = "Satterthwaite")

  est <- ct$beta
  seR <- ct$SE
  df <- ct$df
  tcrit <- qt(0.975, df = df)

  data.frame(
    treatment = rownames(ct),
    logHR = est,
    SE_CR2 = seR,
    df = df,
    HR = exp(est),
    lower = exp(est - tcrit * seR),
    upper = exp(est + tcrit * seR),
    tau2 = tau2,
    row.names = NULL,
    check.names = FALSE
  )
}

fmt_num <- function(x, digits = 2) formatC(x, format = "f", digits = digits)
fmt_ci <- function(hr, lo, hi, digits = 2) {
  paste0(fmt_num(hr, digits), " (", fmt_num(lo, digits), "-", fmt_num(hi, digits), ")")
}

format_against_ref <- function(res, treatment_order, ref, digits = 2) {
  out <- rep(NA_character_, length(treatment_order))
  names(out) <- treatment_order

  for (g in treatment_order) {
    if (g == ref) {
      out[g] <- "Ref."
    } else {
      j <- match(g, gsub("\\s+", "", res$treatment))
      if (!is.na(j)) out[g] <- fmt_ci(res$HR[j], res$lower[j], res$upper[j], digits)
    }
  }
  out
}

split_setting <- function(dat) {
  setting_lower <- tolower(trimws(dat$setting))
  list(
    pop = dat[setting_lower %in% c("bev", "population", "pop"), , drop = FALSE],
    clin = dat[setting_lower %in% c("klinisch", "clinical", "clinic"), , drop = FALSE]
  )
}

append_I2_row <- function(tab, values) {
  r <- as.list(rep(NA_character_, ncol(tab)))
  names(r) <- names(tab)
  r[[1]] <- "Heterogeneity (I2)"
  for (nm in intersect(names(values), names(tab))) {
    r[[nm]] <- paste0(formatC(values[[nm]], format = "f", digits = 1), "%")
  }
  rbind(tab, as.data.frame(r, stringsAsFactors = FALSE, check.names = FALSE))
}

## ============================================================================
## 4. TABLE 3: PRIMARY NMA + SETTING SUBGROUPS + FOLLOW-UP SUBGROUP
## ============================================================================

all_order <- c("normal", "GOLD1", "GOLD2", "GOLD3", "GOLD4")
ss <- split_setting(neu)

fit_all <- fit_CR2(neu, "normal")
fit_pop <- fit_CR2(ss$pop, "normal")
fit_clin <- fit_CR2(ss$clin, "GOLD1")
fit_pop_ref1 <- fit_CR2(ss$pop, "GOLD1")

Table3_main <- data.frame(
  Grade = c("Normal", "GOLD 1", "GOLD 2", "GOLD 3", "GOLD 4"),
  `All studies; ref normal` = format_against_ref(fit_all, all_order, "normal"),
  `Population studies; ref normal` = format_against_ref(fit_pop, all_order, "normal"),
  `Clinical setting; ref GOLD1` = format_against_ref(fit_clin, all_order, "GOLD1"),
  `Population studies; ref GOLD1` = format_against_ref(fit_pop_ref1, all_order, "GOLD1"),
  check.names = FALSE,
  row.names = NULL
)

Table3_main <- append_I2_row(
  Table3_main,
  c(
    `All studies; ref normal` = I2_percent(neu, "normal"),
    `Population studies; ref normal` = I2_percent(ss$pop, "normal"),
    `Clinical setting; ref GOLD1` = I2_percent(ss$clin, "GOLD1"),
    `Population studies; ref GOLD1` = I2_percent(ss$pop, "GOLD1")
  )
)

cat("\nTABLE 3 - MAIN ANALYSIS\n")
print(Table3_main, row.names = FALSE)
write.csv2(Table3_main, file.path(output_dir, "Table3_main.csv"), row.names = FALSE)

# Follow-up split at the study-level median reported in the manuscript (9 years).
fu_dat <- neu[!is.na(neu$follow_up), , drop = FALSE]
fu_by_study <- aggregate(follow_up ~ PMID, data = fu_dat, FUN = function(x) unique(x)[1])
fu_median <- median(fu_by_study$follow_up, na.rm = TRUE)
if (!isTRUE(all.equal(as.numeric(fu_median), 9))) {
  warning("Median follow-up is ", fu_median, ", not 9 years")
}

fu_short <- fu_dat[fu_dat$follow_up < 9, , drop = FALSE]
fu_long <- fu_dat[fu_dat$follow_up >= 9, , drop = FALSE]
fit_fu_short <- fit_CR2(fu_short, "normal")
fit_fu_long <- fit_CR2(fu_long, "normal")
gold_order <- c("GOLD1", "GOLD2", "GOLD3", "GOLD4")

Table3_followup <- data.frame(
  Grade = c("GOLD 1", "GOLD 2", "GOLD 3", "GOLD 4"),
  `<9 years` = format_against_ref(fit_fu_short, gold_order, "normal"),
  `>=9 years` = format_against_ref(fit_fu_long, gold_order, "normal"),
  check.names = FALSE,
  row.names = NULL
)
Table3_followup <- append_I2_row(
  Table3_followup,
  c(
    `<9 years` = I2_percent(fu_short, "normal"),
    `>=9 years` = I2_percent(fu_long, "normal")
  )
)

cat("\nTABLE 3 - FOLLOW-UP SUBGROUP\n")
print(Table3_followup, row.names = FALSE)
write.csv2(Table3_followup, file.path(output_dir, "Table3_followup.csv"), row.names = FALSE)

## ============================================================================
## 5. e-TABLE 4: PRESPECIFIED SENSITIVITY ANALYSES
## ============================================================================

make_sensitivity_block <- function(dat, label) {
  s <- split_setting(dat)
  r_all <- fit_CR2(dat, "normal")
  r_pop <- fit_CR2(s$pop, "normal")
  r_clin <- fit_CR2(s$clin, "GOLD1")
  r_pop1 <- fit_CR2(s$pop, "GOLD1")

  tab <- data.frame(
    Sensitivity = label,
    Grade = c("Normal", "GOLD 1", "GOLD 2", "GOLD 3", "GOLD 4"),
    `All studies; ref normal` = format_against_ref(r_all, all_order, "normal"),
    `Population studies; ref normal` = format_against_ref(r_pop, all_order, "normal"),
    `Clinical setting; ref GOLD1` = format_against_ref(r_clin, all_order, "GOLD1"),
    `Population studies; ref GOLD1` = format_against_ref(r_pop1, all_order, "GOLD1"),
    check.names = FALSE,
    row.names = NULL
  )

  i2 <- data.frame(
    Sensitivity = label,
    Grade = "Heterogeneity (I2)",
    `All studies; ref normal` = paste0(fmt_num(I2_percent(dat, "normal"), 1), "%"),
    `Population studies; ref normal` = paste0(fmt_num(I2_percent(s$pop, "normal"), 1), "%"),
    `Clinical setting; ref GOLD1` = paste0(fmt_num(I2_percent(s$clin, "GOLD1"), 1), "%"),
    `Population studies; ref GOLD1` = paste0(fmt_num(I2_percent(s$pop, "GOLD1"), 1), "%"),
    check.names = FALSE,
    row.names = NULL
  )

  rbind(tab, i2)
}

if (!"include_SA" %in% names(neu)) stop("Column include_SA is missing")

sa_alt_articles <- neu[neu$include_SA == 1, , drop = FALSE]
sa_no_PA <- neu[!neu$author %in% c("Stavem", "He", "Guo", "Ke"), , drop = FALSE]
sa_no_meds <- neu[!neu$author %in% c("Hyung", "Ogata"), , drop = FALSE]

eTable4 <- rbind(
  make_sensitivity_block(
    sa_alt_articles,
    "Bhatta 2020 instead of Leivseth 2013; Zou 2025 instead of Mannino 2012"
  ),
  make_sensitivity_block(
    sa_no_PA,
    "Excluding studies adjusted for physical activity (Stavem, He, Guo, Ke)"
  ),
  make_sensitivity_block(
    sa_no_meds,
    "Excluding studies adjusted for medication (Hyung, Ogata)"
  )
)

cat("\ne-TABLE 4 - SENSITIVITY ANALYSES\n")
print(eTable4, row.names = FALSE)
write.csv2(eTable4, file.path(output_dir, "eTable4_sensitivity.csv"), row.names = FALSE)

## ============================================================================
## 6. e-TABLE 3: LEAVE-ONE-STUDY-OUT, USING THE SAME CR2 METHOD
## ============================================================================

study_meta <- unique(neu[, c("PMID", "author", "year"), drop = FALSE])
loo <- vector("list", length(unique(neu$PMID)))
studies <- unique(neu$PMID)

for (i in seq_along(studies)) {
  p <- studies[i]
  d <- neu[neu$PMID != p, , drop = FALSE]
  
  r <- fit_CR2(d, "normal")
  f <- format_against_ref(r, gold_order, "normal")
  i2 <- I2_percent(d, "normal")
  
  m <- study_meta[study_meta$PMID == p, , drop = FALSE]

  loo[[i]] <- data.frame(
    Study_left_out = paste(m$author[1], m$year[1]),
    PMID = p,
    `GOLD 1` = f["GOLD1"],
    `GOLD 2` = f["GOLD2"],
    `GOLD 3` = f["GOLD3"],
    `GOLD 4` = f["GOLD4"],
    `I2 (%)` = round(i2, 1),
    check.names = FALSE,
    row.names = NULL
  )
}

eTable3 <- do.call(rbind, loo)
year_sort <- suppressWarnings(as.numeric(sub(".* ([0-9]{4})$", "\\1", eTable3$Study_left_out)))
eTable3 <- eTable3[order(year_sort, eTable3$Study_left_out), , drop = FALSE]

cat("\ne-TABLE 3 - LEAVE-ONE-OUT\n")
print(eTable3, row.names = FALSE)
write.csv2(eTable3, file.path(output_dir, "eTable3_leave_one_out.csv"), row.names = FALSE)

## ============================================================================
## 7. NETWORK DIAGNOSTICS AND e-FIGURES 2-6
## ============================================================================

nm <- make_netmeta(neu, "normal")

# e-Figure 2: network graph
pdf(file.path(output_dir, "eFigure2_network.pdf"), width = 8, height = 7)
netgraph(nm)
dev.off()

# e-Figure 3: comparison-adjusted funnel plots.
# funnel.netmeta itself correctly displays Egger, Begg-Mazumdar and
# Thompson-Sharp p-values. For the tabular p-values below, use the adjusted
# data returned by funnel() and meta::metabias(), rather than metabias.netmeta.
ord_all <- c("normal", setdiff(nm$trts, "normal"))

pdf(file.path(output_dir, "eFigure3_funnel.pdf"), width = 8, height = 12)
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(3, 1))
funnel_gold1 <- funnel(
  nm, order = "GOLD1",
  method.bias = c("Egger", "Begg", "Thompson"), digits.pval = 3
)
funnel_normal <- funnel(
  nm, order = "normal",
  method.bias = c("Egger", "Begg", "Thompson"), digits.pval = 3
)
funnel_all <- funnel(
  nm, order = ord_all,
  method.bias = c("Egger", "Begg", "Thompson"), digits.pval = 3
)
par(oldpar)
dev.off()

bias_tests_from_funnel <- function(funnel_data, label) {
  if (!all(c("TE.adj", "seTE") %in% names(funnel_data))) {
    stop("funnel.netmeta return object lacks TE.adj/seTE")
  }
  # Match funnel.netmeta() exactly: its internal comparison-adjusted
  # meta-analysis uses DerSimonian-Laird tau^2 for the asymmetry tests.
  m <- meta::metagen(
    TE = TE.adj,
    seTE = seTE,
    data = funnel_data,
    common = FALSE,
    random = TRUE,
    sm = "HR",
    method.tau = "DL",
    method.tau.ci = ""
  )

  methods <- c(Egger = "Egger", `Begg-Mazumdar` = "Begg", `Thompson-Sharp` = "Thompson")
  out <- lapply(names(methods), function(lbl) {
    z <- suppressWarnings(meta::metabias(m, method = methods[[lbl]], k.min = 10))
    data.frame(
      Funnel = label,
      Test = lbl,
      k = nrow(funnel_data),
      p_value = z$p.value,
      row.names = NULL
    )
  })
  do.call(rbind, out)
}

funnel_tests <- rbind(
  bias_tests_from_funnel(funnel_gold1, "GOLD1 comparator"),
  bias_tests_from_funnel(funnel_normal, "Normal comparator"),
  bias_tests_from_funnel(funnel_all, "All treatments / comparison-adjusted")
)
cat("\nFUNNEL-PLOT ASYMMETRY TESTS\n")
print(funnel_tests, row.names = FALSE)
write.csv2(funnel_tests, file.path(output_dir, "eFigure3_asymmetry_tests.csv"), row.names = FALSE)

# e-Figure 4: netheat inconsistency diagnostic
pdf(file.path(output_dir, "eFigure4_netheat.pdf"), width = 10, height = 10)
netheat(nm)
dev.off()

# e-Figures 5 and 6: direct-comparison meta-analyses used only as approximate
# study-weight displays. They are not the NMA estimates.
neu_plot <- neu[order(neu$year, neu$author), , drop = FALSE]
neu_plot$studlab <- paste(neu_plot$author, neu_plot$year)

make_direct_forest <- function(dat, reference, file) {
  d <- dat[dat$treat1 == reference, , drop = FALSE]
  if (!nrow(d)) stop("No direct comparisons with treat1 = ", reference)

  m <- meta::metagen(
    TE = log_hr,
    seTE = se_log_hr,
    studlab = studlab,
    data = d,
    common = FALSE,
    random = TRUE,
    method.tau = "REML",
    sm = "HR"
  )
  m_sub <- update(m, subgroup = treat2)

  # Let forest.meta determine the required height automatically.
  # The figures are used only to show the approximate direct-comparison
  # study weights, so omit heterogeneity and subgroup-difference text.
  meta::forest(
    m_sub,
    overall = FALSE,
    leftcols = "studlab",
    leftlabs = "Study",
    rightcols = c("effect", "ci", "w.random"),
    rightlabs = c("HR", "95% CI", "Weight"),
    hetstat = FALSE,
    subgroup.hetstat = FALSE,
    test.subgroup = FALSE,
    print.tau2 = FALSE,
    filename = file,
    width = 9
  )
  invisible(m_sub)
}

forest_normal <- make_direct_forest(
  neu_plot, "normal", file.path(output_dir, "eFigure5_forest_ref_normal.pdf")
)
forest_gold1 <- make_direct_forest(
  neu_plot, "GOLD1", file.path(output_dir, "eFigure6_forest_ref_GOLD1.pdf")
)

## ============================================================================
## 8. MORTALITY-PROBABILITY FUNCTIONS (TABLE 4, e-TABLES 5-7)
## ============================================================================

# The pooled all-study HRs and CR2 CIs from the primary analysis are the model
# inputs for all country-specific mortality calculations.

extract_hr_vector <- function(fit) {
  target <- c("GOLD1", "GOLD2", "GOLD3", "GOLD4")
  idx <- match(target, gsub("\\s+", "", fit$treatment))
  if (anyNA(idx)) stop("Could not find all GOLD1-4 estimates")
  list(
    HR = setNames(fit$HR[idx], target),
    LCI = setNames(fit$lower[idx], target),
    UCI = setNames(fit$upper[idx], target)
  )
}

hr_input <- extract_hr_vector(fit_all)
HR <- hr_input$HR
LCI <- hr_input$LCI
UCI <- hr_input$UCI

age_groups_keep <- c(
  "45 - 49", "50 - 54", "55 - 59", "60 - 64", "65 - 69",
  "70 - 74", "75 - 79", "80 - 84", "85+"
)

german_interval_n <- function(label, n100 = 1) {
  label <- trimws(as.character(label))
  if (grepl("^\\d+\\s*-\\s*\\d+$", label)) {
    x <- as.numeric(strsplit(gsub("\\s+", "", label), "-")[[1]])
    return(x[2] - x[1] + 1)
  }
  if (grepl("^\\d+\\+$", label)) return(n100)
  NA_real_
}

read_german_lifetable <- function(file, sheet) {
  if (!file.exists(file)) stop("Life-table file not found: ", file)
  x <- as.data.frame(readxl::read_excel(file, sheet = sheet))

  if (!"altersgruppe" %in% names(x)) {
    age_guess <- grep("alter", names(x), ignore.case = TRUE)
    if (length(age_guess)) names(x)[age_guess[1]] <- "altersgruppe" else names(x)[5] <- "altersgruppe"
  }

  q_guess <- grep("sterbewkeit|sterbewahrscheinlichkeit|^qx$", names(x), ignore.case = TRUE)
  if (!length(q_guess)) stop("Could not identify nqx column in ", file, " / ", sheet)
  if (names(x)[q_guess[1]] != "sterbewkeit") names(x)[q_guess[1]] <- "sterbewkeit"

  if (!"ueberlebende" %in% names(x)) stop("Column 'ueberlebende' missing in ", file, " / ", sheet)
  if (!"gestorbene" %in% names(x)) stop("Column 'gestorbene' missing in ", file, " / ", sheet)

  x$altersgruppe <- trimws(as.character(x$altersgruppe))
  x$sterbewkeit <- as.numeric(x$sterbewkeit)
  x$ueberlebende <- as.numeric(x$ueberlebende)
  x$gestorbene <- as.numeric(x$gestorbene)
  x
}

base_hazards_germany <- function(mort_table, ax = 0.5, n100 = 1, eps = 1e-12) {
  components_85 <- c("85 - 89", "90 - 94", "95 - 99", "100+")
  regular <- age_groups_keep[age_groups_keep != "85+"]
  out <- vector("list", length(regular) + 1)

  for (i in seq_along(regular)) {
    ag <- regular[i]
    row <- mort_table[mort_table$altersgruppe == ag, , drop = FALSE]
    if (nrow(row) != 1) stop("Age group not found exactly once: ", ag)
    n <- german_interval_n(ag, n100)
    qn <- min(as.numeric(row$sterbewkeit), 1 - eps)
    h <- -log1p(-qn) / n
    out[[i]] <- data.frame(age_group = ag, base_hazard = h)
  }

  sub <- mort_table[mort_table$altersgruppe %in% components_85, , drop = FALSE]
  sub <- sub[match(components_85, sub$altersgruppe), , drop = FALSE]
  if (nrow(sub) != length(components_85) || anyNA(sub$altersgruppe)) {
    stop("Could not find all 85+ components: ", paste(components_85, collapse = ", "))
  }

  n_vec <- vapply(components_85, german_interval_n, numeric(1), n100 = n100)
  qn <- pmin(sub$sterbewkeit, 1 - eps)
  h <- -log1p(-qn) / n_vec
  L <- n_vec * (sub$ueberlebende - ax * sub$gestorbene)
  if (any(L <= 0)) stop("Non-positive person-years in 85+ aggregation")
  h85 <- sum(h * L) / sum(L)
  out[[length(out)]] <- data.frame(age_group = "85+", base_hazard = h85)

  ans <- do.call(rbind, out)
  ans$base_probability <- 1 - exp(-ans$base_hazard)
  rownames(ans) <- NULL
  ans
}

make_country_risks <- function(base, sex, HR, LCI, UCI, per = 100000) {
  rows <- vector("list", nrow(base))
  for (i in seq_len(nrow(base))) {
    h0 <- base$base_hazard[i]
    rows[[i]] <- data.frame(
      Age_group = base$age_group[i],
      Sex = sex,
      GOLD = names(HR),
      risk_per100k = (1 - exp(-HR * h0)) * per,
      lci_per100k = (1 - exp(-LCI * h0)) * per,
      uci_per100k = (1 - exp(-UCI * h0)) * per,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

format_country_table <- function(long) {
  long$value <- paste0(
    round(long$risk_per100k), " (",
    round(long$lci_per100k), "-",
    round(long$uci_per100k), ")"
  )
  d <- long[, c("Age_group", "Sex", "GOLD", "value")]
  d$GOLD <- factor(d$GOLD, levels = c("GOLD1", "GOLD2", "GOLD3", "GOLD4"))
  wide <- reshape(d, idvar = c("Age_group", "Sex"), timevar = "GOLD", direction = "wide")
  names(wide) <- sub("^value\\.", "", names(wide))
  names(wide)[names(wide) == "GOLD1"] <- "GOLD 1"
  names(wide)[names(wide) == "GOLD2"] <- "GOLD 2"
  names(wide)[names(wide) == "GOLD3"] <- "GOLD 3"
  names(wide)[names(wide) == "GOLD4"] <- "GOLD 4"
  wide <- wide[
    order(match(wide$Age_group, age_groups_keep), match(wide$Sex, c("Male", "Female"))),
    , drop = FALSE
  ]
  names(wide)[names(wide) == "Age_group"] <- "Age group"
  rownames(wide) <- NULL
  wide
}

calculate_german_period <- function(file, period_label) {
  male <- read_german_lifetable(file, "csv-12613-b03")
  female <- read_german_lifetable(file, "csv-12613-b04")

  base_m <- base_hazards_germany(male)
  base_f <- base_hazards_germany(female)

  long <- rbind(
    make_country_risks(base_m, "Male", HR, LCI, UCI),
    make_country_risks(base_f, "Female", HR, LCI, UCI)
  )
  long$Period <- period_label

  baseline <- rbind(
    data.frame(Age_group = base_m$age_group, Sex = "Male", baseline_per100k = base_m$base_probability * 100000),
    data.frame(Age_group = base_f$age_group, Sex = "Female", baseline_per100k = base_f$base_probability * 100000)
  )
  baseline$Period <- period_label
  baseline <- baseline[
    order(match(baseline$Age_group, age_groups_keep), match(baseline$Sex, c("Male", "Female"))),
    , drop = FALSE
  ]

  list(long = long, table = format_country_table(long), baseline = baseline)
}

## ============================================================================
## 9. TABLE 4 + e-TABLES 5-6: GERMANY
## ============================================================================

main_ger <- calculate_german_period(life_2022_24_file, "2022/24")
Table4 <- main_ger$table
cat("\nTABLE 4 - GERMANY 2022/24\n")
print(Table4, row.names = FALSE)
write.csv2(Table4, file.path(output_dir, "Table4_Germany_2022_24.csv"), row.names = FALSE)

# e-Table 5: baseline mortality of the normal-spirometry reference group.
eTable5 <- main_ger$baseline
eTable5$baseline_per100k <- round(eTable5$baseline_per100k)
names(eTable5)[names(eTable5) == "Age_group"] <- "Age group"
names(eTable5)[names(eTable5) == "baseline_per100k"] <- "Baseline 1-year mortality per 100,000"
cat("\ne-TABLE 5 - GERMAN BASELINE MORTALITY 2022/24\n")
print(eTable5, row.names = FALSE)
write.csv2(eTable5, file.path(output_dir, "eTable5_Germany_2022_24_baseline.csv"), row.names = FALSE)

sens_2023_25 <- calculate_german_period(life_2023_25_file, "2023/25")
sens_2017_19 <- calculate_german_period(life_2017_19_file, "2017/19")

eTable6 <- sens_2017_19$table
cat("\ne-TABLE 6 - GERMANY 2017/19 LIFE-TABLE SENSITIVITY\n")
print(eTable6, row.names = FALSE)
write.csv2(eTable6, file.path(output_dir, "eTable6_Germany_2017_19.csv"), row.names = FALSE)

compare_period <- function(main_long, sens_long, label) {
  key <- c("Age_group", "Sex", "GOLD")
  a <- main_long[, c(key, "risk_per100k")]
  b <- sens_long[, c(key, "risk_per100k")]
  names(a)[4] <- "main"
  names(b)[4] <- "sensitivity"
  x <- merge(a, b, by = key, all = FALSE)
  x$abs_difference <- abs(x$sensitivity - x$main)
  x$rel_difference_pct <- 100 * x$abs_difference / x$main
  data.frame(
    Period = label,
    max_absolute_difference_per100k = max(x$abs_difference, na.rm = TRUE),
    max_relative_difference_pct = max(x$rel_difference_pct, na.rm = TRUE)
  )
}

life_table_sensitivity_summary <- rbind(
  compare_period(main_ger$long, sens_2023_25$long, "2023/25 vs 2022/24"),
  compare_period(main_ger$long, sens_2017_19$long, "2017/19 vs 2022/24")
)
cat("\nGERMAN LIFE-TABLE SENSITIVITY SUMMARY\n")
print(life_table_sensitivity_summary, row.names = FALSE)
write.csv2(
  life_table_sensitivity_summary,
  file.path(output_dir, "eTable6_sensitivity_summary.csv"),
  row.names = FALSE
)

## ============================================================================
## 10. e-TABLE 7: AUSTRIA 2022-2024 WORKED EXAMPLE
## ============================================================================

if (!file.exists(life_austria_file)) stop("Austrian life-table file not found: ", life_austria_file)

aut_raw <- readODS::read_ods(
  path = life_austria_file,
  sheet = life_austria_sheet,
  col_names = FALSE,
  skip = 7
)
aut_raw <- as.data.frame(aut_raw)

to_num <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- trimws(as.character(x))
  x <- gsub("\u00a0", "", x, fixed = TRUE)
  x <- gsub(" ", "", x, fixed = TRUE)
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

age <- to_num(aut_raw[[1]])
aut_raw <- aut_raw[is.finite(age), , drop = FALSE]
age <- as.integer(to_num(aut_raw[[1]]))

# Statistics Austria layout after skip=7:
# col 1 age; cols 2-7 male (qx,lx,dx,Lx,Tx,ex); cols 8-13 female.
austria_male <- data.frame(
  age = age,
  qx = to_num(aut_raw[[2]]),
  Lx = to_num(aut_raw[[5]])
)
austria_female <- data.frame(
  age = age,
  qx = to_num(aut_raw[[8]]),
  Lx = to_num(aut_raw[[11]])
)

austria_male <- austria_male[is.finite(austria_male$age) & is.finite(austria_male$qx) & is.finite(austria_male$Lx), ]
austria_female <- austria_female[is.finite(austria_female$age) & is.finite(austria_female$qx) & is.finite(austria_female$Lx), ]

base_hazards_austria <- function(x, eps = 1e-12) {
  groups <- list(
    "45 - 49" = 45:49,
    "50 - 54" = 50:54,
    "55 - 59" = 55:59,
    "60 - 64" = 60:64,
    "65 - 69" = 65:69,
    "70 - 74" = 70:74,
    "75 - 79" = 75:79,
    "80 - 84" = 80:84
  )
  out <- vector("list", length(groups) + 1)

  k <- 0
  for (nm in names(groups)) {
    ages <- groups[[nm]]
    d <- x[x$age %in% ages, , drop = FALSE]
    if (nrow(d) != length(ages)) stop("Incomplete Austrian age group: ", nm)
    h <- mean(-log1p(-pmin(d$qx, 1 - eps)))
    k <- k + 1
    out[[k]] <- data.frame(age_group = nm, base_hazard = h)
  }

  d85 <- x[x$age >= 85 & is.finite(x$qx) & is.finite(x$Lx) & x$Lx > 0, , drop = FALSE]
  if (!nrow(d85)) stop("No Austrian ages >=85 found")
  h85 <- -log1p(-pmin(d85$qx, 1 - eps))
  h85 <- sum(h85 * d85$Lx) / sum(d85$Lx)
  out[[length(out)]] <- data.frame(age_group = "85+", base_hazard = h85)

  ans <- do.call(rbind, out)
  ans$base_probability <- 1 - exp(-ans$base_hazard)
  ans
}

base_aut_m <- base_hazards_austria(austria_male)
base_aut_f <- base_hazards_austria(austria_female)
aut_long <- rbind(
  make_country_risks(base_aut_m, "Male", HR, LCI, UCI),
  make_country_risks(base_aut_f, "Female", HR, LCI, UCI)
)
eTable7 <- format_country_table(aut_long)

cat("\ne-TABLE 7 - AUSTRIA 2022-2024 WORKED EXAMPLE\n")
print(eTable7, row.names = FALSE)
write.csv2(eTable7, file.path(output_dir, "eTable7_Austria_2022_2024.csv"), row.names = FALSE)

## ============================================================================
## 11. REPRODUCIBILITY INFORMATION
## ============================================================================

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
cat("\nAnalysis complete. Outputs written to:", normalizePath(output_dir), "\n")
