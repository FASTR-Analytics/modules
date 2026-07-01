library(readr)
library(dplyr)

raw <- read_csv(PROJECT_DATA_ICEH, show_col_types = FALSE) # 0-100
coverage <- raw %>%
  mutate(
    estimate = estimate / 100,
    standard_error = standard_error / 100
  ) # 0-1 (existing pass-through)

# Poorest->richest level ordering. Single-sourced from ICEH_STRAT_INFO
# (lib/types/iceh_strats.ts); the R script cannot import lib, so this literal is
# kept verified against it (PLAN_ICEH_DERIVED_INDICES_OPTION_5.md, Refinement #2).
ORDERED <- list(
  wealth_quintiles = c("Q1", "Q2", "Q3", "Q4", "Q5"),
  wealth_deciles = sprintf("D%02d", 1:10)
)

# ---- Inequality measures: HAND-ROLLED v1 (grouped approximations) ------------
# Per ordered strat: y = level estimates (0-1, poorest->richest); w = population
# shares (= 1/n for equal-sized wealth groups); r = fractional-rank midpoints
# (= (i-0.5)/n). Ratio/Difference are EXACT. CIX/SII are APPROXIMATE: the published
# ICEH profile computes them from per-person MICRODATA (logistic SII, covariance
# CIX) that FASTR never has -- it ingests pre-aggregated Retriever stratum
# estimates. Grouped approximations are the fidelity ceiling here.
#
# OPTIONAL: swap in ICEH's official package `ICEHmeasures` (deliberately NOT used).
# Spike 2026-06-18 verdict: do not depend on it for v1. Why:
#   - cixr() on grouped rows is BYTE-IDENTICAL to cix1() below -- no accuracy gain.
#   - siilogit() (grouped logistic SII) is NOT reliably closer to the published
#     numbers than sii_lin() -- sometimes worse (DPT3: lib 43.9 vs published 38.4
#     vs linear 41.0) -- because the dominant error is grouped-vs-microdata, which
#     no method crosses without microdata.
#   - cixr/siilogit are DESIGNED for microdata (ICEH's own examples use per-person
#     rows with survey weight + PSU + a continuous wealth score); feeding them
#     grouped rows is off-label, and cluster_var survey SEs are meaningless here.
#   - it pulls car + survey + lme4 + pbkrtest + doBy into the R container.
# The only reason to adopt it is provenance/credibility ("computed by ICEH's own
# package"), an accepted NON-goal for accuracy. If you ever want that, add
# `ICEHmeasures` to the R container and replace cix1/sii_lin with grouped-mode
# wrappers -- same call sites in ineq_for(), nothing else changes:
#
#   library(ICEHmeasures)
#   cix1 <- function(y, w, r) {
#     df <- tibble::tibble(rank = r, outcome = y, wt = w)
#     as.numeric(ICEHmeasures::cixr(df, rank, outcome, weight_var = wt)$cix[1])
#   }
#   sii_lin <- function(y, w, r) {
#     df <- tibble::tibble(rank = r, outcome = y, wt = w)
#     as.numeric(ICEHmeasures::siilogit(df, rank, outcome, weight_var = wt)$sii[1])
#   }
#   # leave cluster_var = NULL (no PSU in grouped data); siilogit fits glm(binomial)
#   # on a proportion response -- re-verify in the harness before trusting it.
cix1 <- function(y, w, r) {
  mu <- sum(w * y)
  if (mu == 0) NA_real_ else (2 / mu) * sum(w * y * r) - 1
}
sii_lin <- function(y, w, r) {
  if (length(y) < 2) return(NA_real_) # guards empty/degenerate groups (e.g. a wealth ladder absent)
  unname(coef(lm(y ~ r, weights = w))[["r"]])
}

ineq_for <- function(d, dim, lvls) {
  n_expected <- length(lvls)
  r <- (seq_len(n_expected) - 0.5) / n_expected # fractional-rank midpoints
  w <- rep(1 / n_expected, n_expected) # equal wealth-group shares
  d %>%
    filter(strat == dim, level %in% lvls) %>%
    group_by(iceh_indicator, year, source) %>%
    filter(n() == n_expected, !any(is.na(estimate))) %>% # complete ordered sets only
    arrange(match(level, lvls), .by_group = TRUE) %>% # poorest -> richest
    summarise(
      strat = dim,
      ratio = last(estimate) / first(estimate), # Q5 / Q1 (exact)
      difference = 100 * (last(estimate) - first(estimate)), # Q5 - Q1, pp (exact)
      cix = 100 * cix1(estimate, w, r), # concentration index (approx)
      sii = 100 * sii_lin(estimate, w, r), # slope index, linear (approx)
      .groups = "drop"
    )
}

inequality <- bind_rows(
  lapply(names(ORDERED), function(dim) ineq_for(coverage, dim, ORDERED[[dim]]))
)

write_csv(coverage, "M9_iceh_data.csv") # RO-1: coverage + ingested cci (unchanged)
write_csv(inequality, "M9_iceh_inequality.csv") # RO-2: indicator x strat x {ratio,diff,cix,sii}
