COUNTRY_ISO3 <- "ZMB"

SELECTEDCOUNT <- "count_final_outliers"  # count_final_none | count_final_outliers ONLY — completeness-imputed counts (count_final_completeness/_both) are rejected: this module models reporting gaps itself and must see them in the data
FOURIER_K <- 2  # seasonal wave-pairs: 2 = annual + semi-annual; 3-4 for campaign-spike indicators (e.g. vitamin A)
USE_POSTERIOR_DRAWS <- FALSE  # TRUE = properly calibrated subnational CIs, ~10x slower and memory-heavy
ZEROS_REAL <- "all"  # all | detect | auto | strict | none — zero-handling policy, see header
EXCLUDE_PRELAUNCH <- TRUE  # drop grid months before a facility's first-ever report (new facilities / late-rollout indicators)
POSTCLOSURE_GRACE <- 0  # e.g. 12 drops grid months >12m after a facility's last-ever report (presumed closure); 0 = off

PROJECT_DATA_HMIS <- "hmis_ZMB.csv"  # injected per country by the platform; local-testing default
#-------------------------------------------------------------------------------------------------------------
# FASTR PROJECT
# Module: BAYESIAN DISRUPTION DETECTION (LI MODEL)
# Last edit: 2026 Jul 28
#
# Bayesian model author: Mustapha Wasseja
# Module integration:    CB
#
# Two-part Bayesian "Latent Intensity" model fit with INLA that structurally
# separates reporting from service delivery so detected disruptions reflect
# REAL service-rate changes, not artefacts of facilities going dark.
#
#   Part 1 (Reporting): did_report(f,t) ~ Bernoulli(p(f,t))
#     Fitted on the full facility x month grid (1 = reported, 0 = missing).
#
#   Part 2 (Service):   value(f,t) | reported ~ NegBin(mu(f,t), phi)
#     Fitted on reporters only. mu evolves through linear trend, Fourier
#     seasonality, and a facility random effect.
#
#   Detection baseline: sum of mu over the SAME facilities that reported each
#     month (apples-to-apples vs observed). Deviations outside the 95% CI of
#     this quantity are real service-rate changes among active reporters.
#
# Zero handling (ZEROS_REAL): in some HMIS instances an explicit 0 means
# "facility reported, no activity" (a real zero); in others the system
# auto-fills 0 for non-submission (a fake zero that should be missing).
#   "all"  (default) - every 0 is a real report of zero activity. This is the
#                      canonical behaviour and matches Mustapha's scripts as
#                      shipped (his "auto" mode never reclassifies in practice).
#   "auto"           - a 0 is real only if the facility submitted at least one
#                      OTHER indicator that same month; otherwise it is
#                      reclassified as missing (non-reporting). This is the
#                      corrected version of the cross-indicator rule described
#                      in the header of FASTR_BayesianLI_Unified.R.
#   "strict"         - like "auto", but the other indicator must have a value
#                      > 0. For HMIS instances that auto-fill 0s across ALL
#                      indicators at once, where under "auto" the fake zeros
#                      would vouch for each other.
#   "none"           - every 0 is treated as missing.
#
# New-facility handling (EXCLUDE_PRELAUNCH, default TRUE): grid months before
# a facility's first-ever report are dropped, so facilities that opened partway
# through the period don't count their pre-launch months as non-reporting
# (which would deflate the Part 1 reporting probabilities and skew the gap
# decomposition). Disruption flags are unaffected either way.
#
# Closed-facility handling (POSTCLOSURE_GRACE, default 0 = off): the mirror
# problem — a facility that stopped reporting for good keeps counting as
# non-reporting until the end of the period. With grace N > 0, grid months
# more than N months after a facility's LAST-ever report are dropped as
# presumed closure; the first N silent months still count as observable
# non-reporting, so genuine reporting collapse remains visible to Part 1.
# Off by default because trailing silence is ambiguous (closure vs ongoing
# reporting failure). Diagnostic: test_facility_lifespan.R — as of Jul 2026
# only Liberia has a material closure tail (10.7% of facilities dark >=12
# months; 6.4% of the modeled grid).
#
# Gap decomposition (per admin area x period), from the unified script:
#   expected_full  = sum of model mu over ALL facilities (reporters or not)
#   expected_as    = sum of p * mu (reporting-probability-weighted expectation)
#   gap_reporting  = expected_full - expected_as  (shortfall from non-reporting)
#   gap_service    = expected_as - observed       (shortfall from service change)
# The flagging baseline is unchanged: `expected` = sum of mu over the SAME
# facilities that reported that month.
#
# The model fits ONCE at facility x month level per indicator. Aggregation to
# admin_area_1/2/3/4 is a group-by-and-sum on the per-facility posteriors —
# no separate model is run per admin level. The script auto-detects the
# deepest available admin_area level in the raw HMIS and emits rollups up to
# it. Levels not present in the country's data get a schema-only empty CSV
# (matching the FASTR convention).
#
# ------------------------------------- KEY OUTPUTS ----------------------------------------------------------
# FILE: M11_disruptions_analysis_admin_area_1.csv  # National
# FILE: M11_disruptions_analysis_admin_area_2.csv  # admin_area_2
# FILE: M11_disruptions_analysis_admin_area_3.csv  # admin_area_3
# FILE: M11_disruptions_analysis_admin_area_4.csv  # admin_area_4 (or empty if country has no admin_area_4)
#-------------------------------------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(INLA)
})

select <- dplyr::select

TRIM_FROM_PERIOD <- ""  # e.g. 202001 drops data before Jan 2020; "" = use all history
INLA_THREADS <- "1:1"  # e.g. "1:1" forces single-threaded fits (several-fold lower peak memory); "" = INLA default


#-------------------------------------------------------------------------------------------------------------
# STEP 1: LOAD DATA
#-------------------------------------------------------------------------------------------------------------
stopifnot(
  "SELECTEDCOUNT must be count_final_none or count_final_outliers: completeness-imputed counts fill in the very reporting gaps this module exists to model" =
    !SELECTEDCOUNT %in% c("count_final_completeness", "count_final_both")
)

print("Loading data...")
raw_data <- fread(PROJECT_DATA_HMIS)
adj_data <- fread("M2_adjusted_data.csv")

# Auto-detect available admin_area_X columns in raw HMIS
admin_cols <- grep("^admin_area_[0-9]+$", names(raw_data), value = TRUE)
admin_cols <- admin_cols[order(as.integer(sub("admin_area_", "", admin_cols)))]
deepest_level <- as.integer(sub("admin_area_", "", tail(admin_cols, 1)))
print(paste0("[Geo] available admin levels: ",
             paste(admin_cols, collapse = ", "),
             " (deepest = admin_area_", deepest_level, ")"))

# M2_adjusted_data already carries some admin levels (m002 auto-detects them
# too). Pull from raw HMIS only the levels that aren't already on adj_data,
# to avoid column-name collisions on merge.
raw_only_cols <- setdiff(admin_cols, names(adj_data))
admin_lookup <- unique(raw_data[, c("facility_id", raw_only_cols), with = FALSE])
rm(raw_data); gc()

# Trim window (YYYYMM int parameter; blank/NA = use all data)
trim_period <- suppressWarnings(as.integer(TRIM_FROM_PERIOD))
if (!is.na(trim_period) && trim_period > 0L) {
  before_n <- nrow(adj_data)
  adj_data <- adj_data[period_id >= trim_period]
  print(paste0("Trimmed to period_id >= ", trim_period,
               ": ", before_n, " -> ", nrow(adj_data), " rows"))
}

n_zeros <- sum(adj_data[[SELECTEDCOUNT]] == 0, na.rm = TRUE)
if (n_zeros > 0L) {
  print(paste0("[Info] ", n_zeros, " explicit zeros in ", SELECTEDCOUNT,
               " — zero policy: ZEROS_REAL = '", ZEROS_REAL, "'"))
}

# Join admin levels not already on adj_data (admin_area_1 + admin_area_4+)
adj_data <- merge(adj_data, admin_lookup, by = "facility_id", all.x = TRUE)
adj_data[, value := as.numeric(get(SELECTEDCOUNT))]
adj_data <- adj_data[!is.na(value)]
adj_data[, value := pmax(0L, as.integer(round(value)))]
adj_data[, year  := period_id %/% 100L]
adj_data[, month := period_id %%  100L]

# Keep only the columns the model uses — adj_data stays resident for the whole
# run, and M2 files carry several unused count variants (big-country RAM win)
adj_data <- adj_data[, unique(c("facility_id", "indicator_common_id", "period_id",
                                intersect(admin_cols, names(adj_data)),
                                "value", "year", "month")), with = FALSE]
gc()

# ZEROS_REAL = "detect": pick the policy from the data itself, using the same
# heuristic as test_zero_signature.R — zeros rare -> "all" (convention is
# blank-on-non-reporting); wholesale zero-fill signature (>10% of all-zero
# facility-months span 5+ indicators) -> "strict"; otherwise -> "auto".
stopifnot(ZEROS_REAL %in% c("all", "auto", "strict", "none", "detect"))
if (ZEROS_REAL == "detect") {
  zero_share <- mean(adj_data$value == 0)
  if (zero_share < 0.005) {
    ZEROS_REAL <- "all"
  } else {
    fm_det <- adj_data[, .(n_ind = uniqueN(indicator_common_id),
                           n_pos = uniqueN(indicator_common_id[value > 0])),
                       by = .(facility_id, period_id)]
    az <- fm_det[n_pos == 0]
    wholesale <- if (nrow(az) > 0L) az[n_ind >= 5L, .N] / nrow(az) else 0
    ZEROS_REAL <- if (wholesale > 0.10) "strict" else "auto"
    rm(fm_det, az)
  }
  print(paste0("[Zeros] ZEROS_REAL='detect' resolved to '", ZEROS_REAL,
               "' (zero share ", sprintf("%.1f%%", 100 * zero_share), ")"))
}

# For ZEROS_REAL = "auto"/"strict": per facility-month, how many indicators
# were submitted (any non-missing value), and how many had a value > 0?
# A tested 0 contributes 1 to its own n_ind_submitted (so "no OTHER indicator
# submitted" means n_ind_submitted == 1) but nothing to n_ind_pos.
if (ZEROS_REAL %in% c("auto", "strict")) {
  submit_counts <- adj_data[, .(n_ind_submitted = uniqueN(indicator_common_id),
                                n_ind_pos = uniqueN(indicator_common_id[value > 0])),
                            by = .(facility_id, period_id)]
}

#-------------------------------------------------------------------------------------------------------------
# STEP 2: HELPERS
#-------------------------------------------------------------------------------------------------------------
add_fourier <- function(df, K) {
  w <- 2 * pi * df$month / 12
  for (k in seq_len(K)) {
    df[[paste0("sin", k)]] <- sin(k * w)
    df[[paste0("cos", k)]] <- cos(k * w)
  }
  df
}

fourier_terms <- paste(c(paste0("sin", seq_len(FOURIER_K)),
                         paste0("cos", seq_len(FOURIER_K))), collapse = " + ")

# INLA threading: "" = INLA's default. Set e.g. "1:1" to force single-threaded
# fits — several-fold lower peak memory at some speed cost. Use on machines
# where whole-network indicators (30k+ facilities) hit the memory ceiling.
inla_nthreads <- if (nzchar(as.character(INLA_THREADS))) {
  as.character(INLA_THREADS)
} else {
  # getOption probes the INLA binary; fall back to NULL where it is absent
  # (e.g. mocked test runs) — inla() is never reached in that case anyway
  tryCatch(INLA::inla.getOption("num.threads"), error = function(e) NULL)
}

#-------------------------------------------------------------------------------------------------------------
# STEP 3: PER-INDICATOR LI MODEL
#-------------------------------------------------------------------------------------------------------------
out_admin1 <- list()
out_admin2 <- list()
out_admin3 <- list()
out_admin4 <- list()

write_admin <- function(rows_list, file, key_cols, quiet = FALSE) {
  schema_cols <- c(key_cols, "indicator_common_id", "period_id",
                   "observed", "expected", "ppi_lwr", "ppi_upr",
                   "expected_full", "expected_as",
                   "gap_reporting", "gap_service",
                   "flag_deficit", "flag_surplus", "mean_p", "n_facilities")

  if (length(rows_list) == 0L ||
      all(vapply(rows_list, nrow, integer(1)) == 0L)) {
    if (!quiet) print(paste0("  [empty] writing schema-only file: ", file))
    empty <- as.data.frame(setNames(replicate(length(schema_cols),
                                              integer(0), simplify = FALSE),
                                    schema_cols))
    write.csv(empty, file, row.names = FALSE)
    return(invisible(NULL))
  }

  out <- bind_rows(rows_list) %>%
    select(all_of(schema_cols))
  write.csv(out, file, row.names = FALSE)
  if (!quiet) print(paste0("  Saved ", nrow(out), " rows to ", file))
}

flush_outputs <- function(quiet = TRUE) {
  write_admin(out_admin1, "M11_disruptions_analysis_admin_area_1.csv",
              "admin_area_1", quiet = quiet)
  write_admin(out_admin2, "M11_disruptions_analysis_admin_area_2.csv",
              "admin_area_2", quiet = quiet)
  write_admin(out_admin3, "M11_disruptions_analysis_admin_area_3.csv",
              c("admin_area_2", "admin_area_3"), quiet = quiet)
  write_admin(out_admin4, "M11_disruptions_analysis_admin_area_4.csv",
              c("admin_area_2", "admin_area_3", "admin_area_4"), quiet = quiet)
}

indicators <- sort(unique(adj_data$indicator_common_id))
print(paste("Modeling", length(indicators), "indicators..."))

for (ind_name in indicators) {
  print(paste0("======================================================"))
  print(paste("LI model:", ind_name))
  print(paste0("======================================================"))

  fac_data <- adj_data[indicator_common_id == ind_name]
  if (nrow(fac_data) == 0L) { print("  [skip] no data after filter"); next }

  # Apply the zero policy for this indicator
  if (ZEROS_REAL != "all") {
    n_zero <- sum(fac_data$value == 0L)
    if (ZEROS_REAL %in% c("auto", "strict")) {
      fac_data <- merge(fac_data, submit_counts,
                        by = c("facility_id", "period_id"), all.x = TRUE)
      fake <- if (ZEROS_REAL == "auto") {
        fac_data$value == 0L & fac_data$n_ind_submitted == 1L
      } else {  # strict: needs another indicator with a value > 0
        fac_data$value == 0L & fac_data$n_ind_pos == 0L
      }
      print(paste0("  [Zeros] ", n_zero, " zeros | ", n_zero - sum(fake),
                   " real (", ZEROS_REAL, " rule) | ",
                   sum(fake), " fake -> reclassified as missing"))
      fac_data <- fac_data[!fake]
      fac_data[, c("n_ind_submitted", "n_ind_pos") := NULL]
    } else {  # "none"
      print(paste0("  [Zeros] ", n_zero,
                   " zeros treated as missing (ZEROS_REAL='none')"))
      fac_data <- fac_data[value != 0L]
    }
    if (nrow(fac_data) == 0L) { print("  [skip] no data after zero handling"); next }
  }

  # One bad indicator (degenerate data, INLA non-convergence) must not kill a
  # multi-hour multi-indicator run: fail it, log it, move to the next.
  ok <- tryCatch({

  all_facilities <- unique(fac_data[, c("facility_id", admin_cols), with = FALSE])
  all_periods    <- data.table(period_id = sort(unique(adj_data$period_id)))

  # Build full facility x period grid
  full_grid <- as_tibble(CJ(facility_id = all_facilities$facility_id,
                            period_id   = all_periods$period_id, sorted = FALSE)) %>%
    mutate(year  = period_id %/% 100L,
           month = period_id %%  100L) %>%
    left_join(all_facilities, by = "facility_id") %>%
    left_join(fac_data %>% select(facility_id, period_id, value),
              by = c("facility_id", "period_id")) %>%
    mutate(did_report = as.integer(!is.na(value)))

  # Drop pre-launch months: grid cells before a facility's first-ever report
  # would otherwise count as non-reporting and deflate its Part 1 estimate
  if (isTRUE(as.logical(EXCLUDE_PRELAUNCH))) {
    first_report <- fac_data[, .(first_period = min(period_id)), by = facility_id]
    n_before <- nrow(full_grid)
    full_grid <- full_grid %>%
      left_join(as_tibble(first_report), by = "facility_id") %>%
      dplyr::filter(period_id >= first_period) %>%
      select(-first_period)
    if (n_before > nrow(full_grid)) {
      print(paste0("  [Grid] dropped ", n_before - nrow(full_grid),
                   " pre-launch facility-periods (EXCLUDE_PRELAUNCH)"))
    }
  }

  # Optionally drop long trailing silence: grid months more than
  # POSTCLOSURE_GRACE months after a facility's last-ever report are treated
  # as presumed closure rather than ongoing non-reporting. 0 = off.
  pc_grace <- suppressWarnings(as.integer(POSTCLOSURE_GRACE))
  if (!is.na(pc_grace) && pc_grace > 0L) {
    last_report <- fac_data[, .(last_mi = max(period_id %/% 100L * 12L +
                                                period_id %% 100L)),
                            by = facility_id]
    n_before <- nrow(full_grid)
    full_grid <- full_grid %>%
      left_join(as_tibble(last_report), by = "facility_id") %>%
      dplyr::filter(year * 12L + month <= last_mi + pc_grace) %>%
      select(-last_mi)
    if (n_before > nrow(full_grid)) {
      print(paste0("  [Grid] dropped ", n_before - nrow(full_grid),
                   " post-closure facility-periods (POSTCLOSURE_GRACE=",
                   pc_grace, ")"))
    }
  }

  # Time/seasonality covariates (zt scaled on the retained grid)
  full_grid <- full_grid %>%
    mutate(time_index = (year - min(year)) * 12L + month) %>%
    add_fourier(K = FOURIER_K) %>%
    mutate(zt = as.numeric(scale(time_index)))

  fac_levels <- sort(unique(full_grid$facility_id))
  full_grid$fac_idx   <- match(full_grid$facility_id, fac_levels)
  full_grid$fac_idx_p <- full_grid$fac_idx
  n_total <- nrow(full_grid)
  n_rep   <- sum(full_grid$did_report)
  print(paste0("  Grid: ", n_total, " facility-periods | ",
               n_rep, " reported (",
               sprintf("%.1f%%", n_rep / n_total * 100), ") | ",
               length(fac_levels), " facilities"))

  # ── Part 1: Reporting (Bernoulli) ──
  print("  [Part 1] Fitting Bernoulli reporting model...")
  f_report <- as.formula(paste(
    "did_report ~ 1 + zt +", fourier_terms,
    "+ f(fac_idx_p, model = 'iid',",
    "    hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))"
  ))
  fit_p <- inla(f_report,
                family             = "binomial",
                Ntrials            = rep(1, n_total),
                data               = as.data.frame(full_grid),
                control.predictor  = list(compute = TRUE, link = 1),
                control.inla       = list(strategy = "adaptive"),
                num.threads        = inla_nthreads,
                verbose            = FALSE)
  full_grid$p_fitted <- fit_p$summary.fitted.values$`0.5quant`

  # ── Part 2: Service (NegBin on reporters) ──
  print("  [Part 2] Fitting NegBin service model on reporters...")
  reported_data <- full_grid %>% dplyr::filter(did_report == 1)
  reported_data$fac_idx_mu <- reported_data$fac_idx
  f_service <- as.formula(paste(
    "value ~ 1 + zt +", fourier_terms,
    "+ f(fac_idx_mu, model = 'iid',",
    "    hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))"
  ))
  fit_mu <- inla(f_service,
                 family            = "nbinomial",
                 data              = as.data.frame(reported_data),
                 control.predictor = list(compute = TRUE, link = 1),
                 control.inla      = list(strategy = "adaptive"),
                 num.threads       = inla_nthreads,
                 verbose           = FALSE)
  reported_data$mu_fitted <- fit_mu$summary.fitted.values$`0.5quant`
  reported_data$mu_lwr    <- fit_mu$summary.fitted.values$`0.025quant`
  reported_data$mu_upr    <- fit_mu$summary.fitted.values$`0.975quant`

  # Model-expected mu for EVERY grid cell, including non-reporting months —
  # feeds the gap decomposition (expected_full / expected_as). Rebuilt from
  # fixed effects + facility random effect, as in Mustapha's unified script.
  fx <- setNames(fit_mu$summary.fixed[, "mean"], rownames(fit_mu$summary.fixed))
  lin_pred <- fx[["(Intercept)"]] + fx[["zt"]] * full_grid$zt
  for (k in seq_len(FOURIER_K)) {
    lin_pred <- lin_pred +
      fx[[paste0("sin", k)]] * full_grid[[paste0("sin", k)]] +
      fx[[paste0("cos", k)]] * full_grid[[paste0("cos", k)]]
  }
  re_tab <- fit_mu$summary.random$fac_idx_mu
  fac_re <- re_tab$mean[match(full_grid$fac_idx, re_tab$ID)]
  fac_re[is.na(fac_re)] <- 0
  full_grid <- full_grid %>%
    mutate(mu_predicted = exp(lin_pred + fac_re)) %>%
    left_join(reported_data %>% select(facility_id, period_id, mu_fitted),
              by = c("facility_id", "period_id")) %>%
    mutate(mu_final = dplyr::coalesce(mu_fitted, mu_predicted))

  # Optional joint-posterior draws for properly calibrated subnational CI
  draws_mat <- NULL
  if (isTRUE(as.logical(USE_POSTERIOR_DRAWS))) {
    est_gb <- nrow(reported_data) * 500 * 8 / 1e9
    if (est_gb > 4) {
      print(sprintf("  [WARN] posterior-draws matrix will be ~%.1f GB for this indicator — consider USE_POSTERIOR_DRAWS=FALSE on large countries", est_gb))
    }
    print("  [Part 2b] Drawing 500 joint posterior samples...")
    samples <- inla.posterior.sample(500, fit_mu, intern = FALSE)
    pred_idx <- grep("^Predictor:", rownames(samples[[1]]$latent))
    draws_mat <- vapply(samples,
                        function(s) exp(s$latent[pred_idx, 1]),
                        numeric(length(pred_idx)))
    rm(samples)
  }

  #-----------------------------------------------------------------------------
  # STEP 4: ROLL UP TO ADMIN LEVELS
  #-----------------------------------------------------------------------------
  rollup <- function(group_cols) {
    # Full-grid quantities (incl. non-reporters): reporting completeness plus
    # the two counterfactual expectations behind the gap decomposition
    p_part <- full_grid %>%
      group_by(across(all_of(group_cols)), period_id) %>%
      summarise(mean_p        = mean(p_fitted),
                expected_full = sum(mu_final),
                expected_as   = sum(p_fitted * mu_final),
                .groups = "drop")

    if (is.null(draws_mat)) {
      # Marginal-sum CI: cheap, slightly conservative at fine admin levels
      s_part <- reported_data %>%
        group_by(across(all_of(group_cols)), period_id) %>%
        summarise(observed     = sum(value),
                  expected     = sum(mu_fitted),
                  ppi_lwr      = sum(mu_lwr),
                  ppi_upr      = sum(mu_upr),
                  n_facilities = dplyr::n_distinct(facility_id),
                  .groups = "drop")
    } else {
      # Joint-draw CI: aggregate within group across draws, then take quantiles
      rd <- reported_data
      rd$row_idx <- seq_len(nrow(rd))
      s_part <- rd %>%
        group_by(across(all_of(group_cols)), period_id) %>%
        summarise(observed     = sum(value),
                  expected     = sum(mu_fitted),
                  n_facilities = dplyr::n_distinct(facility_id),
                  idx_list     = list(row_idx),
                  .groups = "drop")
      s_part$ppi_lwr <- vapply(s_part$idx_list, function(idx) {
        as.numeric(quantile(colSums(draws_mat[idx, , drop = FALSE]), 0.025))
      }, numeric(1))
      s_part$ppi_upr <- vapply(s_part$idx_list, function(idx) {
        as.numeric(quantile(colSums(draws_mat[idx, , drop = FALSE]), 0.975))
      }, numeric(1))
      s_part$idx_list <- NULL
    }

    s_part %>%
      left_join(p_part, by = c(group_cols, "period_id")) %>%
      mutate(gap_reporting       = expected_full - expected_as,
             gap_service         = expected_as - observed,
             flag_deficit        = as.integer(observed < ppi_lwr),
             flag_surplus        = as.integer(observed > ppi_upr),
             indicator_common_id = ind_name)
  }

  out_admin1[[ind_name]] <- rollup("admin_area_1")
  out_admin2[[ind_name]] <- rollup("admin_area_2")
  out_admin3[[ind_name]] <- rollup(c("admin_area_2", "admin_area_3"))
  if ("admin_area_4" %in% admin_cols) {
    out_admin4[[ind_name]] <- rollup(c("admin_area_2", "admin_area_3", "admin_area_4"))
  }

  rm(full_grid, reported_data, fit_p, fit_mu, draws_mat); gc()

  # Crash-safety: persist accumulated results after every indicator, so a
  # mid-run failure (e.g. OOM on a later indicator) loses nothing already done
  flush_outputs(quiet = TRUE)
  TRUE
  }, error = function(e) {
    print(paste0("  [ERROR] ", ind_name, " failed: ", conditionMessage(e),
                 " — skipping this indicator"))
    gc()
    FALSE
  })
  if (!ok) next
}

#-------------------------------------------------------------------------------------------------------------
# STEP 5: WRITE OUTPUTS (final flush — files are also written incrementally
# after every indicator, so a mid-run crash loses at most the indicator in
# progress)
#-------------------------------------------------------------------------------------------------------------
print("Saving outputs...")
flush_outputs(quiet = FALSE)
print("=== ANALYSIS COMPLETE ===")
