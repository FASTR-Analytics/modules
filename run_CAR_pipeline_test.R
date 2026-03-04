# Ad-hoc script: Run full pipeline for CAR (CAF) to test MICS Sub-national data
# Usage: source("run_CAR_pipeline_test.R")
# Pipeline: m001 → m002 → m004 (all-in-one) AND m005 → m006, then compare
# -------------------------------------------------------------------

cat("\n",
    "=================================================================\n",
    "  CAR (CAF) Pipeline Test — MICS Sub-national Verification\n",
    "=================================================================\n\n")

# ── Helper: read module, replace param lines via regex, eval in global env ──
run_module <- function(file, overrides) {
  code <- readLines(file)
  for (o in overrides) {
    code <- sub(o$pattern, o$replacement, code)
  }
  eval(parse(text = code), envir = globalenv())
}

# ── Common overrides ────────────────────────────────────────────────

# ISO3 and HMIS file
iso3_sub <- list(
  pattern     = '^COUNTRY_ISO3\\s*<-.*',
  replacement = 'COUNTRY_ISO3 <- "CAF"'
)
hmis_sub <- list(
  pattern     = '^PROJECT_DATA_HMIS\\s*<-.*',
  replacement = 'PROJECT_DATA_HMIS <- "hmis_RCA.csv"'
)

# M001-specific: no opd in CAR HMIS, add delivery consistency pair
dqa_ind_sub <- list(
  pattern     = '^DQA_INDICATORS\\s*<-.*',
  replacement = 'DQA_INDICATORS <- c("penta1", "anc1", "delivery")'
)
consistency_sub <- list(
  pattern     = '^CONSISTENCY_PAIRS_USED\\s*<-.*',
  replacement = 'CONSISTENCY_PAIRS_USED <- c("penta", "anc", "delivery")'
)

# M004/M005 count variable and analysis level
count_var_sub <- list(
  pattern     = '^SELECTED_COUNT_VARIABLE\\s*<-.*',
  replacement = 'SELECTED_COUNT_VARIABLE <- "count_final_both"'
)
analysis_level_sub <- list(
  pattern     = '^ANALYSIS_LEVEL\\s*<-.*',
  replacement = 'ANALYSIS_LEVEL <- "NATIONAL_PLUS_AA2"'
)

# CAR demographic rates (UN IGME / UNWPP estimates)
nmr_sub <- list(
  pattern     = '^P1_NMR\\s*<-.*',
  replacement = 'P1_NMR <- 0.033'
)
pnmr_sub <- list(
  pattern     = '^P2_PNMR\\s*<-.*',
  replacement = 'P2_PNMR <- 0.030'
)
imr_sub <- list(
  pattern     = '^INFANT_MORTALITY_RATE\\s*<-.*',
  replacement = 'INFANT_MORTALITY_RATE <- 0.092'
)
u5mr_sub <- list(
  pattern     = '^UNDER5_MORTALITY_RATE\\s*<-.*',
  replacement = 'UNDER5_MORTALITY_RATE <- 0.130'
)
stillbirth_sub <- list(
  pattern     = '^STILLBIRTH_RATE\\s*<-.*',
  replacement = 'STILLBIRTH_RATE <- 0.023'
)
preg_loss_sub <- list(
  pattern     = '^PREGNANCY_LOSS_RATE\\s*<-.*',
  replacement = 'PREGNANCY_LOSS_RATE <- 0.03'
)
twin_sub <- list(
  pattern     = '^TWIN_RATE\\s*<-.*',
  replacement = 'TWIN_RATE <- 0.015'
)

# M006 denominator selection — all "best" (auto-select)
denom_pregnancy_sub <- list(
  pattern     = '^DENOM_PREGNANCY\\s*<-.*',
  replacement = 'DENOM_PREGNANCY <- "best"'
)
denom_livebirth_sub <- list(
  pattern     = '^DENOM_LIVEBIRTH\\s*<-.*',
  replacement = 'DENOM_LIVEBIRTH <- "best"'
)
denom_dpt_sub <- list(
  pattern     = '^DENOM_DPT\\s*<-.*',
  replacement = 'DENOM_DPT <- "best"'
)
denom_measles1_sub <- list(
  pattern     = '^DENOM_MEASLES1\\s*<-.*',
  replacement = 'DENOM_MEASLES1 <- "best"'
)
denom_measles2_sub <- list(
  pattern     = '^DENOM_MEASLES2\\s*<-.*',
  replacement = 'DENOM_MEASLES2 <- "best"'
)
denom_vitamina_sub <- list(
  pattern     = '^DENOM_VITAMINA\\s*<-.*',
  replacement = 'DENOM_VITAMINA <- "best"'
)
denom_fullimm_sub <- list(
  pattern     = '^DENOM_FULLIMM\\s*<-.*',
  replacement = 'DENOM_FULLIMM <- "best"'
)

# ── Grouped override lists per module ───────────────────────────────

overrides_m001 <- list(iso3_sub, hmis_sub, dqa_ind_sub, consistency_sub)

overrides_m002 <- list(iso3_sub, hmis_sub)

overrides_m004 <- list(iso3_sub, count_var_sub, analysis_level_sub,
                       nmr_sub, pnmr_sub, imr_sub, u5mr_sub,
                       stillbirth_sub, preg_loss_sub, twin_sub)

overrides_m005 <- list(iso3_sub, count_var_sub, analysis_level_sub,
                       nmr_sub, pnmr_sub, imr_sub, u5mr_sub,
                       stillbirth_sub, preg_loss_sub, twin_sub)

overrides_m006 <- list(iso3_sub,
                       denom_pregnancy_sub, denom_livebirth_sub,
                       denom_dpt_sub, denom_measles1_sub,
                       denom_measles2_sub, denom_vitamina_sub,
                       denom_fullimm_sub)

# ═══════════════════════════════════════════════════════════════════
#  STEP 1: Run M001 (Data Quality Assessment)
# ═══════════════════════════════════════════════════════════════════
message("\n===== STEP 1: Running M001 (Data Quality Assessment) for CAF =====\n")
run_module("m001_module_data_quality_assessment.R", overrides_m001)
message("\n✓ M001 complete.\n")

# ═══════════════════════════════════════════════════════════════════
#  STEP 2: Run M002 (Data Quality Adjustments)
# ═══════════════════════════════════════════════════════════════════
message("\n===== STEP 2: Running M002 (Data Quality Adjustments) for CAF =====\n")
run_module("m002_module_data_quality_adjustments.R", overrides_m002)
message("\n✓ M002 complete.\n")

# ═══════════════════════════════════════════════════════════════════
#  STEP 3: Run M004 (All-in-one Coverage Estimates)
# ═══════════════════════════════════════════════════════════════════
message("\n===== STEP 3: Running M004 (All-in-One Coverage) for CAF =====\n")
run_module("m004_module_coverage_estimates.R", overrides_m004)
message("\n✓ M004 complete.\n")

# ═══════════════════════════════════════════════════════════════════
#  STEP 4: Back up M004 outputs (prevent confusion with M006)
# ═══════════════════════════════════════════════════════════════════
message("\n===== STEP 4: Backing up M004 outputs =====\n")

m4_files <- c(
  "M4_coverage_estimation.csv",
  "M4_coverage_estimation_admin_area_2.csv",
  "M4_coverage_estimation_admin_area_3.csv",
  "M4_selected_denominator_per_indicator.csv"
)

for (f in m4_files) {
  if (file.exists(f)) {
    backup <- sub("^M4_", "m004_backup_", f)
    file.copy(f, backup, overwrite = TRUE)
    message("  Backed up: ", f, " -> ", backup)
  }
}
message("\n✓ M004 backups complete.\n")

# ═══════════════════════════════════════════════════════════════════
#  STEP 5: Run M005 (Coverage Part 1 — Denominators)
# ═══════════════════════════════════════════════════════════════════
message("\n===== STEP 5: Running M005 (Coverage Part 1 — Denominators) for CAF =====\n")
run_module("m005_module_coverage_estimates_part1.R", overrides_m005)
message("\n✓ M005 complete.\n")

# ═══════════════════════════════════════════════════════════════════
#  STEP 6: Run M006 (Coverage Part 2 — Selection & Projection)
# ═══════════════════════════════════════════════════════════════════
message("\n===== STEP 6: Running M006 (Coverage Part 2 — Selection) for CAF =====\n")
run_module("m006_module_coverage_estimates_part2.R", overrides_m006)
message("\n✓ M006 complete.\n")

# ═══════════════════════════════════════════════════════════════════
#  STEP 7: Comparison — M004 vs M006 outputs
# ═══════════════════════════════════════════════════════════════════
message("\n===== STEP 7: Comparing M004 vs M006 outputs =====\n")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ── 7a. Load M4 backups and M6 outputs ──

# National level
m4_nat <- if (file.exists("m004_backup_coverage_estimation.csv")) {
  fread("m004_backup_coverage_estimation.csv")
} else {
  message("  WARNING: m004_backup_coverage_estimation.csv not found")
  NULL
}

m6_nat <- if (file.exists("M6_coverage_estimation_national.csv")) {
  fread("M6_coverage_estimation_national.csv")
} else {
  message("  WARNING: M6_coverage_estimation_national.csv not found")
  NULL
}

# Admin2 level
m4_aa2 <- if (file.exists("m004_backup_coverage_estimation_admin_area_2.csv")) {
  fread("m004_backup_coverage_estimation_admin_area_2.csv")
} else {
  message("  WARNING: m004_backup_coverage_estimation_admin_area_2.csv not found")
  NULL
}

m6_aa2 <- if (file.exists("M6_coverage_estimation_admin2.csv")) {
  fread("M6_coverage_estimation_admin2.csv")
} else {
  message("  WARNING: M6_coverage_estimation_admin2.csv not found")
  NULL
}

# ── 7b. Row count comparison ──
cat("\n--- Row Counts ---\n")
cat(sprintf("  M4 national:  %s rows\n", if (!is.null(m4_nat)) nrow(m4_nat) else "N/A"))
cat(sprintf("  M6 national:  %s rows\n", if (!is.null(m6_nat)) nrow(m6_nat) else "N/A"))
cat(sprintf("  M4 admin2:    %s rows\n", if (!is.null(m4_aa2)) nrow(m4_aa2) else "N/A"))
cat(sprintf("  M6 admin2:    %s rows\n", if (!is.null(m6_aa2)) nrow(m6_aa2) else "N/A"))

# ── 7c. Coverage value comparison (admin2 level) ──
if (!is.null(m4_aa2) && !is.null(m6_aa2)) {

  # Identify common join columns
  # M4 uses indicator_common_id, year, admin_area_2
  # M6 may use indicator or indicator_common_id — detect
  m4_ind_col <- intersect(c("indicator_common_id", "indicator"), names(m4_aa2))[1]
  m6_ind_col <- intersect(c("indicator_common_id", "indicator"), names(m6_aa2))[1]

  m4_year_col <- intersect(c("year"), names(m4_aa2))[1]
  m6_year_col <- intersect(c("year"), names(m6_aa2))[1]

  m4_aa2_col <- intersect(c("admin_area_2"), names(m4_aa2))[1]
  m6_aa2_col <- intersect(c("admin_area_2"), names(m6_aa2))[1]

  # Coverage column — M4/M6 use "coverage_cov"
  m4_cov_col <- intersect(c("coverage_cov", "coverage", "coverage_pct", "coverage_estimate"), names(m4_aa2))[1]
  m6_cov_col <- intersect(c("coverage_cov", "coverage", "coverage_pct", "coverage_estimate"), names(m6_aa2))[1]

  cat(sprintf("\n  M4 columns: %s\n", paste(names(m4_aa2), collapse = ", ")))
  cat(sprintf("  M6 columns: %s\n", paste(names(m6_aa2), collapse = ", ")))

  if (!is.na(m4_ind_col) && !is.na(m6_ind_col) &&
      !is.na(m4_year_col) && !is.na(m6_year_col) &&
      !is.na(m4_aa2_col) && !is.na(m6_aa2_col) &&
      !is.na(m4_cov_col) && !is.na(m6_cov_col)) {

    # Standardize column names for join
    m4_join <- m4_aa2 %>%
      select(indicator = !!sym(m4_ind_col),
             year      = !!sym(m4_year_col),
             admin_area_2 = !!sym(m4_aa2_col),
             coverage_m4  = !!sym(m4_cov_col)) %>%
      mutate(across(c(year), as.integer),
             coverage_m4 = as.numeric(coverage_m4))

    m6_join <- m6_aa2 %>%
      select(indicator = !!sym(m6_ind_col),
             year      = !!sym(m6_year_col),
             admin_area_2 = !!sym(m6_aa2_col),
             coverage_m6  = !!sym(m6_cov_col)) %>%
      mutate(across(c(year), as.integer),
             coverage_m6 = as.numeric(coverage_m6))

    merged <- inner_join(m4_join, m6_join,
                         by = c("indicator", "year", "admin_area_2"))

    cat(sprintf("\n--- Coverage Comparison (admin2) ---\n"))
    cat(sprintf("  Matched rows: %d\n", nrow(merged)))

    if (nrow(merged) > 0) {
      merged <- merged %>%
        mutate(abs_diff = abs(coverage_m4 - coverage_m6))

      cat(sprintf("  Max absolute difference: %.4f\n", max(merged$abs_diff, na.rm = TRUE)))
      cat(sprintf("  Mean absolute difference: %.4f\n", mean(merged$abs_diff, na.rm = TRUE)))

      # Show rows with discrepancies > 0.01
      discrepancies <- merged %>% filter(abs_diff > 0.01)
      if (nrow(discrepancies) > 0) {
        cat(sprintf("\n  WARNING: %d rows with abs diff > 0.01:\n", nrow(discrepancies)))
        print(head(as.data.frame(discrepancies), 20))
      } else {
        cat("  PASS: All matched rows agree within 0.01\n")
      }
    }
  } else {
    cat("\n  Could not align columns for comparison.\n")
    cat(sprintf("  M4 detected: ind=%s, year=%s, aa2=%s, cov=%s\n",
                m4_ind_col, m4_year_col, m4_aa2_col, m4_cov_col))
    cat(sprintf("  M6 detected: ind=%s, year=%s, aa2=%s, cov=%s\n",
                m6_ind_col, m6_year_col, m6_aa2_col, m6_cov_col))
  }
} else {
  cat("\n  Skipping coverage comparison — missing admin2 output(s).\n")
}

# ═══════════════════════════════════════════════════════════════════
#  STEP 8: MICS Sub-national Source Verification
# ═══════════════════════════════════════════════════════════════════
message("\n===== STEP 8: MICS Sub-national Source Verification =====\n")

check_mics_source <- function(df, label) {
  if (is.null(df)) {
    cat(sprintf("  %s: SKIPPED (file not found)\n", label))
    return(invisible(NULL))
  }

  # Look for source column — may be "survey_source", "survey_raw_source", "source", etc.
  source_cols <- intersect(
    c("survey_source", "survey_raw_source", "source", "denom_source", "survey_type",
      "survey_raw_source_detail", "survey_source_detail"),
    names(df)
  )

  # Look for admin_area_2 subnational rows
  aa2_col <- intersect(c("admin_area_2"), names(df))

  if (length(source_cols) > 0) {
    src_col <- source_cols[1]
    all_sources <- unique(df[[src_col]])
    cat(sprintf("  %s — source column '%s': %s\n",
                label, src_col, paste(all_sources, collapse = ", ")))

    has_mics <- any(grepl("mics", tolower(as.character(all_sources))))
    if (has_mics) {
      cat(sprintf("  %s — PASS: 'mics' source detected\n", label))
    } else {
      cat(sprintf("  %s — INFO: 'mics' not found in source column\n", label))
    }

    # Check subnational MICS rows
    if (length(aa2_col) > 0) {
      subnational_rows <- df[grepl("^RS ", df[[aa2_col]]), ]
      if (nrow(subnational_rows) > 0) {
        sub_sources <- unique(subnational_rows[[src_col]])
        cat(sprintf("  %s — RS 1-7 rows: %d, sources: %s\n",
                    label, nrow(subnational_rows),
                    paste(sub_sources, collapse = ", ")))
      } else {
        cat(sprintf("  %s — No RS 1-7 rows found\n", label))
      }
    }
  } else {
    cat(sprintf("  %s — No source column found (columns: %s)\n",
                label, paste(names(df), collapse = ", ")))
  }

  # Check for non-NA survey values at subnational level
  survey_cols <- intersect(
    c("survey_value", "survey_coverage", "survey_estimate",
      "coverage_avgsurveyprojection", "coverage_cov"),
    names(df)
  )
  if (length(survey_cols) > 0 && length(aa2_col) > 0) {
    sv_col <- survey_cols[1]
    sub_with_survey <- df[grepl("^RS ", df[[aa2_col]]) &
                            !is.na(df[[sv_col]]), ]
    cat(sprintf("  %s — RS 1-7 rows with non-NA '%s': %d\n",
                label, sv_col, nrow(sub_with_survey)))
    if (nrow(sub_with_survey) > 0) {
      cat(sprintf("  %s — PASS: MICS subnational survey data present\n", label))
    } else {
      cat(sprintf("  %s — WARN: No subnational survey values found\n", label))
    }
  }
}

# Check all available outputs
check_mics_source(m4_aa2, "M4 admin2")
check_mics_source(m6_aa2, "M6 admin2")

# Also check M5 combined results if available
if (file.exists("M5_combined_results_admin2.csv")) {
  m5_aa2 <- fread("M5_combined_results_admin2.csv")
  check_mics_source(m5_aa2, "M5 admin2")
}

# ═══════════════════════════════════════════════════════════════════
#  Summary
# ═══════════════════════════════════════════════════════════════════
cat("\n",
    "=================================================================\n",
    "  CAR (CAF) Pipeline Test Complete\n",
    "=================================================================\n",
    "  Outputs produced:\n",
    "    M1: M1_output_*.csv\n",
    "    M2: M2_adjusted_data*.csv\n",
    "    M4: M4_coverage_estimation*.csv (backups: m004_backup_*.csv)\n",
    "    M5: M5_*.csv\n",
    "    M6: M6_coverage_estimation_*.csv\n",
    "=================================================================\n\n")
