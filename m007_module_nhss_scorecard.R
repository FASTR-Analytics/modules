COUNTRY_ISO3 <- "NGA"
SELECTED_COUNT_VARIABLE <- "count_final_none"



#-------------------------------------------------------------------------------------------------------------
# CB - R code FASTR PROJECT
# Module: NATIONAL HEALTH SECTOR SCORECARD (NHSS)
# Last edit: 2026 Feb 18
#
# Calculates NHSS scorecard indicators at multiple geographic levels.
# All indicator data comes from the M2 adjusted data pipeline.
# Population is the only separate asset file (LGA-level denominator from DHIS2).
#
# OUTPUTS (one per geo level):
#   M7_output_scorecard_{level}.csv         — Scorecard indicators (long format, areas + national)
#   M7_output_scorecard_summary_{level}.csv — Raw quarterly counts by area (debug)
#
# Load Required Libraries -----------------------------------------------------------------------------------
library(dplyr)
library(readr)
library(tidyr)

# Load Data ------------------------------------------------------------------------------------------------
ADJUSTED_DATA_FILE <- "M2_adjusted_data.csv"
POPULATION_FILE <- "total_population.csv"
adjusted_data <- read_csv(ADJUSTED_DATA_FILE, show_col_types = FALSE)
population <- read_csv(POPULATION_FILE, show_col_types = FALSE)

# Scorecard Parameters
SCORECARD_YEAR <- 2025
SCORECARD_QUARTER <- 4



# Population assumptions (quarterly rates)
PREGNANT_WOMEN_PCT <- 0.05 * 0.25
BIRTHS_PCT <- 0.04 * 0.25
WOMEN_15_49_PCT <- 0.22 * 0.25
CHILDREN_U5_PCT <- 0.12 * 0.25
INFANTS_0_6M_PCT <- 0.02

# Helper: quarter months from quarter number
quarter_months <- function(q) {
  switch(as.character(q), "1" = 1:3, "2" = 4:6, "3" = 7:9, "4" = 10:12)
}

# Define Functions ------------------------------------------------------------------------------------------

# Helper: get all admin_area columns for a given level (includes parents)
# e.g., "admin_area_4" -> c("admin_area_2", "admin_area_3", "admin_area_4")
geo_columns <- function(geo_level) {
  level_num <- as.integer(sub("admin_area_", "", geo_level))
  paste0("admin_area_", 2:level_num)
}

aggregate_to_quarter <- function(data, geo_level) {
  geo_cols <- geo_columns(geo_level)
  target_quarter_id <- as.numeric(paste0(SCORECARD_YEAR, sprintf("%02d", SCORECARD_QUARTER)))

  if (!geo_level %in% names(data)) {
    stop(sprintf("ERROR: Column '%s' not found in data. Available: %s",
                 geo_level, paste(grep("admin_area", names(data), value = TRUE), collapse = ", ")))
  }

  if (!"year" %in% names(data)) {
    data <- data %>% mutate(
      year = as.integer(substr(as.character(period_id), 1, 4)),
      month = as.integer(substr(as.character(period_id), 5, 6)),
      quarter_id = as.integer(paste0(year, sprintf("%02d", ceiling(month / 3))))
    )
  }

  data %>%
    filter(year == SCORECARD_YEAR, quarter_id == target_quarter_id) %>%
    group_by(across(all_of(geo_cols)), indicator_common_id) %>%
    summarise(count = sum(.data[[SELECTED_COUNT_VARIABLE]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = indicator_common_id, values_from = count, values_fill = 0)
}

merge_population <- function(area_data, geo_level) {
  area_names <- unique(area_data[[geo_level]])

  # Aggregate population to the target geo level
  pop <- population %>%
    group_by(.data[[geo_level]], period_id) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

  # Try current year (latest month first)
  current_year_periods <- sprintf("%04d%02d", SCORECARD_YEAR, 12:1)
  pop_data <- NULL

  for (period in current_year_periods) {
    pop_data <- pop %>% filter(period_id == period, .data[[geo_level]] %in% area_names)
    if (nrow(pop_data) > 0) {
      message(sprintf("  Using population from period: %s", period))
      break
    }
  }

  # Fallback: average of last 3 months of previous year
  if (is.null(pop_data) || nrow(pop_data) == 0) {
    message("  No current year population. Trying Q4 of previous year...")
    prev_periods <- sprintf("%04d%02d", SCORECARD_YEAR - 1, 10:12)
    prev_data <- list()
    for (period in prev_periods) {
      period_data <- pop %>% filter(period_id == period, .data[[geo_level]] %in% area_names)
      if (nrow(period_data) > 0) prev_data[[period]] <- period_data
    }
    if (length(prev_data) > 0) {
      pop_data <- bind_rows(prev_data) %>%
        group_by(.data[[geo_level]]) %>%
        summarise(count = mean(count, na.rm = TRUE), .groups = "drop")
      message(sprintf("  Using average population from %d months of previous year", length(prev_data)))
    } else {
      stop("ERROR: No population data found for current or previous year!")
    }
  }

  # Check coverage
  missing <- setdiff(area_names, unique(pop_data[[geo_level]]))
  if (length(missing) > 0) {
    warning(sprintf("Population missing for %d areas: %s", length(missing), paste(head(missing, 5), collapse = ", ")))
  }

  pop_summary <- pop_data %>%
    group_by(.data[[geo_level]]) %>%
    summarise(total_population = first(count), .groups = "drop")

  result <- area_data %>% left_join(pop_summary, by = geo_level)
  message(sprintf("  Population merged for %d/%d areas", sum(!is.na(result$total_population)), nrow(area_data)))
  result
}

calculate_scorecard <- function(data, geo_cols) {
  has_col <- function(col_name) col_name %in% names(data)

  # Pre-compute stockout sum (can't use across() with local vars inside mutate)
  stockout_cols <- c("stockout_bcg", "stockout_hepb", "stockout_opv", "stockout_penta",
                     "stockout_rotavirus", "stockout_measles", "stockout_yellowfever",
                     "stockout_pcv", "stockout_ipv", "stockout_mena", "stockout_tt", "stockout_any")
  available_stockout <- intersect(stockout_cols, names(data))
  if (length(available_stockout) > 0) {
    data$.stockout_sum <- rowSums(data[available_stockout], na.rm = TRUE)
  }

  data %>%
    mutate(
      # B: ANC Coverage (1 Visit)
      anc_coverage_1_visit = if(has_col("anc1")) (anc1 / (total_population * PREGNANT_WOMEN_PCT)) * 100 else NA,

      # C: ANC4/ANC1 Ratio
      anc4_anc1_ratio = if(has_col("anc1") & has_col("anc4")) {
        ifelse(anc1 == 0, NA, (anc4 / anc1) * 100)
      } else NA,

      # D: Skilled Birth Attendance
      skilled_birth_attendance = if(has_col("delivery") & has_col("sba")) {
        ifelse(delivery == 0, NA, (sba / delivery) * 100)
      } else NA,

      # E: Uterotonics Coverage
      uterotonics_coverage = if(has_col("delivery") & has_col("uterotonics")) {
        ifelse(delivery == 0, NA, (uterotonics / delivery) * 100)
      } else NA,

      # F: Fistula per 10000 Deliveries
      fistula_per_1000_deliveries = if(has_col("delivery") & has_col("obstetric_fistula")) {
        ifelse(delivery == 0, NA, (obstetric_fistula / (delivery / 10000)))
      } else NA,

      # G: Newborn Resuscitation
      newborn_resuscitation = if(has_col("birth_asphyxia") & has_col("neonatal_resuscitation")) {
        ifelse(birth_asphyxia == 0, NA, (neonatal_resuscitation / birth_asphyxia) * 100)
      } else NA,

      # H: Postnatal Visits (Within 3 days) — sum of newborn 1d + 2-3d visits
      postnatal_visits_3d = if(has_col("live_births") & (has_col("pnc_newborn_1d") | has_col("pnc_newborn_2_3d") | has_col("pnc"))) {
        pnc_total <- if(has_col("pnc_newborn_1d") & has_col("pnc_newborn_2_3d")) {
          pnc_newborn_1d + pnc_newborn_2_3d
        } else if(has_col("pnc")) pnc else 0
        ifelse(live_births == 0, NA, (pnc_total / live_births) * 100)
      } else NA,

      # I: LBW KMC Coverage — uses actual LBW counts if available, falls back to 15% estimate
      lbw_kmc_coverage = if(has_col("kmc")) {
        lbw_total <- if(has_col("lbw_female") & has_col("lbw_male")) {
          lbw_female + lbw_male
        } else if(has_col("live_births")) live_births * 0.15 else NA
        ifelse(is.na(lbw_total) | lbw_total == 0, NA, (kmc / lbw_total) * 100)
      } else NA,

      # J: Birth Registration
      birth_registration = if(has_col("child_registration")) {
        (child_registration / (total_population * BIRTHS_PCT)) * 100
      } else NA,

      # K: Modern Contraceptive Use
      modern_contraceptive_use = if(has_col("modern_fp")) {
        (modern_fp / (total_population * WOMEN_15_49_PCT)) * 100
      } else NA,

      # L: Pneumonia Treatment
      pneumonia_antibiotic_treatment = if(has_col("pneumonia") & has_col("pneumonia_treatment")) {
        ifelse(pneumonia == 0, NA, (pneumonia_treatment / pneumonia) * 100)
      } else NA,

      # M: Diarrhea Treatment
      diarrhea_ors_zinc_treatment = if(has_col("diarrhoea") & has_col("ors_zinc")) {
        ifelse(diarrhoea == 0, NA, (ors_zinc / diarrhoea) * 100)
      } else NA,

      # N: IPTp3 Coverage
      iptp3_coverage = if(has_col("anc1") & has_col("iptp3")) {
        ifelse(anc1 == 0, NA, (iptp3 / anc1) * 100)
      } else NA,

      # O: Under-5 LLIN Coverage (uses fully_immunized as denominator per legacy Excel)
      under5_llin_coverage = if(has_col("llin") && has_col("fully_immunized")) {
        (llin / fully_immunized) * 100
      } else NA,

      # P: BCG Coverage
      bcg_coverage = if(has_col("bcg")) {
        (bcg / (total_population * BIRTHS_PCT)) * 100
      } else NA,

      # Q: Fully Immunized Coverage
      fully_immunized_coverage = if(has_col("fully_immunized")) {
        (fully_immunized / (total_population * BIRTHS_PCT)) * 100
      } else NA,

      # R: Vaccine Stockout — stockout flags / reporting facilities * 100
      vaccine_stockout_percentage = if(has_col(".stockout_sum") & has_col("nhmis_actual_reports")) {
        stockout_den <- nhmis_actual_reports + if(has_col("nhmis_v2013_actual_reports")) nhmis_v2013_actual_reports else 0
        ifelse(stockout_den == 0, NA, (.stockout_sum / stockout_den) * 100)
      } else NA,

      # S: Exclusive Breastfeeding
      exclusive_breastfeeding_rate = if(has_col("ebf")) {
        (ebf / (total_population * INFANTS_0_6M_PCT)) * 100
      } else NA,

      # T: Growth Monitoring
      growth_monitoring_coverage = if(has_col("nutrition_screening")) {
        (nutrition_screening / (total_population * 0.20 * 0.25)) * 100
      } else NA,

      # U: GBV Care Coverage
      gbv_care_coverage = if(has_col("gbv_cases") & has_col("gbv_care")) {
        ifelse(gbv_cases == 0, NA, (gbv_care / gbv_cases) * 100)
      } else NA,

      # V: Facility Utilisation (OPD per 100 person-years)
      facility_utilisation_opd = if(has_col("opd")) {
        (opd / (total_population * 0.25)) * 100
      } else NA,

      # W: Under-1 Fully Immunised (MCV1 proxy)
      under1_fully_immunised_mcv1 = if(has_col("measles1")) {
        (measles1 / (total_population * BIRTHS_PCT)) * 100
      } else NA,

      # X: NHMIS Reporting Rate — actual / expected * 100
      nhmis_reporting_rate_final = if(has_col("nhmis_actual_reports") & has_col("nhmis_expected_reports")) {
        ifelse(nhmis_expected_reports == 0, NA, (nhmis_actual_reports / nhmis_expected_reports) * 100)
      } else NA,

      # Y: NHMIS Timeliness — on_time / expected * 100
      nhmis_data_timeliness_final = if(has_col("nhmis_actual_reports_ontime") & has_col("nhmis_expected_reports")) {
        ifelse(nhmis_expected_reports == 0, NA, (nhmis_actual_reports_ontime / nhmis_expected_reports) * 100)
      } else NA
    ) %>%
    select(all_of(geo_cols), anc_coverage_1_visit:nhmis_data_timeliness_final) %>%
    mutate(across(where(is.numeric), ~round(.x, 2)))
}

create_total_row <- function(merged_data, geo_cols) {
  total_counts <- merged_data %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))
  for (gc in geo_cols) total_counts[[gc]] <- paste0("__", COUNTRY_ISO3)

  calculate_scorecard(total_counts, geo_cols) %>%
    mutate(across(where(is.numeric), ~round(.x, 2)))
}

convert_scorecard_to_long <- function(scorecard_wide, geo_cols) {
  quarter_id <- as.numeric(paste0(SCORECARD_YEAR, sprintf("%02d", SCORECARD_QUARTER)))

  scorecard_wide %>%
    pivot_longer(cols = -all_of(geo_cols), names_to = "indicator_common_id", values_to = "value") %>%
    mutate(quarter_id = quarter_id, year = SCORECARD_YEAR) %>%
    select(all_of(geo_cols), quarter_id, year, indicator_common_id, value) %>%
    arrange(.data[[tail(geo_cols, 1)]], indicator_common_id)
}

create_summary_table <- function(adjusted_data, population, geo_level) {
  geo_cols <- geo_columns(geo_level)
  q_months <- quarter_months(SCORECARD_QUARTER)
  target_periods <- as.integer(sprintf("%04d%02d", SCORECARD_YEAR, q_months))

  # Previous quarter for population fallback
  prev_q <- if (SCORECARD_QUARTER == 1) 4 else SCORECARD_QUARTER - 1
  prev_year <- if (SCORECARD_QUARTER == 1) SCORECARD_YEAR - 1 else SCORECARD_YEAR
  prev_periods <- sprintf("%04d%02d", prev_year, quarter_months(prev_q))

  # Aggregate population to the target geo level
  pop <- population %>%
    group_by(.data[[geo_level]], period_id) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

  # Get indicator data for target quarter
  q_indicators <- adjusted_data %>%
    filter(period_id %in% target_periods) %>%
    group_by(across(all_of(geo_cols)), indicator_common_id) %>%
    summarise(total_quarter = sum(.data[[SELECTED_COUNT_VARIABLE]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = indicator_common_id, values_from = total_quarter, values_fill = 0)

  # Get previous quarter population average
  prev_pop_data <- list()
  for (period in prev_periods) {
    period_data <- pop %>% filter(period_id == period)
    if (nrow(period_data) > 0) prev_pop_data[[period]] <- period_data
  }

  if (length(prev_pop_data) == 0) {
    warning("No previous quarter population found! Using latest available.")
    latest_period <- pop %>% filter(!is.na(count)) %>%
      summarise(max_period = max(period_id)) %>% pull(max_period)
    if (!is.na(latest_period)) {
      prev_pop_data[[as.character(latest_period)]] <- pop %>% filter(period_id == latest_period)
    }
  }

  if (length(prev_pop_data) == 0) {
    stop("ERROR: No population data available!")
  }

  pop_avg <- bind_rows(prev_pop_data) %>%
    group_by(.data[[geo_level]]) %>%
    summarise(avg_population_prev_quarter = mean(count, na.rm = TRUE), .groups = "drop")

  # Merge indicators with population
  summary_table <- q_indicators %>%
    left_join(pop_avg, by = geo_level) %>%
    select(all_of(geo_cols), avg_population_prev_quarter, everything())

  # Compute vaccine stockout percentage for debug
  stockout_cols <- c("stockout_bcg", "stockout_hepb", "stockout_opv", "stockout_penta",
                     "stockout_rotavirus", "stockout_measles", "stockout_yellowfever",
                     "stockout_pcv", "stockout_ipv", "stockout_mena", "stockout_tt", "stockout_any")
  avail_stockout <- intersect(stockout_cols, names(summary_table))
  if (length(avail_stockout) > 0 && "nhmis_actual_reports" %in% names(summary_table)) {
    summary_table$stockout_sum <- rowSums(summary_table[avail_stockout], na.rm = TRUE)
    den <- summary_table$nhmis_actual_reports +
      if ("nhmis_v2013_actual_reports" %in% names(summary_table)) summary_table$nhmis_v2013_actual_reports else 0
    summary_table$vaccine_stockout_pct <- ifelse(den == 0, NA, (summary_table$stockout_sum / den) * 100)
  }

  # Create total row
  total_row <- summary_table %>%
    summarise(
      avg_population_prev_quarter = sum(avg_population_prev_quarter, na.rm = TRUE),
      across(where(is.numeric) & !matches("avg_population_prev_quarter"), \(x) sum(x, na.rm = TRUE))
    )
  for (gc in geo_cols) total_row[[gc]] <- "__NATIONAL"

  bind_rows(summary_table, total_row) %>% arrange(.data[[tail(geo_cols, 1)]])
}

# ------------------- Main Execution ------------------------------------------------------------------------
message(sprintf("Calculating NHSS Scorecard for %d Q%d using %s...",
                SCORECARD_YEAR, SCORECARD_QUARTER, SELECTED_COUNT_VARIABLE))

# ===== ADMIN_AREA_2 =====
message("\n========== Processing admin_area_2 ==========")
if ("admin_area_2" %in% names(adjusted_data)) {
  geo_cols_a2 <- geo_columns("admin_area_2")

  agg_a2 <- aggregate_to_quarter(adjusted_data, "admin_area_2")
  merged_a2 <- merge_population(agg_a2, "admin_area_2")

  scorecard_wide_a2 <- calculate_scorecard(merged_a2, geo_cols_a2)
  total_wide_a2 <- create_total_row(merged_a2, geo_cols_a2)

  scorecard_a2 <- bind_rows(
    convert_scorecard_to_long(scorecard_wide_a2, geo_cols_a2),
    convert_scorecard_to_long(total_wide_a2, geo_cols_a2)
  ) %>% arrange(admin_area_2, indicator_common_id)

  summary_a2 <- create_summary_table(adjusted_data, population, "admin_area_2")

  write_csv(scorecard_a2, "M7_output_scorecard_admin_area_2.csv")
  write_csv(summary_a2, "M7_output_scorecard_summary_admin_area_2.csv")
  message(sprintf("  Wrote M7_output_scorecard_admin_area_2.csv (%d rows)", nrow(scorecard_a2)))
} else {
  dummy_scorecard_a2 <- data.frame(
    admin_area_2 = character(), quarter_id = integer(), year = integer(),
    indicator_common_id = character(), value = numeric()
  )
  dummy_summary_a2 <- data.frame(
    admin_area_2 = character(), avg_population_prev_quarter = numeric()
  )
  write_csv(dummy_scorecard_a2, "M7_output_scorecard_admin_area_2.csv")
  write_csv(dummy_summary_a2, "M7_output_scorecard_summary_admin_area_2.csv")
  message("  No admin_area_2 data — saved empty files")
}

# ===== ADMIN_AREA_3 =====
message("\n========== Processing admin_area_3 ==========")
if ("admin_area_3" %in% names(adjusted_data)) {
  geo_cols_a3 <- geo_columns("admin_area_3")

  agg_a3 <- aggregate_to_quarter(adjusted_data, "admin_area_3")
  merged_a3 <- merge_population(agg_a3, "admin_area_3")

  scorecard_wide_a3 <- calculate_scorecard(merged_a3, geo_cols_a3)
  total_wide_a3 <- create_total_row(merged_a3, geo_cols_a3)

  scorecard_a3 <- bind_rows(
    convert_scorecard_to_long(scorecard_wide_a3, geo_cols_a3),
    convert_scorecard_to_long(total_wide_a3, geo_cols_a3)
  ) %>% arrange(admin_area_3, indicator_common_id)

  summary_a3 <- create_summary_table(adjusted_data, population, "admin_area_3")

  write_csv(scorecard_a3, "M7_output_scorecard_admin_area_3.csv")
  write_csv(summary_a3, "M7_output_scorecard_summary_admin_area_3.csv")
  message(sprintf("  Wrote M7_output_scorecard_admin_area_3.csv (%d rows)", nrow(scorecard_a3)))
} else {
  dummy_scorecard_a3 <- data.frame(
    admin_area_2 = character(), admin_area_3 = character(),
    quarter_id = integer(), year = integer(),
    indicator_common_id = character(), value = numeric()
  )
  dummy_summary_a3 <- data.frame(
    admin_area_2 = character(), admin_area_3 = character(),
    avg_population_prev_quarter = numeric()
  )
  write_csv(dummy_scorecard_a3, "M7_output_scorecard_admin_area_3.csv")
  write_csv(dummy_summary_a3, "M7_output_scorecard_summary_admin_area_3.csv")
  message("  No admin_area_3 data — saved empty files")
}

# ===== ADMIN_AREA_4 =====
message("\n========== Processing admin_area_4 ==========")
if ("admin_area_4" %in% names(adjusted_data)) {
  geo_cols_a4 <- geo_columns("admin_area_4")

  agg_a4 <- aggregate_to_quarter(adjusted_data, "admin_area_4")
  merged_a4 <- merge_population(agg_a4, "admin_area_4")

  scorecard_wide_a4 <- calculate_scorecard(merged_a4, geo_cols_a4)
  total_wide_a4 <- create_total_row(merged_a4, geo_cols_a4)

  scorecard_a4 <- bind_rows(
    convert_scorecard_to_long(scorecard_wide_a4, geo_cols_a4),
    convert_scorecard_to_long(total_wide_a4, geo_cols_a4)
  ) %>% arrange(admin_area_4, indicator_common_id)

  summary_a4 <- create_summary_table(adjusted_data, population, "admin_area_4")

  write_csv(scorecard_a4, "M7_output_scorecard_admin_area_4.csv")
  write_csv(summary_a4, "M7_output_scorecard_summary_admin_area_4.csv")
  message(sprintf("  Wrote M7_output_scorecard_admin_area_4.csv (%d rows)", nrow(scorecard_a4)))
} else {
  dummy_scorecard_a4 <- data.frame(
    admin_area_2 = character(), admin_area_3 = character(), admin_area_4 = character(),
    quarter_id = integer(), year = integer(),
    indicator_common_id = character(), value = numeric()
  )
  dummy_summary_a4 <- data.frame(
    admin_area_2 = character(), admin_area_3 = character(), admin_area_4 = character(),
    avg_population_prev_quarter = numeric()
  )
  write_csv(dummy_scorecard_a4, "M7_output_scorecard_admin_area_4.csv")
  write_csv(dummy_summary_a4, "M7_output_scorecard_summary_admin_area_4.csv")
  message("  No admin_area_4 data — saved empty files")
}

message("\nDone. All geographic levels processed.")
