COUNTRY_ISO3 <- "NGA"
SELECTED_COUNT_VARIABLE <- "count_final_none"

# Population assumptions (quarterly rates)
BIRTHS_PCT <- 0.04 * 0.25
WOMEN_15_49_PCT <- 0.22 * 0.25
INTERPOLATE_POPULATION <- FALSE


#-------------------------------------------------------------------------------------------------------------
# CB - R code FASTR PROJECT
# Module: NATIONAL HEALTH SECTOR SCORECARD (NHSS)
# Last edit: 2026 Apr 10
#
# Calculates NHSS scorecard indicators at multiple geographic levels.
# All indicator data comes from the M2 adjusted data pipeline.
# Population is the only separate asset file (LGA-level denominator from DHIS2).
#
# Automatically detects which quarters have BOTH indicator data (all 3 months)
# AND population data, then produces scorecards for every such quarter.
#
# OUTPUTS (one per geo level):
#   M7_output_scorecard_{level}.csv         — Scorecard indicators (long format, areas + national)
#
# Load Required Libraries -----------------------------------------------------------------------------------
library(dplyr)
library(readr)
library(tidyr)

# Load Data ------------------------------------------------------------------------------------------------
ADJUSTED_DATA_FILE <- "M2_adjusted_data.csv"
POPULATION_FILE <- "total_population_NGA.csv"
adjusted_data <- read_csv(ADJUSTED_DATA_FILE, show_col_types = FALSE)
population <- read_csv(POPULATION_FILE, show_col_types = FALSE)
if ("indicator_common_id" %in% names(population)) {
  population <- population %>% filter(indicator_common_id == "total_population") %>% select(-indicator_common_id)
}
if ("admin_area_1" %in% names(population)) {
  population <- population %>% select(-admin_area_1)
}

# Derive nhmis_timely_and_data: on-time reports that also contain data
# John's definition: on-time AND (General Attendance > 0 OR OPD > 0 OR Penta3 > 0)
content_indicators <- intersect(c("gnl_attendance", "opd", "penta3"), unique(adjusted_data$indicator_common_id))
if (length(content_indicators) > 0 && "nhmis_actual_reports_ontime" %in% unique(adjusted_data$indicator_common_id)) {
  message("Deriving nhmis_timely_and_data (on-time + content check)...")
  cv <- SELECTED_COUNT_VARIABLE
  has_content <- adjusted_data %>%
    filter(indicator_common_id %in% content_indicators, !is.na(.data[[cv]]), .data[[cv]] > 0) %>%
    distinct(facility_id, period_id)
  ontime_rows <- adjusted_data %>%
    filter(indicator_common_id == "nhmis_actual_reports_ontime", !is.na(.data[[cv]]), .data[[cv]] > 0)
  timely_data_rows <- ontime_rows %>%
    semi_join(has_content, by = c("facility_id", "period_id")) %>%
    mutate(indicator_common_id = "nhmis_timely_and_data")
  message(sprintf("  %d/%d on-time reports also have content (%.0f%% empty filtered out)",
                  nrow(timely_data_rows), nrow(ontime_rows),
                  (1 - nrow(timely_data_rows) / nrow(ontime_rows)) * 100))
  adjusted_data <- bind_rows(adjusted_data, timely_data_rows)
}

# Detect Available Quarters ---------------------------------------------------------------------------------
# A quarter is "complete" for indicators if all 3 months have at least one row of data.
# A quarter has population if we can find pop data for that year (or fall back to previous year).

detect_indicator_quarters <- function(data) {
  periods <- unique(data$period_id)
  year_month <- data.frame(
    year = as.integer(substr(as.character(periods), 1, 4)),
    month = as.integer(substr(as.character(periods), 5, 6))
  )
  year_month$quarter <- ceiling(year_month$month / 3)
  # Count how many distinct months per year-quarter
  quarter_coverage <- year_month %>%
    distinct(year, quarter, month) %>%
    group_by(year, quarter) %>%
    summarise(n_months = n(), .groups = "drop") %>%
    filter(n_months == 3)
  quarter_coverage %>%
    mutate(quarter_id = as.integer(paste0(year, sprintf("%02d", quarter)))) %>%
    select(year, quarter, quarter_id)
}

detect_population_years <- function(pop) {
  periods <- unique(pop$period_id)
  unique(as.integer(substr(as.character(periods), 1, 4)))
}

indicator_quarters <- detect_indicator_quarters(adjusted_data)
pop_years <- detect_population_years(population)

# For each indicator quarter, check if population is available (same year OR previous year as fallback)
available_quarters <- indicator_quarters %>%
  filter(year %in% pop_years | (year - 1) %in% pop_years) %>%
  arrange(year, quarter)

if (nrow(available_quarters) == 0) {
  stop("ERROR: No quarters found with both complete indicator data and population data!")
}

message(sprintf("Found %d quarter(s) with both indicator and population data:", nrow(available_quarters)))
for (i in seq_len(nrow(available_quarters))) {
  message(sprintf("  %d Q%d", available_quarters$year[i], available_quarters$quarter[i]))
}

# Define Functions ------------------------------------------------------------------------------------------

# Helper: get all admin_area columns for a given level (includes parents)
# e.g., "admin_area_4" -> c("admin_area_2", "admin_area_3", "admin_area_4")
geo_columns <- function(geo_level) {
  level_num <- as.integer(sub("admin_area_", "", geo_level))
  paste0("admin_area_", 2:level_num)
}

aggregate_to_quarter <- function(data, geo_level, sc_year, sc_quarter) {
  geo_cols <- geo_columns(geo_level)
  target_quarter_id <- as.numeric(paste0(sc_year, sprintf("%02d", sc_quarter)))

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
    filter(year == sc_year, quarter_id == target_quarter_id) %>%
    group_by(across(all_of(geo_cols)), indicator_common_id) %>%
    summarise(count = sum(.data[[SELECTED_COUNT_VARIABLE]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = indicator_common_id, values_from = count, values_fill = 0)
}

interpolate_population <- function(pop_data, geo_level, sc_year, sc_quarter) {
  # Get one population value per area per year
  yearly_pop <- pop_data %>%
    mutate(pop_year = as.integer(substr(as.character(period_id), 1, 4))) %>%
    group_by(.data[[geo_level]], pop_year) %>%
    summarise(pop_value = first(count), .groups = "drop")

  # Quarter midpoint as fractional year (Q1=0.125, Q2=0.375, Q3=0.625, Q4=0.875)
  quarter_midpoints <- c(0.125, 0.375, 0.625, 0.875)
  t_target <- sc_year + quarter_midpoints[sc_quarter]

  # Anchor = mid-year (year + 0.5)
  # Q1-Q2: interpolate between (year-1) and year anchors
  # Q3-Q4: interpolate between year and (year+1) anchors
  if (sc_quarter <= 2) {
    year_a <- sc_year - 1
    year_b <- sc_year
  } else {
    year_a <- sc_year
    year_b <- sc_year + 1
  }

  available_years <- sort(unique(yearly_pop$pop_year))
  has_a <- year_a %in% available_years
  has_b <- year_b %in% available_years

  if (has_a && has_b) {
    anchor_a <- year_a + 0.5
    weight_b <- t_target - anchor_a  # anchor_b - anchor_a = 1.0

    pop_a <- yearly_pop %>% filter(pop_year == year_a) %>% select(.data[[geo_level]], pop_a = pop_value)
    pop_b <- yearly_pop %>% filter(pop_year == year_b) %>% select(.data[[geo_level]], pop_b = pop_value)

    result <- pop_a %>%
      inner_join(pop_b, by = geo_level) %>%
      mutate(total_population = pop_a + (pop_b - pop_a) * weight_b) %>%
      select(.data[[geo_level]], total_population)

    message(sprintf("  Interpolated population: %.1f%% weight on %d, %.1f%% on %d",
                    (1 - weight_b) * 100, year_a, weight_b * 100, year_b))
  } else if (has_b) {
    result <- yearly_pop %>% filter(pop_year == year_b) %>%
      select(.data[[geo_level]], total_population = pop_value)
    message(sprintf("  Edge quarter (no year %d data): using %d population flat", year_a, year_b))
  } else if (has_a) {
    result <- yearly_pop %>% filter(pop_year == year_a) %>%
      select(.data[[geo_level]], total_population = pop_value)
    message(sprintf("  Edge quarter (no year %d data): using %d population flat", year_b, year_a))
  } else {
    stop(sprintf("ERROR: No population data for years %d or %d!", year_a, year_b))
  }

  result
}

merge_population <- function(area_data, geo_level, sc_year, sc_quarter = NULL) {
  area_names <- unique(area_data[[geo_level]])

  # Aggregate population to the target geo level
  pop <- population %>%
    group_by(.data[[geo_level]], period_id) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

  # --- Interpolation path ---
  if (INTERPOLATE_POPULATION && !is.null(sc_quarter)) {
    pop_interp <- interpolate_population(pop, geo_level, sc_year, sc_quarter) %>%
      filter(.data[[geo_level]] %in% area_names)

    missing <- setdiff(area_names, unique(pop_interp[[geo_level]]))
    if (length(missing) > 0) {
      warning(sprintf("Population missing for %d areas: %s", length(missing), paste(head(missing, 5), collapse = ", ")))
    }

    result <- area_data %>% left_join(pop_interp, by = geo_level)
    message(sprintf("  Population (interpolated) merged for %d/%d areas",
                    sum(!is.na(result$total_population)), nrow(area_data)))
    return(result)
  }
  # --- End interpolation path ---

  # Try current year (latest month first)
  current_year_periods <- sprintf("%04d%02d", sc_year, 12:1)
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
    prev_periods <- sprintf("%04d%02d", sc_year - 1, 10:12)
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
      stop(sprintf("ERROR: No population data found for year %d or previous year!", sc_year))
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

  data %>%
    mutate(
      # 1: ANC4/ANC1 <20wks Ratio
      anc4_anc1_before20_ratio = if(has_col("anc1_before20") & has_col("anc4")) {
        ifelse(anc1_before20 == 0, NA, (anc4 / anc1_before20))
      } else NA,

      # 2: ANC4/ANC1 Ratio
      anc4_anc1_ratio = if(has_col("anc1") & has_col("anc4")) {
        ifelse(anc1 == 0, NA, (anc4 / anc1))
      } else NA,

      # 3: Skilled Birth Attendance
      skilled_birth_attendance = if(has_col("delivery") & has_col("sba")) {
        ifelse(delivery == 0, NA, (sba / delivery))
      } else NA,

      # 4: New FP Acceptors / Women of Reproductive Age
      new_fp_acceptors_rate = if(has_col("new_fp")) {
        (new_fp / (total_population * WOMEN_15_49_PCT))
      } else NA,

      # 5: ACT for Uncomplicated Malaria
      #    mal_treatment (ouzURM9c1FI) / mal_confirmed_uncomplicated (HdtaLx63988)
      #    Fallback: mal_treatment / mal_positive
      act_malaria_treatment = if(has_col("mal_treatment") & has_col("mal_confirmed_uncomplicated")) {
        ifelse(mal_confirmed_uncomplicated == 0, NA, (mal_treatment / mal_confirmed_uncomplicated))
      } else if(has_col("mal_treatment") & has_col("mal_positive")) {
        ifelse(mal_positive == 0, NA, (mal_treatment / mal_positive))
      } else NA,

      # 6: Penta3 Coverage
      penta3_coverage = if(has_col("penta3")) {
        (penta3 / (total_population * BIRTHS_PCT))
      } else NA,

      # 7: Fully Immunized Coverage
      fully_immunized_coverage = if(has_col("fully_immunized")) {
        (fully_immunized / (total_population * BIRTHS_PCT))
      } else NA,

      # 8: HTN New Cases per 10,000 person-years
      htn_new_per_10000 = if(has_col("hypertension_new")) {
        (hypertension_new / (total_population * 0.25)) * 10000
      } else NA,

      # 9: Diabetes New Cases per 10,000 person-years
      diabetes_new_per_10000 = if(has_col("diabetes_new")) {
        (diabetes_new / (total_population * 0.25)) * 10000
      } else NA,

      # 10: NHMIS Timeliness (on-time + content check / expected)
      nhmis_data_timeliness_final = if(has_col("nhmis_timely_and_data") & has_col("nhmis_expected_reports")) {
        ifelse(nhmis_expected_reports == 0, NA, (nhmis_timely_and_data / nhmis_expected_reports))
      } else if(has_col("nhmis_actual_reports_ontime") & has_col("nhmis_expected_reports")) {
        ifelse(nhmis_expected_reports == 0, NA, (nhmis_actual_reports_ontime / nhmis_expected_reports))
      } else NA
    ) %>%
    select(all_of(geo_cols), anc4_anc1_before20_ratio:nhmis_data_timeliness_final) %>%
    mutate(across(where(is.numeric), ~round(.x, 2)))
}


convert_scorecard_to_long <- function(scorecard_wide, geo_cols, sc_year, sc_quarter) {
  qid <- as.numeric(paste0(sc_year, sprintf("%02d", sc_quarter)))

  scorecard_wide %>%
    pivot_longer(cols = -all_of(geo_cols), names_to = "indicator_common_id", values_to = "value") %>%
    mutate(quarter_id = qid) %>%
    select(all_of(geo_cols), quarter_id, indicator_common_id, value) %>%
    arrange(.data[[tail(geo_cols, 1)]], indicator_common_id)
}

# ------------------- Main Execution ------------------------------------------------------------------------
# Process each geo level across ALL available quarters, then stack into one output file per level.

process_geo_level <- function(geo_level, data, empty_cols) {
  if (!geo_level %in% names(data)) {
    message(sprintf("  No %s data — saved empty file", geo_level))
    empty_df <- do.call(data.frame, setNames(
      lapply(empty_cols, function(type) if (type == "character") character() else if (type == "integer") integer() else numeric()),
      names(empty_cols)
    ))
    return(empty_df)
  }

  geo_cols <- geo_columns(geo_level)
  all_results <- list()

  for (i in seq_len(nrow(available_quarters))) {
    sc_year <- available_quarters$year[i]
    sc_quarter <- available_quarters$quarter[i]
    message(sprintf("  --- %d Q%d ---", sc_year, sc_quarter))

    agg <- aggregate_to_quarter(data, geo_level, sc_year, sc_quarter)
    if (nrow(agg) == 0) {
      message(sprintf("    No data for %s in %d Q%d, skipping", geo_level, sc_year, sc_quarter))
      next
    }
    merged <- merge_population(agg, geo_level, sc_year, sc_quarter)
    scorecard_wide <- calculate_scorecard(merged, geo_cols)
    scorecard_long <- convert_scorecard_to_long(scorecard_wide, geo_cols, sc_year, sc_quarter)
    all_results[[length(all_results) + 1]] <- scorecard_long
  }

  if (length(all_results) == 0) return(data.frame())
  bind_rows(all_results) %>% arrange(.data[[tail(geo_cols, 1)]], quarter_id, indicator_common_id)
}

message(sprintf("\nCalculating NHSS Scorecard for %d quarter(s) using %s...",
                nrow(available_quarters), SELECTED_COUNT_VARIABLE))

# ===== ADMIN_AREA_2 =====
message("\n========== Processing admin_area_2 ==========")
scorecard_a2 <- process_geo_level("admin_area_2", adjusted_data,
  c(admin_area_2 = "character", quarter_id = "integer",
    indicator_common_id = "character", value = "numeric"))
write_csv(scorecard_a2, "M7_output_scorecard_admin_area_2.csv")
message(sprintf("  Wrote M7_output_scorecard_admin_area_2.csv (%d rows)", nrow(scorecard_a2)))

# ===== ADMIN_AREA_3 =====
message("\n========== Processing admin_area_3 ==========")
scorecard_a3 <- process_geo_level("admin_area_3", adjusted_data,
  c(admin_area_2 = "character", admin_area_3 = "character", quarter_id = "integer",
    indicator_common_id = "character", value = "numeric"))
write_csv(scorecard_a3, "M7_output_scorecard_admin_area_3.csv")
message(sprintf("  Wrote M7_output_scorecard_admin_area_3.csv (%d rows)", nrow(scorecard_a3)))

# ===== ADMIN_AREA_4 =====
message("\n========== Processing admin_area_4 ==========")
scorecard_a4 <- process_geo_level("admin_area_4", adjusted_data,
  c(admin_area_2 = "character", admin_area_3 = "character", admin_area_4 = "character",
    quarter_id = "integer", indicator_common_id = "character", value = "numeric"))
write_csv(scorecard_a4, "M7_output_scorecard_admin_area_4.csv")
message(sprintf("  Wrote M7_output_scorecard_admin_area_4.csv (%d rows)", nrow(scorecard_a4)))

message(sprintf("\nDone. Processed %d quarter(s) across all geographic levels.", nrow(available_quarters)))
