# Scorecard Parameters
SCORECARD_YEAR <- 2025
SCORECARD_QUARTER <- 1
COUNTRY_ISO3 <- "NGA"
SELECTED_COUNT_VARIABLE <- "count_final_none" #adjusted for outliers and completeness

# Geographic level parameters
# SCORECARD_GEO_LEVEL: which admin column becomes each scorecard row (the breakdown level)
# SCORECARD_FILTER_COL / _VALUE: optional — restrict to areas within a parent area
# Example: national scorecard broken down by state (default)
#   SCORECARD_GEO_LEVEL <- "admin_area_3"
#   SCORECARD_FILTER_COL <- NULL
# Example: Kano state scorecard broken down by LGA
#   SCORECARD_GEO_LEVEL <- "admin_area_4"
#   SCORECARD_FILTER_COL <- "admin_area_3"
#   SCORECARD_FILTER_VALUE <- "Kano"
#   (requires asset files with admin_area_2 values matching LGA names)
SCORECARD_GEO_LEVEL <- "admin_area_3"
SCORECARD_FILTER_COL <- NULL
SCORECARD_FILTER_VALUE <- NULL

# Population assumptions (quarterly rates)
PREGNANT_WOMEN_PCT <- 0.05 * 0.25
BIRTHS_PCT <- 0.04 * 0.25
WOMEN_15_49_PCT <- 0.22 * 0.25
CHILDREN_U5_PCT <- 0.12 * 0.25
INFANTS_0_6M_PCT <- 0.02

# Asset file paths
VACCINE_STOCKOUT_ASSET <- "vaccine_stockout_pct.csv"
NHMIS_REPORTING_ASSET <- "nhmis_2019_reporting_rate.csv"
NHMIS_TIMELINESS_ASSET <- "nhmis_data_timeliness.csv"
POPULATION_ASSET <- "total_population.csv"

#-------------------------------------------------------------------------------------------------------------
# CB - R code FASTR PROJECT
# Module: NATIONAL HEALTH SECTOR SCORECARD (NHSS)
# Last edit: 2026 Feb 17 (multi-level geo support)

# This module calculates the NHSS scorecard indicators for quarterly performance dialogues

# -------------------------- KEY OUTPUTS ----------------------------------------------------------------------
# FILE: M7_scorecard_combined.csv        # Combined state + national scorecard in long format
# FILE: M7_scorecard_summary_Q{Q}_{YEAR}.csv # Quarterly summary with population

# Load Required Libraries -----------------------------------------------------------------------------------
library(dplyr)
library(readr)
library(tidyr)

# Load Data ------------------------------------------------------------------------------------------------
adjusted_data <- read_csv("M2_adjusted_data_admin_area.csv")
vaccine_stockout <- read_csv(VACCINE_STOCKOUT_ASSET)
nhmis_reporting <- read_csv(NHMIS_REPORTING_ASSET)
nhmis_timeliness <- read_csv(NHMIS_TIMELINESS_ASSET)
population <- read_csv(POPULATION_ASSET)

# Define Functions ------------------------------------------------------------------------------------------
aggregate_to_state_quarter <- function(data) {
  # Convert quarter to YYYYQQ format
  target_quarter_id <- as.numeric(paste0(SCORECARD_YEAR, sprintf("%02d", SCORECARD_QUARTER)))

  message(sprintf("Looking for quarter_id: %d", target_quarter_id))
  message(sprintf("Scorecard geo level: %s", SCORECARD_GEO_LEVEL))

  # Validate geo level column exists
  if (!SCORECARD_GEO_LEVEL %in% names(data)) {
    stop(sprintf("ERROR: Column '%s' not found in data. Available: %s",
                 SCORECARD_GEO_LEVEL, paste(grep("admin_area", names(data), value = TRUE), collapse = ", ")))
  }

  # Apply parent area filter if set
  if (!is.null(SCORECARD_FILTER_COL) && !is.null(SCORECARD_FILTER_VALUE)) {
    if (!SCORECARD_FILTER_COL %in% names(data)) {
      stop(sprintf("ERROR: Filter column '%s' not found in data.", SCORECARD_FILTER_COL))
    }
    n_before <- nrow(data)
    data <- data %>% filter(.data[[SCORECARD_FILTER_COL]] == SCORECARD_FILTER_VALUE)
    message(sprintf("Filtered to %s = '%s': %d -> %d rows",
                    SCORECARD_FILTER_COL, SCORECARD_FILTER_VALUE, n_before, nrow(data)))
    if (nrow(data) == 0) {
      stop(sprintf("ERROR: No data after filtering %s = '%s'. Check spelling / available values.",
                   SCORECARD_FILTER_COL, SCORECARD_FILTER_VALUE))
    }
  }

  # Derive year and quarter_id from period_id if not already present
  if (!"year" %in% names(data)) {
    data <- data %>% mutate(
      year = as.integer(substr(as.character(period_id), 1, 4)),
      month = as.integer(substr(as.character(period_id), 5, 6)),
      quarter_id = as.integer(paste0(year, sprintf("%02d", ceiling(month / 3))))
    )
  }

  data %>%
    filter(year == SCORECARD_YEAR, quarter_id == target_quarter_id) %>%
    group_by(.data[[SCORECARD_GEO_LEVEL]], indicator_common_id) %>%
    summarise(count = sum(.data[[SELECTED_COUNT_VARIABLE]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = indicator_common_id, values_from = count, values_fill = 0) %>%
    rename(admin_area_2 = !!sym(SCORECARD_GEO_LEVEL))
}

merge_assets <- function(state_data) {
  # Calculate quarter periods dynamically
  quarter_months <- case_when(
    SCORECARD_QUARTER == 1 ~ c(1, 2, 3),
    SCORECARD_QUARTER == 2 ~ c(4, 5, 6),
    SCORECARD_QUARTER == 3 ~ c(7, 8, 9),
    SCORECARD_QUARTER == 4 ~ c(10, 11, 12)
  )
  
  quarter_periods <- sprintf("%04d%02d", SCORECARD_YEAR, quarter_months)
  message(sprintf("Looking for asset data in periods: %s", paste(quarter_periods, collapse = ", ")))
  
  # Get states from adjusted data
  adjusted_states <- unique(state_data$admin_area_2)
  message(sprintf("States in adjusted data (%d): %s", 
                  length(adjusted_states), paste(head(adjusted_states, 5), collapse = ", ")))
  
  # Early check for country mismatch
  all_asset_states <- unique(c(
    vaccine_stockout$admin_area_2,
    nhmis_reporting$admin_area_2, 
    nhmis_timeliness$admin_area_2,
    population$admin_area_2
  ))
  
  overlap_states <- intersect(adjusted_states, all_asset_states)
  
  if(length(overlap_states) == 0) {
    stop(sprintf("ERROR: No state overlap between adjusted data and asset files!\n  Adjusted data states: %s\n  Asset file states: %s\n  Check if you're using asset files for the correct country.",
                 paste(head(adjusted_states, 3), collapse = ", "),
                 paste(head(all_asset_states, 3), collapse = ", ")))
  }
  
  if(length(overlap_states) < length(adjusted_states) * 0.5) {
    warning(sprintf("WARNING: Low state overlap (%d/%d). Possible country mismatch?", 
                    length(overlap_states), length(adjusted_states)))
  }
  
  # Validate asset files against adjusted data
  validate_asset <- function(asset_data, asset_name) {
    if(nrow(asset_data) == 0) {
      stop(sprintf("ERROR: %s file is empty!", asset_name))
    }
    
    # Check periods
    asset_periods <- unique(asset_data$period_id)
    available_periods <- intersect(quarter_periods, asset_periods)
    
    if(length(available_periods) == 0) {
      warning(sprintf("WARNING: %s: No data for target periods %s. Available: %s", 
                      asset_name, paste(quarter_periods, collapse = ", "), 
                      paste(head(asset_periods, 5), collapse = ", ")))
      return(data.frame(admin_area_2 = character(0), count = numeric(0)))
    }
    
    # Check states
    asset_states <- unique(asset_data$admin_area_2)
    missing_states <- setdiff(adjusted_states, asset_states)
    extra_states <- setdiff(asset_states, adjusted_states)
    
    if(length(missing_states) > 0) {
      warning(sprintf("WARNING: %s: Missing states - %s", 
                      asset_name, paste(head(missing_states, 5), collapse = ", ")))
    }
    
    if(length(extra_states) > 0) {
      message(sprintf("INFO: %s: Extra states (will be ignored) - %s", 
                      asset_name, paste(head(extra_states, 3), collapse = ", ")))
    }
    
    # Return filtered data
    asset_data %>% 
      filter(period_id %in% available_periods, admin_area_2 %in% adjusted_states) %>%
      group_by(admin_area_2) %>% 
      summarise(count = mean(count, na.rm = TRUE), .groups = "drop")
  }
  
  # Validate and process each asset
  vaccine_data <- validate_asset(vaccine_stockout, "Vaccine Stockout")
  reporting_data <- validate_asset(nhmis_reporting, "NHMIS Reporting") 
  timeliness_data <- validate_asset(nhmis_timeliness, "NHMIS Timeliness")
  
  # Handle population with fallback logic
  current_year_periods <- sprintf("%04d%02d", SCORECARD_YEAR, c(12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1))
  
  pop_data <- NULL
  found_current_year <- FALSE
  
  # Try current year first
  for(period in current_year_periods) {
    pop_data <- population %>% filter(period_id == period, admin_area_2 %in% adjusted_states)
    if(nrow(pop_data) > 0) {
      message(sprintf("Using population data from current year period: %s", period))
      found_current_year <- TRUE
      break
    }
  }
  
  # If no current year data found, use average of last 3 months of previous year
  if(!found_current_year || is.null(pop_data) || nrow(pop_data) == 0) {
    message("No current year population data found. Trying last 3 months of previous year...")
    
    prev_year_periods <- sprintf("%04d%02d", SCORECARD_YEAR - 1, c(10, 11, 12))
    prev_year_data_list <- list()
    
    for(period in prev_year_periods) {
      period_data <- population %>% filter(period_id == period, admin_area_2 %in% adjusted_states)
      if(nrow(period_data) > 0) {
        prev_year_data_list[[period]] <- period_data
        message(sprintf("Found population data for period: %s", period))
      }
    }
    
    if(length(prev_year_data_list) > 0) {
      combined_prev_year <- bind_rows(prev_year_data_list, .id = "period_source")
      
      pop_data <- combined_prev_year %>%
        group_by(admin_area_2) %>%
        summarise(
          count = mean(count, na.rm = TRUE),
          period_id = paste("AVG", paste(unique(period_source), collapse = "-"), sep = "_"),
          .groups = "drop"
        )
      
      message(sprintf("Using average population from %d months of previous year (%s)", 
                      length(prev_year_data_list), 
                      paste(names(prev_year_data_list), collapse = ", ")))
    } else {
      stop("ERROR: No population data found for current year or last 3 months of previous year!")
    }
  }
  
  if(is.null(pop_data) || nrow(pop_data) == 0) {
    stop("ERROR: No population data found!")
  }
  
  # Check population coverage
  pop_states <- unique(pop_data$admin_area_2)
  missing_pop_states <- setdiff(adjusted_states, pop_states)
  
  if(length(missing_pop_states) > 0) {
    stop(sprintf("ERROR: Population data missing for states: %s", 
                 paste(missing_pop_states, collapse = ", ")))
  }
  
  pop_summary <- pop_data %>% 
    group_by(admin_area_2) %>% 
    summarise(total_population = first(count), .groups = "drop")
  
  # Create asset list with renamed columns
  assets_final <- list(
    vaccine_data %>% rename(vaccine_stockout_pct = count),
    reporting_data %>% rename(nhmis_reporting_rate = count),
    timeliness_data %>% rename(nhmis_timeliness = count),
    pop_summary
  )
  
  # Join everything
  result <- Reduce(function(x, y) left_join(x, y, by = "admin_area_2"), c(list(state_data), assets_final))
  
  # Final validation
  states_with_all_assets <- result %>%
    filter(!is.na(total_population)) %>%
    nrow()
  
  message(sprintf("Successfully merged assets for %d/%d states", 
                  states_with_all_assets, nrow(state_data)))
  
  if(states_with_all_assets < nrow(state_data)) {
    warning("WARNING: Some states missing asset data - scorecard may be incomplete")
  }
  
  return(result)
}

calculate_scorecard <- function(data) {
  # Helper function to check if column exists
  has_col <- function(col_name) col_name %in% names(data)
  
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
      
      # H: Postnatal Visits (Within 3 days)
      postnatal_visits_3d = if(has_col("live_births") & has_col("pnc")) {
        ifelse(live_births == 0, NA, (pnc / live_births) * 100)
      } else NA,
      
      # I: LBW KMC Coverage - ESTIMATION (lbw female and male not available)
      lbw_kmc_coverage = if(has_col("kmc") & has_col("live_births")) {
        (kmc / (live_births * 0.15)) * 100
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
      
      # N: IPTP3 Coverage
      iptp3_coverage = if(has_col("anc1") & has_col("iptp3")) {
        ifelse(anc1 == 0, NA, (iptp3 / anc1) * 100)
      } else NA,

      # O: Under-5 LLIN Coverage (reproducing Excel formula - uses fully_immunized as denominator)
      # NOTE: Logically should be (llin / (total_population * CHILDREN_U5_PCT)) * 100 for actual coverage rate
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
      
      # R: Vaccine Stockout (from assets)
      vaccine_stockout_percentage = if(has_col("vaccine_stockout_pct")) vaccine_stockout_pct else NA,
      
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

      # X: NHMIS Reporting Rate (from assets)
      nhmis_reporting_rate_final = if(has_col("nhmis_reporting_rate")) nhmis_reporting_rate else NA,

      # Y: NHMIS Timeliness (from assets)
      nhmis_data_timeliness_final = if(has_col("nhmis_timeliness")) nhmis_timeliness else NA
    ) %>%
    select(admin_area_2, anc_coverage_1_visit:nhmis_data_timeliness_final) %>%
    mutate(across(where(is.numeric), ~round(.x, 2)))
}

create_total_row <- function(state_data, merged_data) {
  # Label: parent area value if filtered, otherwise country ISO3
  total_label <- if (!is.null(SCORECARD_FILTER_VALUE)) SCORECARD_FILTER_VALUE else COUNTRY_ISO3

  # Aggregate raw counts across all areas
  total_counts <- merged_data %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
    mutate(admin_area_2 = total_label)

  # Calculate indicators from totals
  calculate_scorecard(total_counts) %>%
    mutate(
      # Override asset indicators with means (can't sum percentages)
      vaccine_stockout_percentage = mean(state_data$vaccine_stockout_percentage, na.rm = TRUE),
      nhmis_reporting_rate_final = mean(state_data$nhmis_reporting_rate_final, na.rm = TRUE),
      nhmis_data_timeliness_final = mean(state_data$nhmis_data_timeliness_final, na.rm = TRUE)
    ) %>%
    mutate(across(where(is.numeric), ~round(.x, 2)))
}

convert_scorecard_to_long <- function(scorecard_wide) {
  quarter_id <- as.numeric(paste0(SCORECARD_YEAR, sprintf("%02d", SCORECARD_QUARTER)))
  
  scorecard_wide %>%
    pivot_longer(
      cols = -admin_area_2,
      names_to = "indicator_common_id", 
      values_to = "value"
    ) %>%
    mutate(
      quarter_id = quarter_id,
      year = SCORECARD_YEAR
    ) %>%
    select(admin_area_2, quarter_id, year, indicator_common_id, value) %>%
    arrange(admin_area_2, indicator_common_id)
}

create_summary_table <- function(adjusted_data, population) {
  # Compute quarter and previous-quarter periods dynamically
  quarter_months <- switch(as.character(SCORECARD_QUARTER),
    "1" = c(1, 2, 3), "2" = c(4, 5, 6), "3" = c(7, 8, 9), "4" = c(10, 11, 12)
  )
  target_periods_num <- as.integer(sprintf("%04d%02d", SCORECARD_YEAR, quarter_months))
  target_periods_chr <- sprintf("%04d%02d", SCORECARD_YEAR, quarter_months)

  # Previous quarter for population fallback
  prev_q <- if (SCORECARD_QUARTER == 1) 4 else SCORECARD_QUARTER - 1
  prev_year <- if (SCORECARD_QUARTER == 1) SCORECARD_YEAR - 1 else SCORECARD_YEAR
  prev_months <- switch(as.character(prev_q),
    "1" = c(1, 2, 3), "2" = c(4, 5, 6), "3" = c(7, 8, 9), "4" = c(10, 11, 12)
  )
  prev_periods <- sprintf("%04d%02d", prev_year, prev_months)

  quarter_label <- sprintf("Q%d %d", SCORECARD_QUARTER, SCORECARD_YEAR)
  prev_quarter_label <- sprintf("Q%d %d", prev_q, prev_year)

  message(sprintf("Creating %s summary table...", quarter_label))
  message("Available columns in adjusted_data: ", paste(names(adjusted_data), collapse = ", "))

  # Apply parent area filter if set
  if (!is.null(SCORECARD_FILTER_COL) && !is.null(SCORECARD_FILTER_VALUE)) {
    adjusted_data <- adjusted_data %>% filter(.data[[SCORECARD_FILTER_COL]] == SCORECARD_FILTER_VALUE)
  }

  # Get indicator data for target quarter
  q_indicators <- adjusted_data %>%
    filter(period_id %in% target_periods_num) %>%
    group_by(.data[[SCORECARD_GEO_LEVEL]], indicator_common_id) %>%
    summarise(
      total_quarter = sum(.data[[SELECTED_COUNT_VARIABLE]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(admin_area_2 = !!sym(SCORECARD_GEO_LEVEL)) %>%
    pivot_wider(
      names_from = indicator_common_id,
      values_from = total_quarter,
      values_fill = 0
    )

  message(sprintf("Found %s data for %d states and %d indicators",
                  quarter_label, nrow(q_indicators), ncol(q_indicators) - 1))

  # Get previous quarter population average
  prev_pop_data <- list()

  for(period in prev_periods) {
    period_data <- population %>% filter(period_id == period)
    if(nrow(period_data) > 0) {
      prev_pop_data[[period]] <- period_data
      message(sprintf("Found %s population data for period: %s", prev_quarter_label, period))
    }
  }

  if(length(prev_pop_data) == 0) {
    warning(sprintf("No %s population data found! Using latest available.", prev_quarter_label))
    latest_pop <- population %>%
      filter(!is.na(count)) %>%
      arrange(desc(period_id)) %>%
      slice_head(n = 1) %>%
      pull(period_id) %>%
      unique()

    if(length(latest_pop) > 0) {
      prev_pop_data[[latest_pop[1]]] <- population %>% filter(period_id == latest_pop[1])
      message(sprintf("Using fallback population data from: %s", latest_pop[1]))
    }
  }

  # Calculate average population
  if(length(prev_pop_data) > 0) {
    pop_avg <- bind_rows(prev_pop_data, .id = "period_source") %>%
      group_by(admin_area_2) %>%
      summarise(
        avg_population_prev_quarter = mean(count, na.rm = TRUE),
        .groups = "drop"
      )

    message(sprintf("Calculated %s population average from %d months",
                    prev_quarter_label, length(prev_pop_data)))
  } else {
    stop(sprintf("ERROR: No population data available for %s or fallback!", prev_quarter_label))
  }

  # Merge indicators with population
  summary_table <- q_indicators %>%
    left_join(pop_avg, by = "admin_area_2") %>%
    select(admin_area_2, avg_population_prev_quarter, everything())

  # Create total row (national if unfiltered, parent area if filtered)
  total_label <- if (!is.null(SCORECARD_FILTER_VALUE)) paste0("__TOTAL_", SCORECARD_FILTER_VALUE) else "__NATIONAL"
  total_row <- summary_table %>%
    summarise(
      admin_area_2 = total_label,
      avg_population_prev_quarter = sum(avg_population_prev_quarter, na.rm = TRUE),
      across(where(is.numeric) & !matches("avg_population_prev_quarter"), \(x) sum(x, na.rm = TRUE))
    )

  # Combine state and national data
  final_summary_table <- bind_rows(summary_table, total_row) %>%
    arrange(admin_area_2)

  message(sprintf("Summary table created with %d states + national and %d indicators",
                  nrow(summary_table),
                  ncol(summary_table) - 2))

  return(final_summary_table)
}

# Function to generate CREATE TABLE SQL
generate_sql_schema <- function(table_name, columns) {
  sql_lines <- paste0("  ", columns)
  sql <- c(
    paste0("CREATE TABLE ", table_name, " ("),
    paste(sql_lines, collapse = ",\n"),
    ");\n"
  )
  return(paste(sql, collapse = "\n"))
}

# ------------------- Main Execution ------------------------------------------------------------------------
geo_desc <- if (!is.null(SCORECARD_FILTER_VALUE)) {
  sprintf("%s within %s = '%s'", SCORECARD_GEO_LEVEL, SCORECARD_FILTER_COL, SCORECARD_FILTER_VALUE)
} else {
  sprintf("%s (national)", SCORECARD_GEO_LEVEL)
}
message(sprintf("Calculating NHSS Scorecard for %d Q%d at %s using %s...",
                SCORECARD_YEAR, SCORECARD_QUARTER, geo_desc, SELECTED_COUNT_VARIABLE))

# Aggregate facility data to state-quarter
state_aggregated <- aggregate_to_state_quarter(adjusted_data)

if(nrow(state_aggregated) == 0) {
  stop(sprintf("ERROR: No facility data found for %d Q%d!", SCORECARD_YEAR, SCORECARD_QUARTER))
}

if(ncol(state_aggregated) <= 1) {
  stop(sprintf("ERROR: No indicators found in facility data for %d Q%d!", SCORECARD_YEAR, SCORECARD_QUARTER))
}

message(sprintf("Aggregated %d indicators across %d areas (%s)",
                ncol(state_aggregated) - 1, nrow(state_aggregated), SCORECARD_GEO_LEVEL))

# Merge with assets
merged_data <- merge_assets(state_aggregated)

if(!"total_population" %in% names(merged_data)) {
  stop("ERROR: Population data not found - check POPULATION_ASSET file")
}

if(sum(!is.na(merged_data$total_population)) == 0) {
  stop("ERROR: All population values are NA - check population data period")
}

message("Asset data merged successfully")

# Calculate scorecards
state_scorecard_wide <- calculate_scorecard(merged_data)
if(nrow(state_scorecard_wide) == 0) {
  stop("ERROR: Scorecard calculation failed - no output generated")
}
message("State scorecard calculated")

total_scorecard_wide <- create_total_row(state_scorecard_wide, merged_data)
message("Total row scorecard calculated")

# Convert to long format and combine
state_scorecard <- convert_scorecard_to_long(state_scorecard_wide)
total_label <- if (!is.null(SCORECARD_FILTER_VALUE)) paste0("__TOTAL_", SCORECARD_FILTER_VALUE) else "__NATIONAL"
total_scorecard <- convert_scorecard_to_long(total_scorecard_wide) %>%
  mutate(admin_area_2 = total_label)

combined_scorecard <- bind_rows(state_scorecard, total_scorecard) %>%
  arrange(admin_area_2, indicator_common_id)

# Create quarterly summary table
q_summary <- create_summary_table(adjusted_data, population)

# Save output files — append filter value to filename if running sub-national
file_suffix <- if (!is.null(SCORECARD_FILTER_VALUE)) paste0("_", gsub(" ", "_", SCORECARD_FILTER_VALUE)) else ""
write_csv(combined_scorecard, sprintf("M7_scorecard_combined%s.csv", file_suffix))
write_csv(q_summary, sprintf("M7_scorecard_summary_Q%d_%d%s.csv", SCORECARD_QUARTER, SCORECARD_YEAR, file_suffix))
