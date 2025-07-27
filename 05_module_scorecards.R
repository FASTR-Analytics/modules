# Scorecard Parameters
SCORECARD_YEAR <- 2025
SCORECARD_QUARTER <- 1
SELECTED_COUNT_VARIABLE <- "count_final_none" #adjusted for outliers and completeness

# Population assumptions (quarterly rates)
PREGNANT_WOMEN_PCT <- 0.05 * 0.25
BIRTHS_PCT <- 0.04 * 0.25  
WOMEN_15_49_PCT <- 0.22 * 0.25
CHILDREN_U5_PCT <- 0.12 * 0.25
INFANTS_0_6M_PCT <- 0.015

# Asset file paths
VACCINE_STOCKOUT_ASSET <- "/Users/claireboulange/Desktop/unicef/vaccine_stockout_pct.csv"
NHMIS_REPORTING_ASSET <- "/Users/claireboulange/Desktop/unicef/nhmis_2019_reporting_rate.csv"
NHMIS_TIMELINESS_ASSET <- "/Users/claireboulange/Desktop/unicef/nhmis_data_timeliness.csv"
POPULATION_ASSET <- "/Users/claireboulange/Desktop/unicef//total_population.csv"

#-------------------------------------------------------------------------------------------------------------
# CB - R code FASTR PROJECT
# Module: RMNCAH SCORECARD CALCULATION
# Last edit: 2025 July 27

# This module calculates the 25 RMNCAH scorecard indicators (as per the Excel template)

# -------------------------- KEY OUTPUTS ----------------------------------------------------------------------
# FILE: M5_scorecard_combined.csv        # Combined state + national scorecard (25 indicators) in long format
# FILE: M5_scorecard_data_summary.csv # Q1 2025 summary with 2024 population
# FILE: M5_sql_schema_output.txt         # SQL schemas for both output files

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
  
  data %>%
    filter(year == SCORECARD_YEAR, quarter_id == target_quarter_id) %>%
    group_by(admin_area_3, indicator_common_id) %>%
    summarise(count = sum(.data[[SELECTED_COUNT_VARIABLE]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = indicator_common_id, values_from = count, values_fill = 0) %>%
    rename(admin_area_2 = admin_area_3)
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
  
  # Helper function to safely get column or return 0
  safe_col <- function(col_name) {
    if(has_col(col_name)) data[[col_name]] else 0
  }
  
  # Check which indicators we can calculate
  missing_cols <- c()
  if(!has_col("pnc_1d")) missing_cols <- c(missing_cols, "pnc_1d")
  if(!has_col("pnc_2_3d")) missing_cols <- c(missing_cols, "pnc_2_3d")
  
  if(length(missing_cols) > 0) {
    message("WARNING: Missing columns for scorecard calculation: ", paste(missing_cols, collapse = ", "))
    message("Affected indicators will be set to NA")
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
      
      # F: Fistula per 1000 Deliveries
      fistula_per_1000_deliveries = if(has_col("delivery") & has_col("obstetric_fistula")) {
        ifelse(delivery == 0, NA, (obstetric_fistula / (delivery / 10000)))
      } else NA,
      
      # G: Newborn Resuscitation
      newborn_resuscitation = if(has_col("birth_asphyxia") & has_col("neonatal_resuscitation")) {
        ifelse(birth_asphyxia == 0, NA, (neonatal_resuscitation / birth_asphyxia) * 100)
      } else NA,
      
      # H: Postnatal Visits (3 days)
      postnatal_visits_3d = if(has_col("live_births") & (has_col("pnc_1d") | has_col("pnc_2_3d"))) {
        pnc_total <- safe_col("pnc_1d") + safe_col("pnc_2_3d")
        ifelse(live_births == 0, NA, (pnc_total / live_births) * 100)
      } else NA,
      
      # I: LBW KMC Coverage - ESTIMATION
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
      
      # O: Malaria ACT Treatment
      malaria_act_treatment_rate = if(has_col("malaria") & has_col("act_treatment")) {
        ifelse(malaria == 0, NA, (act_treatment / malaria) * 100)
      } else NA,
      
      # P: Under-5 LLIN Coverage
      under5_llin_coverage = if(has_col("llin")) {
        (llin / (total_population * CHILDREN_U5_PCT)) * 100
      } else NA,
      
      # Q: BCG Coverage
      bcg_coverage = if(has_col("bcg")) {
        (bcg / (total_population * BIRTHS_PCT)) * 100
      } else NA,
      
      # R: Penta3 Coverage
      penta3_coverage = if(has_col("penta3")) {
        (penta3 / (total_population * BIRTHS_PCT)) * 100
      } else NA,
      
      # S: Fully Immunized Coverage
      fully_immunized_coverage = if(has_col("fully_immunized")) {
        (fully_immunized / (total_population * BIRTHS_PCT)) * 100
      } else NA,
      
      # T: Vaccine Stockout (from assets)
      vaccine_stockout_percentage = if(has_col("vaccine_stockout_pct")) vaccine_stockout_pct else NA,
      
      # U: Exclusive Breastfeeding
      exclusive_breastfeeding_rate = if(has_col("ebf")) {
        (ebf / (total_population * INFANTS_0_6M_PCT)) * 100
      } else NA,
      
      # V: Growth Monitoring
      growth_monitoring_coverage = if(has_col("nutrition_screening")) {
        (nutrition_screening / (total_population * CHILDREN_U5_PCT)) * 100
      } else NA,
      
      # W: GBV Care Coverage
      gbv_care_coverage = if(has_col("gbv_cases") & has_col("gbv_care")) {
        ifelse(gbv_cases == 0, NA, (gbv_care / gbv_cases) * 100)
      } else NA,
      
      # X: NHMIS Reporting Rate (from assets)
      nhmis_reporting_rate_final = if(has_col("nhmis_reporting_rate")) nhmis_reporting_rate else NA,
      
      # Y: NHMIS Timeliness (from assets)
      nhmis_data_timeliness_final = if(has_col("nhmis_timeliness")) nhmis_timeliness else NA
    ) %>%
    select(admin_area_2, anc_coverage_1_visit:nhmis_data_timeliness_final) %>%
    mutate(across(where(is.numeric), ~round(.x, 2)))
}

create_national <- function(state_data, merged_data) {
  # Aggregate raw counts nationally
  national_counts <- merged_data %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
    mutate(admin_area_2 = "Nigeria")
  
  # Calculate national indicators
  calculate_scorecard(national_counts) %>%
    mutate(
      # Override asset indicators with means
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

create_q1_summary_table <- function(adjusted_data, population) {
  message("Creating Q1 2025 summary table...")
  
  # Define Q1 2025 periods (Jan, Feb, March)
  q1_2025_periods <- c("202501", "202502", "202503")
  q4_2024_periods <- c("202410", "202411", "202412")
  
  message("Available columns in adjusted_data: ", paste(names(adjusted_data), collapse = ", "))
  
  # Get Q1 2025 indicator data
  q1_indicators <- adjusted_data %>%
    filter(period_id %in% c(202501, 202502, 202503)) %>%
    group_by(admin_area_3, indicator_common_id) %>%
    summarise(
      total_q1_2025 = sum(.data[[SELECTED_COUNT_VARIABLE]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(admin_area_2 = admin_area_3) %>%
    pivot_wider(
      names_from = indicator_common_id, 
      values_from = total_q1_2025, 
      values_fill = 0
    )
  
  message(sprintf("Found Q1 2025 data for %d states and %d indicators", 
                  nrow(q1_indicators), ncol(q1_indicators) - 1))
  
  # Get Q4 2024 population average
  q4_2024_pop_data <- list()
  
  for(period in q4_2024_periods) {
    period_data <- population %>% filter(period_id == period)
    if(nrow(period_data) > 0) {
      q4_2024_pop_data[[period]] <- period_data
      message(sprintf("Found Q4 2024 population data for period: %s", period))
    }
  }
  
  if(length(q4_2024_pop_data) == 0) {
    warning("No Q4 2024 population data found! Using latest available.")
    latest_pop <- population %>%
      filter(!is.na(count)) %>%
      arrange(desc(period_id)) %>%
      slice_head(n = 1) %>%
      pull(period_id) %>%
      unique()
    
    if(length(latest_pop) > 0) {
      q4_2024_pop_data[[latest_pop[1]]] <- population %>% filter(period_id == latest_pop[1])
      message(sprintf("Using fallback population data from: %s", latest_pop[1]))
    }
  }
  
  # Calculate average population
  if(length(q4_2024_pop_data) > 0) {
    pop_avg_q4_2024 <- bind_rows(q4_2024_pop_data, .id = "period_source") %>%
      group_by(admin_area_2) %>%
      summarise(
        avg_population_q4_2024 = mean(count, na.rm = TRUE),
        .groups = "drop"
      )
    
    message(sprintf("Calculated Q4 2024 population average from %d months", 
                    length(q4_2024_pop_data)))
  } else {
    stop("ERROR: No population data available for Q4 2024 or fallback!")
  }
  
  # Merge indicators with population
  summary_table <- q1_indicators %>%
    left_join(pop_avg_q4_2024, by = "admin_area_2") %>%
    select(admin_area_2, avg_population_q4_2024, everything())
  
  # Create national totals
  national_totals <- summary_table %>%
    summarise(
      admin_area_2 = "NATIONAL",
      avg_population_q4_2024 = sum(avg_population_q4_2024, na.rm = TRUE),
      across(where(is.numeric) & !matches("avg_population_q4_2024"), \(x) sum(x, na.rm = TRUE))
    )
  
  # Combine state and national data
  final_summary_table <- bind_rows(summary_table, national_totals) %>%
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
message(sprintf("Calculating RMNCAH Scorecard for %d Q%d using %s...", 
                SCORECARD_YEAR, SCORECARD_QUARTER, SELECTED_COUNT_VARIABLE))

# Aggregate facility data to state-quarter
state_aggregated <- aggregate_to_state_quarter(adjusted_data)

if(nrow(state_aggregated) == 0) {
  stop(sprintf("ERROR: No facility data found for %d Q%d!", SCORECARD_YEAR, SCORECARD_QUARTER))
}

if(ncol(state_aggregated) <= 1) {
  stop(sprintf("ERROR: No indicators found in facility data for %d Q%d!", SCORECARD_YEAR, SCORECARD_QUARTER))
}

message(sprintf("Aggregated %d indicators across %d states", 
                ncol(state_aggregated) - 1, nrow(state_aggregated)))

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

national_scorecard_wide <- create_national(state_scorecard_wide, merged_data)
message("National scorecard calculated")

# Convert to long format and combine
state_scorecard <- convert_scorecard_to_long(state_scorecard_wide)
national_scorecard <- convert_scorecard_to_long(national_scorecard_wide) %>%
  mutate(admin_area_2 = "NATIONAL")

combined_scorecard <- bind_rows(state_scorecard, national_scorecard) %>%
  arrange(admin_area_2, indicator_common_id)

# Create Q1 summary table
q1_summary <- create_q1_summary_table(adjusted_data, population)

# Save output files
write_csv(combined_scorecard, "M5_scorecard_combined.csv")
write_csv(q1_summary, "M5_Q1_2025_summary_with_2024_population.csv")

# -------------------------------- Generate SQL Schemas --------------------------------
# M5 scorecard schema
scorecard_cols <- c(
  "admin_area_2 TEXT NOT NULL",
  "quarter_id INTEGER NOT NULL", 
  "year INTEGER NOT NULL",
  "indicator_common_id TEXT NOT NULL",
  "value NUMERIC"
)

scorecard_sql <- generate_sql_schema("ro_m5_scorecard_combined_csv", scorecard_cols)

# Q1 summary schema
q1_summary_base_cols <- c(
  "admin_area_2 TEXT NOT NULL",
  "avg_population_q4_2024 NUMERIC NOT NULL"
)

if(exists("q1_summary") && nrow(q1_summary) > 0) {
  indicator_cols <- names(q1_summary)[!names(q1_summary) %in% c("admin_area_2", "avg_population_q4_2024")]
  indicator_sql_cols <- paste0(indicator_cols, " NUMERIC")
  q1_summary_all_cols <- c(q1_summary_base_cols, indicator_sql_cols)
} else {
  q1_summary_all_cols <- c(q1_summary_base_cols, "-- Add indicator columns as NUMERIC based on your data")
}

q1_summary_sql <- generate_sql_schema("ro_m5_scorecard_data_summary_csv", q1_summary_all_cols)

# Write combined SQL schemas
combined_sql_output <- c(
  "-- M5 Scorecard Combined Table",
  scorecard_sql,
  "",
  "-- M5 Scorecard Data Summary Table", 
  q1_summary_sql
)

writeLines(combined_sql_output, "M5_sql_schema_output.txt")

message(sprintf("Scorecard completed: %d states + national level saved", 
                length(unique(state_scorecard$admin_area_2))))
message("Files created:")
message("- M5_scorecard_combined.csv (states + national)")
message("- M5_scorecard_data_summary.csv (Q1 summary)")
message("- M5_sql_schema_output.txt (SQL schemas)")