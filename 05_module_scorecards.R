# Scorecard Parameters
SCORECARD_YEAR <- 2025
SCORECARD_QUARTER <- 1
SELECTED_COUNT_VARIABLE <- "count_final_both" #adjusted for outliers and completeness

# Population assumptions (quarterly rates)
PREGNANT_WOMEN_PCT <- 0.05 * 0.25
BIRTHS_PCT <- 0.04 * 0.25  
WOMEN_15_49_PCT <- 0.22 * 0.25
CHILDREN_U5_PCT <- 0.12 * 0.25
INFANTS_0_6M_PCT <- 0.015

# Asset file paths
VACCINE_STOCKOUT_ASSET <- "/Users/claireboulange/Desktop/Nigeria_RMNCAH_Additional_Data/scorecard_assets/vaccine_stockout_pct.csv"
NHMIS_REPORTING_ASSET <- "/Users/claireboulange/Desktop/Nigeria_RMNCAH_Additional_Data/scorecard_assets/nhmis_2019_reporting_rate.csv"
NHMIS_TIMELINESS_ASSET <- "/Users/claireboulange/Desktop/Nigeria_RMNCAH_Additional_Data/scorecard_assets/nhmis_data_timeliness.csv"
POPULATION_ASSET <- "/Users/claireboulange/Desktop/Nigeria_RMNCAH_Additional_Data/scorecard_assets/total_population.csv"

#-------------------------------------------------------------------------------------------------------------
# CB - R code FASTR PROJECT
# Module: RMNCAH SCORECARD CALCULATION
# Last edit: 2025 July 16

# This module calculates the 25 RMNCAH scorecard indicators (as per the Excel template)

# -------------------------- KEY OUTPUTS ----------------------------------------------------------------------
# FILE: M5_scorecard_state.csv           # State-level scorecard (25 indicators) in long format
# FILE: M5_scorecard_national.csv        # National scorecard (25 indicators) in long format

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
    group_by(admin_area_2, indicator_common_id) %>%
    summarise(count = sum(.data[[SELECTED_COUNT_VARIABLE]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = indicator_common_id, values_from = count, values_fill = 0)
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
  pop_periods_try <- c(
    sprintf("%04d%02d", SCORECARD_YEAR, c(12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1)),
    sprintf("%04d12", SCORECARD_YEAR - 1)
  )
  
  pop_data <- NULL
  for(period in pop_periods_try) {
    pop_data <- population %>% filter(period_id == period, admin_area_2 %in% adjusted_states)
    if(nrow(pop_data) > 0) {
      message(sprintf("Using population data from period: %s", period))
      break
    }
  }
  
  if(is.null(pop_data) || nrow(pop_data) == 0) {
    stop("ERROR: No population data found for any period or states!")
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
  data %>%
    mutate(
      # B: ANC Coverage (1 Visit)
      anc_coverage_1_visit = (anc1 / (total_population * PREGNANT_WOMEN_PCT)) * 100,
      # C: ANC4/ANC1 Ratio  
      anc4_anc1_ratio = ifelse(anc1 == 0, NA, (anc4 / anc1) * 100),
      # D: Skilled Birth Attendance
      skilled_birth_attendance = ifelse(deliveries == 0, NA, (sba / deliveries) * 100),
      # E: Uterotonics Coverage
      uterotonics_coverage = ifelse(deliveries == 0, NA, (uterotonics / deliveries) * 100),
      # F: Fistula per 1000 Deliveries
      fistula_per_1000_deliveries = ifelse(deliveries == 0, NA, (obstetric_fistula / (deliveries / 10000))),
      # G: Newborn Resuscitation
      newborn_resuscitation = ifelse(birth_asphyxia == 0, NA, (neonatal_resuscitation / birth_asphyxia) * 100),
      # H: Postnatal Visits (3 days)
      postnatal_visits_3d = ifelse(live_births == 0, NA, ((pnc_1d + pnc_2_3d) / live_births) * 100),
      # I: LBW KMC Coverage
      lbw_kmc_coverage = (kmc / (live_births * 0.15)) * 100,
      # J: Birth Registration
      birth_registration = (child_registration / (total_population * BIRTHS_PCT)) * 100,
      # K: Modern Contraceptive Use
      modern_contraceptive_use = (modern_fp / (total_population * WOMEN_15_49_PCT)) * 100,
      # L: Pneumonia Treatment
      pneumonia_antibiotic_treatment = ifelse(pneumonia == 0, NA, (pneumonia_treatment / pneumonia) * 100),
      # M: Diarrhea Treatment
      diarrhea_ors_zinc_treatment = ifelse(diarrhoea == 0, NA, (ors_zinc / diarrhoea) * 100),
      # N: IPTP3 Coverage
      iptp3_coverage = ifelse(anc1 == 0, NA, (iptp3 / anc1) * 100),
      # O: Malaria ACT Treatment
      malaria_act_treatment_rate = ifelse(malaria == 0, NA, (act_treatment / malaria) * 100),
      # P: Under-5 LLIN Coverage
      under5_llin_coverage = (llin / (total_population * CHILDREN_U5_PCT)) * 100,
      # Q: BCG Coverage
      bcg_coverage = (bcg / (total_population * BIRTHS_PCT)) * 100,
      # R: Penta3 Coverage
      penta3_coverage = (penta3 / (total_population * BIRTHS_PCT)) * 100,
      # S: Fully Immunized Coverage
      fully_immunized_coverage = (fully_immunized / (total_population * BIRTHS_PCT)) * 100,
      # T: Vaccine Stockout (from assets)
      vaccine_stockout_percentage = vaccine_stockout_pct,
      # U: Exclusive Breastfeeding
      exclusive_breastfeeding_rate = (ebf / (total_population * INFANTS_0_6M_PCT)) * 100,
      # V: Growth Monitoring
      growth_monitoring_coverage = (nutrition_screening / (total_population * CHILDREN_U5_PCT)) * 100,
      # W: GBV Care Coverage
      gbv_care_coverage = ifelse(gbv_cases == 0, NA, (gbv_care / gbv_cases) * 100),
      # X: NHMIS Reporting Rate (from assets)
      nhmis_reporting_rate_final = nhmis_reporting_rate,
      # Y: NHMIS Timeliness (from assets)
      nhmis_data_timeliness_final = nhmis_timeliness
    ) %>%
    select(admin_area_2, anc_coverage_1_visit:nhmis_data_timeliness_final) %>%
    mutate(across(where(is.numeric), ~round(.x, 2)))
}

create_national <- function(state_data, merged_data) {
  # Aggregate raw counts nationally
  national_counts <- merged_data %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
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
  # Create quarter_id in YYYYQQ format (e.g., 202501 for 2025 Q1)
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

# ------------------- Main Execution ------------------------------------------------------------------------
message(sprintf("Calculating RMNCAH Scorecard for %d Q%d using %s...", 
                SCORECARD_YEAR, SCORECARD_QUARTER, SELECTED_COUNT_VARIABLE))

# Aggregate facility data to state-quarter
state_aggregated <- aggregate_to_state_quarter(adjusted_data)

# Validate aggregated data
if(nrow(state_aggregated) == 0) {
  stop(sprintf("ERROR: No facility data found for %d Q%d! Check if data exists for this period.", 
               SCORECARD_YEAR, SCORECARD_QUARTER))
}

if(ncol(state_aggregated) <= 1) {
  stop(sprintf("ERROR: No indicators found in facility data for %d Q%d!", 
               SCORECARD_YEAR, SCORECARD_QUARTER))
}

message(sprintf("Aggregated %d indicators across %d states", 
                ncol(state_aggregated) - 1, nrow(state_aggregated)))

# Merge with assets
merged_data <- merge_assets(state_aggregated)

# Validate merged data
if(!"total_population" %in% names(merged_data)) {
  stop("ERROR: Population data not found - check POPULATION_ASSET file")
}

if(sum(!is.na(merged_data$total_population)) == 0) {
  stop("ERROR: All population values are NA - check population data period")
}

message("Asset data merged successfully")

# Calculate state scorecard (wide format)
state_scorecard_wide <- calculate_scorecard(merged_data)

# Validate scorecard
if(nrow(state_scorecard_wide) == 0) {
  stop("ERROR: Scorecard calculation failed - no output generated")
}

message("State scorecard calculated")

# Create national scorecard (wide format)
national_scorecard_wide <- create_national(state_scorecard_wide, merged_data)
message("National scorecard calculated")

# Convert to long format
state_scorecard <- convert_scorecard_to_long(state_scorecard_wide)
national_scorecard <- convert_scorecard_to_long(national_scorecard_wide)

# Save outputs in long format
write_csv(state_scorecard, "M5_scorecard_state.csv")
write_csv(national_scorecard, "M5_scorecard_national.csv")

message(sprintf("Scorecard completed: %d states + national level saved in long format", 
                length(unique(state_scorecard$admin_area_2))))