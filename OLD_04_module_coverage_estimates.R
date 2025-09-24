SELECTED_COUNT_VARIABLE <- "count_final_both"  # Options: "count_final_none", "count_final_outlier", "count_final_completeness", "count_final_both"

CURRENT_YEAR <- as.numeric(format(Sys.Date(), "%Y"))  # Dynamically get current year
MIN_YEAR <- 2000  # Set a fixed minimum year for filtering

PREGNANCY_LOSS_RATE <- 0.03 
TWIN_RATE <- 0.015       
STILLBIRTH_RATE <- 0.02
P1_NMR <- 0.041      #Default = 0.03
P2_PNMR <- 0.022
INFANT_MORTALITY_RATE <- 0.063  #Default = 0.05


ANALYSIS_LEVEL <- "NATIONAL_PLUS_AA2"      # Options: "NATIONAL_ONLY", "NATIONAL_PLUS_AA2", "NATIONAL_PLUS_AA2_AA3"

#-------------------------------------------------------------------------------------------------------------
# CB - R code FASTR PROJECT
# Last edit: 2025 Sep 8
# Module: COVERAGE ESTIMATES
#
# ------------------------------ Load Required Libraries -----------------------------------------------------
library(dplyr)
library(tidyr)       
library(zoo)       
library(stringr)     
library(purrr)

# ------------------------------ Define File Paths -----------------------------
# Direct read from GitHub
PROJECT_DATA_COVERAGE <- "https://raw.githubusercontent.com/FASTR-Analytics/modules/refs/heads/main/survey_data_unified.csv"
PROJECT_DATA_POPULATION <- "https://raw.githubusercontent.com/FASTR-Analytics/modules/refs/heads/main/population_estimates_only.csv"



# Input Datasets
adjusted_volume_data <- read.csv("M2_adjusted_data_national.csv", fileEncoding = "UTF-8")
adjusted_volume_data_subnational <- read.csv("M2_adjusted_data_admin_area.csv", fileEncoding = "UTF-8")
survey_data_unified <- read.csv(PROJECT_DATA_COVERAGE, fileEncoding = "UTF-8")
population_estimates_only <- read.csv(PROJECT_DATA_POPULATION, fileEncoding = "UTF-8")

# ------------------------------ Prepare Data for Analysis -------------------------
# A flag to track if pnc1 was renamed to pnc1_mother
pnc1_renamed_to_mother <- FALSE

message("Analysis mode: ", ANALYSIS_LEVEL)

# Always run national analysis
survey_data_national <- survey_data_unified %>% filter(admin_area_2 == "NATIONAL")

# Initialize subnational variables
hmis_data_subnational <- NULL
survey_data_subnational <- NULL
combined_admin2_export <- NULL
combined_admin3_export <- NULL

# Validate and prepare subnational data based on ANALYSIS_LEVEL
if (ANALYSIS_LEVEL %in% c("NATIONAL_PLUS_AA2", "NATIONAL_PLUS_AA2_AA3")) {
  
  # First check: Do we have any subnational survey data?
  survey_data_subnational <- survey_data_unified %>% filter(admin_area_2 != "NATIONAL")
  
  if (nrow(survey_data_subnational) == 0) {
    warning("SAFEGUARD: No subnational survey data found. Falling back to: NATIONAL_ONLY")
    original_level <- ANALYSIS_LEVEL
    ANALYSIS_LEVEL <- "NATIONAL_ONLY"
    hmis_data_subnational <- NULL
    survey_data_subnational <- NULL
  } else {
    
    # Check if admin_area_3 data can be used for NATIONAL_PLUS_AA2_AA3
    if (ANALYSIS_LEVEL == "NATIONAL_PLUS_AA2_AA3") {
      has_admin3_hmis <- "admin_area_3" %in% names(adjusted_volume_data_subnational)
      
      if (has_admin3_hmis) {
        hmis_admin3_values <- adjusted_volume_data_subnational %>%
          filter(!is.na(admin_area_3) & admin_area_3 != "" & admin_area_3 != "ZONE") %>%
          distinct(admin_area_3) %>% pull(admin_area_3)
        
        survey_admin2_values <- survey_data_subnational %>%
          distinct(admin_area_2) %>% pull(admin_area_2)
        
        matching_areas <- intersect(hmis_admin3_values, survey_admin2_values)
        has_usable_admin3 <- length(matching_areas) > 0
        
        if (length(matching_areas) > 0) {
          message("✓ admin_area_3 validation passed: ", length(matching_areas), "/", length(hmis_admin3_values), " areas match")
        }
      } else {
        has_usable_admin3 <- FALSE
      }
      
      if (!has_usable_admin3) {
        warning("SAFEGUARD: admin_area_3 data not usable. Falling back to: NATIONAL_PLUS_AA2")
        original_level <- ANALYSIS_LEVEL
        ANALYSIS_LEVEL <- "NATIONAL_PLUS_AA2"
      }
    }
    
    # Prepare HMIS data based on final analysis level
    hmis_data_subnational <- adjusted_volume_data_subnational
  }
}

message("Final analysis level: ", ANALYSIS_LEVEL)

# ------------------------------ Rename for Test Instance -------------------------------
adjusted_volume_data <- adjusted_volume_data %>%
  mutate(admin_area_1 = case_when(
    admin_area_1 %in% c("Pays 001", "Country 001") ~ "Federal Govt of Somalia",
    TRUE ~ admin_area_1
  ))

# ------------------------------ Define Parameters --------------------------------
# Coverage Estimation Parameters
coverage_params <- list(
  indicators = c(
    # Core
    "anc1", "anc4", "delivery", "bcg",
    "penta1", "penta3",
    "measles1", "measles2",
    "rota1", "rota2",
    "opv1", "opv2", "opv3",
    "pnc1_mother",
    "nmr", "imr"
  )
)

# List of survey variables to carry forward (for forward-fill and projections)
survey_vars <- c(
  "avgsurvey_anc1", "avgsurvey_anc4", "avgsurvey_delivery",
  "avgsurvey_bcg",
  "avgsurvey_penta1", "avgsurvey_penta3",
  "avgsurvey_measles1", "avgsurvey_measles2",
  "avgsurvey_rota1", "avgsurvey_rota2",
  "avgsurvey_opv1", "avgsurvey_opv2", "avgsurvey_opv3",
  "avgsurvey_pnc1_mother",
  "avgsurvey_nmr", "avgsurvey_imr", "postnmr"
)

# ------------------------------ Define Functions --------------------------------
# Part 1 - prepare hmis data
process_hmis_adjusted_volume <- function(adjusted_volume_data, count_col = SELECTED_COUNT_VARIABLE) {
  
  expected_indicators <- c(
    # Core RMNCH indicators
    "anc1", "anc4", "delivery", "bcg", "penta1", "penta3", "nmr", "imr",
    "measles1", "measles2", "rota1", "rota2", "opv1", "opv2", "opv3", "pnc1_mother"
  )
  
  message("Loading and mapping adjusted HMIS volume...")
  
  # Check if both pnc1_mother and pnc1 exist in the data
  has_pnc1_mother <- "pnc1_mother" %in% adjusted_volume_data$indicator_common_id
  has_pnc1 <- "pnc1" %in% adjusted_volume_data$indicator_common_id
  
  if (!has_pnc1_mother && has_pnc1) {
    # If pnc1_mother doesn't exist but pnc1 does, rename pnc1 to pnc1_mother
    adjusted_volume_data$indicator_common_id[adjusted_volume_data$indicator_common_id == "pnc1"] <- "pnc1_mother"
    pnc1_renamed_to_mother <<- TRUE # Set flag to TRUE
  }
  
  
  has_admin2 <- "admin_area_2" %in% names(adjusted_volume_data)
  
  # Ensure year and month exist
  if (!all(c("year", "month") %in% names(adjusted_volume_data))) {
    adjusted_volume_data <- adjusted_volume_data %>%
      mutate(
        year = as.integer(substr(period_id, 1, 4)),
        month = as.integer(substr(period_id, 5, 6))
      )
  }
  
  group_vars <- if (has_admin2) c("admin_area_1", "admin_area_2", "year") else c("admin_area_1", "year")
  
  adjusted_volume <- adjusted_volume_data %>%
    mutate(count = .data[[count_col]]) %>%
    select(any_of(c("admin_area_1", "admin_area_2", "year", "month", "indicator_common_id", "count"))) %>%
    arrange(across(any_of(c("admin_area_1", "admin_area_2", "year", "month", "indicator_common_id"))))
  
  missing <- setdiff(expected_indicators, unique(adjusted_volume$indicator_common_id))
  if (length(missing) > 0) {
    warning("The following indicators are not available in the HMIS data: ", paste(missing, collapse = ", "))
  }
  
  hmis_countries <- unique(adjusted_volume$admin_area_1)
  message("HMIS data for country: ", paste(hmis_countries, collapse = ", "))
  
  nummonth_data <- adjusted_volume %>%
    distinct(across(all_of(c(group_vars, "month")))) %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(nummonth = n_distinct(month, na.rm = TRUE), .groups = "drop")
  
  message("Aggregating HMIS volume to annual level...")
  
  annual_hmis <- adjusted_volume %>%
    group_by(across(all_of(c(group_vars, "indicator_common_id")))) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = indicator_common_id,
      values_from = count,
      names_prefix = "count",
      values_fill = list(count = 0)
    ) %>%
    left_join(nummonth_data, by = group_vars) %>%
    arrange(across(all_of(group_vars)))
  
  list(
    annual_hmis = annual_hmis,
    hmis_countries = hmis_countries
  )
}

# Part 2 - prepare survey data - UPDATED HARMONIZATION
process_survey_data <- function(survey_data, hmis_countries, min_year = MIN_YEAR, max_year = CURRENT_YEAR) {
  
  # Harmonize indicator names used in survey to match HMIS format
  survey_data <- survey_data %>%
    mutate(indicator_common_id = recode(indicator_common_id,
                                        "polio1" = "opv1",
                                        "polio2" = "opv2",
                                        "polio3" = "opv3",
                                        "pnc1" = "pnc1_mother"
    ))
  
  is_national <- all(unique(survey_data$admin_area_2) == "NATIONAL")
  
  indicators <- c("anc1", "anc4", "delivery", "bcg", "penta1", "penta3", "measles1", "measles2",
                  "rota1", "rota2", "opv1", "opv2", "opv3", "pnc1_mother", "nmr", "imr")
  
  survey_filtered <- if (is_national) {
    survey_data %>% filter(admin_area_2 == "NATIONAL")
  } else {
    survey_data %>% filter(admin_area_2 != "NATIONAL")
  }
  
  survey_filtered <- survey_filtered %>%
    mutate(source = case_when(
      str_detect(tolower(source), "dhs")   ~ "dhs",
      str_detect(tolower(source), "mics")  ~ "mics",
      str_detect(tolower(source), "unwpp") ~ "unwpp",
      TRUE ~ tolower(source)
    ))
  
  raw_survey_values <- survey_filtered %>%
    filter(source %in% c("dhs", "mics")) %>%
    group_by(admin_area_1, admin_area_2, year, indicator_common_id) %>%
    arrange(factor(source, levels = c("dhs", "mics"))) %>%
    summarise(survey_value = first(survey_value), .groups = "drop") %>%
    filter(year >= min_year & year <= max_year) %>%
    pivot_wider(names_from = indicator_common_id,
                values_from = survey_value,
                names_glue = "rawsurvey_{indicator_common_id}")
  
  full_years <- seq(min_year, max_year)
  group_keys <- if (is_national) {
    c("admin_area_1", "indicator_common_id", "source")
  } else {
    c("admin_area_1", "admin_area_2", "indicator_common_id", "source")
  }
  
  survey_extended <- survey_filtered %>%
    filter(year %in% full_years) %>%
    group_by(across(all_of(group_keys)), .drop = FALSE) %>%
    group_modify(~ {
      if (nrow(.x) == 0) return(tibble())
      .x %>% complete(year = full_years) %>% arrange(year) %>%
        mutate(survey_value_carry = zoo::na.locf(survey_value, na.rm = FALSE))
    }) %>%
    ungroup()
  
  survey_wide <- survey_extended %>%
    select(all_of(c("admin_area_1", if (!is_national) "admin_area_2")),
           year, indicator_common_id, source, survey_value_carry) %>%
    pivot_wider(names_from = c(source, indicator_common_id),
                values_from = survey_value_carry,
                names_glue = "{indicator_common_id}_{source}",
                values_fn = mean)
  
  for (ind in indicators) {
    dhs_col <- paste0(ind, "_dhs")
    mics_col <- paste0(ind, "_mics")
    avg_col <- paste0("avgsurvey_", ind)
    if (dhs_col %in% names(survey_wide) || mics_col %in% names(survey_wide)) {
      survey_wide[[avg_col]] <- coalesce(survey_wide[[dhs_col]], survey_wide[[mics_col]])
    }
  }
  
  survey_wide <- survey_wide %>%
    mutate(postnmr = ifelse("avgsurvey_imr" %in% names(.) & "avgsurvey_nmr" %in% names(.),
                            avgsurvey_imr - avgsurvey_nmr, NA_real_))
  
  carry_group <- if (is_national) "admin_area_1" else c("admin_area_1", "admin_area_2")
  
  survey_carried <- survey_wide %>%
    group_by(across(all_of(carry_group))) %>%
    complete(year = full_seq(year, 1)) %>%
    arrange(across(all_of(carry_group)), year) %>%
    mutate(across(everything(), ~ zoo::na.locf(.x, na.rm = FALSE))) %>%
    ungroup()
  
  for (ind in c(indicators, "postnmr")) {
    avg_col <- paste0("avgsurvey_", ind)
    carry_col <- paste0(ind, "carry")
    if (avg_col %in% names(survey_carried)) {
      survey_carried[[carry_col]] <- survey_carried[[avg_col]]
    }
  }
  
  if (!is_national) {
    for (ind in c("bcg", "penta1", "penta3")) {
      carry_col <- paste0(ind, "carry")
      if (!(carry_col %in% names(survey_carried))) {
        survey_carried[[carry_col]] <- 0.70
      } else {
        survey_carried[[carry_col]] <- ifelse(is.na(survey_carried[[carry_col]]), 0.70,
                                              survey_carried[[carry_col]])
      }
    }
  }
  
  if (is_national) {
    survey_carried <- survey_carried %>% mutate(admin_area_2 = "NATIONAL")
    raw_survey_values <- raw_survey_values %>% mutate(admin_area_2 = "NATIONAL")
  }
  
  return(list(
    carried = survey_carried %>% arrange(across(any_of(c("admin_area_1", if (!is_national) "admin_area_2", "year")))),
    raw = raw_survey_values
  ))
}

#Part 2b - prepare unwpp data
process_national_population_data <- function(population_data, hmis_countries) {
  
  population_data %>%
    filter(admin_area_2 == "NATIONAL") %>%
    filter(admin_area_1 %in% hmis_countries) %>%
    mutate(source = tolower(source)) %>%
    select(admin_area_1, year, indicator_common_id, survey_value, source) %>%
    pivot_wider(
      names_from  = c(indicator_common_id, source),
      values_from = survey_value,
      names_glue  = "{indicator_common_id}_{source}",
      values_fn   = mean
    )
}

#Part 3 - calculate denominators
calculate_denominators <- function(hmis_data, survey_data, population_data = NULL) {
  if (!"nmrcarry" %in% names(survey_data)) {
    message("`nmrcarry` not found in survey data – filling with default from P1_NMR = ", P1_NMR)
    survey_data$nmrcarry <- P1_NMR
  }
  
  has_admin_area_2 <- "admin_area_2" %in% names(hmis_data)
  
  if (has_admin_area_2) {
    # Standard join admin_area_2 to admin_area_2 (survey data is restructured)
    data <- hmis_data %>%
      full_join(survey_data, by = c("admin_area_1", "admin_area_2", "year"))
  } else {
    data <- hmis_data %>%
      full_join(survey_data,     by = c("admin_area_1", "year")) %>%
      full_join(population_data, by = c("admin_area_1", "year"))
  }
  
  indicator_vars <- list(
    anc1      = c("countanc1", "anc1carry"),
    anc4      = c("countanc4", "anc4carry"),
    delivery  = c("countdelivery", "deliverycarry"),
    penta1    = c("countpenta1", "penta1carry"),
    penta2    = c("countpenta2", "penta2carry"),
    penta3    = c("countpenta3", "penta3carry"),
    opv1      = c("countopv1", "opv1carry"),
    opv2      = c("countopv2", "opv2carry"),
    opv3      = c("countopv3", "opv3carry"),
    measles1  = c("countmeasles1", "measles1carry"),
    measles2  = c("countmeasles2", "measles2carry"),
    bcg       = c("countbcg", "bcgcarry"),
    livebirth = c("countlivebirth", "livebirthcarry"),
    pnc1_mother = c("countpnc1_mother", "pnc1_mothercarry"),
    nmr       = c("countnmr", "nmrcarry")
  )
  
  available_vars <- names(data)
  
  safe_mutate <- function(var_name, formula) {
    required_vars <- indicator_vars[[var_name]]
    if (all(required_vars %in% available_vars)) formula else NA_real_
  }
  
  safe_calc <- function(expr) {
    tryCatch(expr, error = function(e) NA_real_)
  }
  
  if (all(indicator_vars$livebirth %in% available_vars)) {
    data <- data %>% mutate(
      dlivebirths_livebirth   = safe_mutate("livebirth", countlivebirth),
      dlivebirths_pregnancy   = safe_calc(dlivebirths_livebirth * (1 - 0.5 * TWIN_RATE) / ((1 - STILLBIRTH_RATE) * (1 - PREGNANCY_LOSS_RATE))),
      dlivebirths_delivery    = safe_calc(dlivebirths_pregnancy * (1 - PREGNANCY_LOSS_RATE)),
      dlivebirths_birth       = safe_calc(dlivebirths_livebirth / (1 - STILLBIRTH_RATE)),
      dlivebirths_dpt         = safe_calc(dlivebirths_livebirth * (1 - P1_NMR)),
      dlivebirths_measles1    = safe_calc(dlivebirths_dpt * (1 - P2_PNMR)),
      dlivebirths_measles2    = safe_calc(dlivebirths_dpt * (1 - 2 * P2_PNMR))
    )
  }
  
  if (all(indicator_vars$anc1 %in% available_vars)) {
    data <- data %>% mutate(
      danc1_pregnancy         = safe_mutate("anc1", countanc1 / anc1carry),
      danc1_delivery          = safe_calc(danc1_pregnancy * (1 - PREGNANCY_LOSS_RATE)),
      danc1_birth             = safe_calc(danc1_delivery / (1 - 0.5 * TWIN_RATE)),
      danc1_livebirth         = safe_calc(danc1_birth * (1 - STILLBIRTH_RATE)),
      danc1_dpt               = safe_calc(danc1_livebirth * (1 - P1_NMR)),
      danc1_measles1          = safe_calc(danc1_dpt * (1 - P2_PNMR)),
      danc1_measles2          = safe_calc(danc1_dpt * (1 - 2 * P2_PNMR))
    )
  }
  
  if (all(indicator_vars$delivery %in% available_vars)) {
    data <- data %>% mutate(
      ddelivery_livebirth     = safe_mutate("delivery", countdelivery / deliverycarry),
      ddelivery_birth         = safe_calc(ddelivery_livebirth / (1 - STILLBIRTH_RATE)),
      ddelivery_pregnancy     = safe_calc(ddelivery_birth * (1 - 0.5 * TWIN_RATE) / (1 - PREGNANCY_LOSS_RATE)),
      ddelivery_dpt           = safe_calc(ddelivery_livebirth * (1 - P1_NMR)),
      ddelivery_measles1      = safe_calc(ddelivery_dpt * (1 - P2_PNMR)),
      ddelivery_measles2      = safe_calc(ddelivery_dpt * (1 - 2 * P2_PNMR))
    )
  }
  
  if (all(indicator_vars$penta1 %in% available_vars)) {
    data <- data %>% mutate(
      dpenta1_dpt             = safe_mutate("penta1", countpenta1 / penta1carry),
      dpenta1_measles1        = safe_calc(dpenta1_dpt * (1 - P2_PNMR)),
      dpenta1_measles2        = safe_calc(dpenta1_dpt * (1 - 2 * P2_PNMR))
    )
  }
  
  if (!has_admin_area_2 && all(indicator_vars$bcg %in% available_vars)) {
    data <- data %>% mutate(
      dbcg_pregnancy = safe_mutate("bcg", (countbcg / bcgcarry) / (1 - PREGNANCY_LOSS_RATE) / (1 + TWIN_RATE) / (1 - STILLBIRTH_RATE)),
      dbcg_livebirth = safe_mutate("bcg", countbcg / bcgcarry),
      dbcg_dpt = safe_mutate("bcg", (countbcg / bcgcarry) * (1 - P1_NMR)),
      dbcg_mcv = safe_mutate("bcg", (countbcg / bcgcarry) * (1 - P1_NMR) * (1 - P2_PNMR))
    )
  }
  
  if (!has_admin_area_2) {
    data <- data %>%
      mutate(
        nummonth = if_else(is.na(nummonth) | nummonth == 0, 12, nummonth),
        
        dwpp_pregnancy = if_else(!is.na(crudebr_unwpp) & !is.na(poptot_unwpp),
                                 (crudebr_unwpp / 1000) * poptot_unwpp / (1 + TWIN_RATE), NA_real_),
        dwpp_livebirth = if_else(!is.na(crudebr_unwpp) & !is.na(poptot_unwpp),
                                 (crudebr_unwpp / 1000) * poptot_unwpp, NA_real_),
        dwpp_dpt = if_else(!is.na(totu1pop_unwpp), totu1pop_unwpp, NA_real_),
        dwpp_measles1 = if_else(!is.na(totu1pop_unwpp) & !is.na(nmrcarry),
                                totu1pop_unwpp * (1 - (nmrcarry / 100)), NA_real_),
        dwpp_measles2 = if_else(!is.na(totu1pop_unwpp) & !is.na(nmrcarry) & !is.na(postnmr),
                                totu1pop_unwpp * (1 - (nmrcarry / 100)) * (1 - (2 * postnmr / 100)), NA_real_)
      ) %>%
      mutate(
        dwpp_pregnancy = if_else(nummonth < 12, dwpp_pregnancy * (nummonth / 12), dwpp_pregnancy),
        dwpp_livebirth = if_else(nummonth < 12, dwpp_livebirth * (nummonth / 12), dwpp_livebirth),
        dwpp_dpt       = if_else(nummonth < 12, dwpp_dpt * (nummonth / 12), dwpp_dpt),
        dwpp_measles1  = if_else(nummonth < 12, dwpp_measles1 * (nummonth / 12), dwpp_measles1),
        dwpp_measles2  = if_else(nummonth < 12, dwpp_measles2 * (nummonth / 12), dwpp_measles2)
      )
  }
  
  return(data)
}

#Part 4 - calculate coverage and compare all denominators
evaluate_coverage_by_denominator <- function(data) {
  # Determine if this is national-level data
  has_admin_area_2 <- "admin_area_2" %in% names(data)
  is_national_level <- has_admin_area_2 && all(data$admin_area_2 == "NATIONAL")
  
  geo_keys <- if (has_admin_area_2) {
    c("admin_area_1", "admin_area_2", "year")
  } else {
    c("admin_area_1", "year")
  }
  
  # Numerators
  numerator_long <- data %>%
    select(all_of(geo_keys), starts_with("count")) %>%
    pivot_longer(
      cols = -all_of(geo_keys),
      names_to = "numerator_col",
      values_to = "numerator"
    ) %>%
    filter(numerator_col != "count") %>%
    mutate(indicator_common_id = str_remove(numerator_col, "^count")) %>%
    select(-numerator_col) %>%
    distinct()
  
  # Denominator pattern: match all relevant *_suffix style names
  denom_pattern <- "^(d.*)_(pregnancy|livebirth|dpt|measles1|measles2)$"
  
  # Denominator-to-indicator map (based on suffix only)
  suffix_indicator_map <- tribble(
    ~suffix,       ~indicators,
    
    # Used for indicators related to pregnancy services and ANC
    "pregnancy",   c("anc1", "anc4"),
    
    # Used for indicators that apply to newborns or children under 5
    "livebirth",   c("delivery", "bcg", "pnc1_mother"),
    
    # Used for infant immunization indicators (0–1 year)
    "dpt",         c("penta1", "penta2", "penta3", "opv1", "opv2", "opv3",
                     "pcv1", "pcv2", "pcv3", "rota1", "rota2", "ipv1", "ipv2"),
    
    # Used for coverage of first measles dose
    "measles1",    c("measles1"),
    
    # Used for coverage of second measles dose
    "measles2",    c("measles2")
  )
  
  #Denominators
  denominator_long <- data %>%
    select(all_of(geo_keys), matches(denom_pattern)) %>%
    pivot_longer(
      cols = -all_of(geo_keys),
      names_to = "denominator",
      values_to = "denominator_value"
    ) %>%
    mutate(
      denominator_type = str_extract(denominator, "(pregnancy|livebirth|dpt|measles1|measles2)"),
      indicator_common_id = purrr::map(denominator_type, ~ {
        matched <- suffix_indicator_map %>% filter(suffix == .x)
        if (nrow(matched) == 0) NA_character_ else matched$indicators[[1]]
      })
    ) %>%
    unnest_longer(indicator_common_id) %>%
    filter(!is.na(indicator_common_id)) %>%
    distinct()
  
  numerator_long <- distinct(numerator_long)
  denominator_long <- distinct(denominator_long)
  
  # Join numerator and denominator
  coverage_data <- full_join(
    numerator_long,
    denominator_long,
    by = c(geo_keys, "indicator_common_id")
  ) %>%
    mutate(coverage = numerator / denominator_value) %>%
    drop_na(coverage)
  
  # Reference values
  carry_cols <- grep("carry$", names(data), value = TRUE)
  carry_values <- data %>%
    select(all_of(geo_keys), all_of(carry_cols)) %>%
    pivot_longer(
      cols = all_of(carry_cols),
      names_to = "indicator_common_id",
      names_pattern = "(.*)carry$",
      values_to = "reference_value"
    ) %>%
    drop_na(reference_value) %>%
    group_by(across(all_of(c(geo_keys, "indicator_common_id")))) %>%
    summarise(reference_value = mean(reference_value, na.rm = TRUE), .groups = "drop")
  
  print(unique(carry_values$indicator_common_id))
  
  # Calculate error
  coverage_with_error <- left_join(
    coverage_data,
    carry_values,
    by = c(geo_keys, "indicator_common_id")
  ) %>%
    mutate(
      squared_error = (coverage - reference_value)^2,
      source_type = case_when(
        str_starts(denominator, "danc1_")      & indicator_common_id == "anc1"     ~ "reference_based",
        str_starts(denominator, "ddelivery_")  & indicator_common_id == "delivery" ~ "reference_based",
        str_starts(denominator, "dpenta1_")    & indicator_common_id == "penta1"   ~ "reference_based",
        str_starts(denominator, "dbcg_")       & indicator_common_id == "bcg"      ~ "reference_based",
        str_starts(denominator, "dwpp_")                                       ~ "unwpp_based",
        TRUE ~ "independent"
      )
    )
  
  # Rank by error
  ranked <- coverage_with_error %>%
    filter(!is.na(squared_error)) %>%
    group_by(across(all_of(geo_keys)), indicator_common_id) %>%
    arrange(squared_error) %>%
    mutate(rank = row_number()) %>%
    ungroup()
  
  # Best-only output
  best <- ranked %>%
    filter(rank == 1) %>%
    select(all_of(geo_keys), indicator_common_id,
           coverage, reference_value, denominator, denominator_type, squared_error)
  
  list(
    full_ranking = ranked,
    best_only = best
  )
}

#Part 5 - run projections
project_coverage_from_all <- function(ranked_coverage) {
  message("Projecting survey coverage forward using HMIS deltas...")
  
  if (!"reference_value" %in% names(ranked_coverage)) {
    stop("ERROR!! 'reference_value' column not found in ranked_coverage.")
  }
  
  has_admin_area_2 <- "admin_area_2" %in% names(ranked_coverage)
  geo_keys <- if (has_admin_area_2) {
    c("admin_area_1", "admin_area_2", "indicator_common_id", "denominator")
  } else {
    c("admin_area_1", "indicator_common_id", "denominator")
  }
  
  ranked_with_delta <- ranked_coverage %>%
    arrange(across(all_of(c(geo_keys, "year")))) %>%
    group_by(across(all_of(geo_keys))) %>%
    mutate(
      coverage_delta = if_else(
        !is.na(coverage) & !is.na(lag(coverage)),
        coverage - lag(coverage),
        0
      )
    ) %>%
    ungroup()
  
  all_projected <- ranked_with_delta %>%
    group_by(across(all_of(geo_keys))) %>%
    arrange(year) %>%
    mutate(
      avgsurveyprojection = first(reference_value) + cumsum(coverage_delta),
      projection_source = paste0("avgsurveyprojection_", denominator)
    ) %>%
    ungroup()
  
  return(all_projected)
}

#Part 6 - prepare outputs
prepare_combined_coverage_from_projected <- function(projected_data, raw_survey_wide) {
  has_admin_area_2 <- "admin_area_2" %in% names(projected_data)
  
  join_keys <- if (has_admin_area_2) {
    c("admin_area_1", "admin_area_2", "year", "indicator_common_id")
  } else {
    c("admin_area_1", "year", "indicator_common_id")
  }
  
  raw_survey_long <- raw_survey_wide %>%
    pivot_longer(
      cols = starts_with("rawsurvey_"),
      names_to = "indicator_common_id",
      names_prefix = "rawsurvey_",
      values_to = "coverage_original_estimate"
    ) %>%
    filter(!is.na(coverage_original_estimate)) %>%
    select(all_of(join_keys), coverage_original_estimate) %>%
    distinct()
  
  min_years <- raw_survey_long %>%
    filter(!is.na(year)) %>%
    group_by(across(setdiff(join_keys, "year"))) %>%
    summarise(min_year = min(year), .groups = "drop") %>%
    filter(!is.na(min_year) & is.finite(min_year))
  
  max_year <- max(projected_data$year, na.rm = TRUE)
  
  # Updated valid suffix-to-indicator map
  valid_suffix_map <- list(
    pregnancy  = c("anc1", "anc4"),
    livebirth  = c("bcg", "delivery", "pnc1_mother"),
    dpt        = c("penta1", "penta2", "penta3", "opv1", "opv2", "opv3",
                   "pcv1", "pcv2", "pcv3", "rota1", "rota2", "ipv1", "ipv2"),
    measles1   = c("measles1"),
    measles2   = c("measles2")
  )
  
  # Filter projected_data to only keep valid denominator-indicator pairs
  valid_denominator_map <- projected_data %>%
    select(
      admin_area_1,
      admin_area_2 = if (has_admin_area_2) "admin_area_2" else NULL,
      indicator_common_id,
      denominator
    ) %>%
    distinct() %>%
    mutate(
      suffix = str_extract(denominator, "(pregnancy|livebirth|dpt|measles1|measles2)")
    ) %>%
    filter(map2_lgl(indicator_common_id, suffix, ~ .x %in% valid_suffix_map[[.y]])) %>%
    select(-suffix)
  
  expansion_grid <- min_years %>%
    inner_join(valid_denominator_map, by = setdiff(join_keys, "year")) %>%
    rowwise() %>%
    mutate(year = list(seq.int(min_year, max_year))) %>%
    unnest(year) %>%
    ungroup() %>%
    select(-min_year)
  
  survey_expanded <- left_join(
    expansion_grid,
    raw_survey_long,
    by = join_keys
  )
  
  combined <- full_join(
    projected_data,
    survey_expanded,
    by = c(join_keys, "denominator")
  )
  
  is_national <- all(is.na(combined$admin_area_2)) || all(combined$admin_area_2 == "NATIONAL")
  
  combined <- combined %>%
    mutate(
      coverage_original_estimate = ifelse(is.nan(coverage_original_estimate), NA_real_, coverage_original_estimate),
      admin_area_2 = if (!has_admin_area_2) "NATIONAL" else admin_area_2
    )
  
  if (is_national) {
    if ("coverage_original_estimate" %in% names(combined)) {
      combined <- combined %>%
        group_by(across(all_of(c(setdiff(join_keys, "year"), "denominator")))) %>%
        mutate(
          last_survey_year = if (all(is.na(coverage_original_estimate))) NA_integer_ else max(year[!is.na(coverage_original_estimate)], na.rm = TRUE),
          avgsurveyprojection = case_when(
            year == last_survey_year ~ coverage_original_estimate,
            year > last_survey_year ~ avgsurveyprojection,
            TRUE ~ NA_real_
          ),
          coverage_original_estimate = ifelse(
            year > last_survey_year,
            NA_real_,
            coverage_original_estimate
          )
        ) %>%
        ungroup() %>%
        select(-last_survey_year)
    }
  }
  
  combined <- combined %>%
    transmute(
      admin_area_1,
      admin_area_2,
      year,
      indicator_common_id,
      denominator,
      coverage_original_estimate,
      coverage_avgsurveyprojection = avgsurveyprojection,
      coverage_cov = coverage,
      rank,
      source_type
    )
  
  combined <- combined %>%
    select(
      admin_area_1,
      admin_area_2,
      indicator_common_id,
      year,
      denominator,
      everything()
    )
  
  return(combined)
}

# ------------------------------ Main Execution ---------------------------------------------------------------
# 1 - prepare the hmis data
hmis_processed <- process_hmis_adjusted_volume(adjusted_volume_data)

# 2 - prepare the survey data
survey_processed_national <- process_survey_data(
  survey_data = survey_data_national,
  hmis_countries = hmis_processed$hmis_countries
)

national_population_processed <- process_national_population_data(
  population_data = population_estimates_only,
  hmis_countries = hmis_processed$hmis_countries
)

# 3 - calculate the denominators
denominators_national <- calculate_denominators(
  hmis_data = hmis_processed$annual_hmis,
  survey_data = survey_processed_national$carried,
  population_data = national_population_processed
)

# 4 - calculate coverage and compare the denominators
national_coverage_eval <- evaluate_coverage_by_denominator(denominators_national)

# 5 - project survey coverage forward using HMIS deltas
national_coverage_projected <- project_coverage_from_all(national_coverage_eval$full_ranking)

# 6 - prepare results and save
combined_national <- prepare_combined_coverage_from_projected(
  projected_data = national_coverage_projected,
  raw_survey_wide = survey_processed_national$raw
)

# 7 - detect the best single denominator per indicator
best_denom_per_indicator <- national_coverage_eval$full_ranking %>%
  filter(source_type == "independent") %>%
  group_by(admin_area_1, indicator_common_id, denominator) %>%
  summarise(total_error = sum(squared_error, na.rm = TRUE), .groups = "drop") %>%
  group_by(admin_area_1, indicator_common_id) %>%
  slice_min(order_by = total_error, n = 1) %>%
  ungroup()

# Clean print to console
message("Selected denominator per indicator:")
best_denom_per_indicator %>%
  arrange(indicator_common_id) %>%
  distinct(indicator_common_id, denominator) %>%
  mutate(msg = sprintf("  - %s → %s", indicator_common_id, denominator)) %>%
  pull(msg) %>%
  walk(message)

main_export <- combined_national %>%
  inner_join(
    best_denom_per_indicator,
    by = c("admin_area_1", "indicator_common_id", "denominator")
  ) %>%
  select(admin_area_1, indicator_common_id, year, denominator,
         coverage_original_estimate, coverage_avgsurveyprojection, coverage_cov)

early_survey <- combined_national %>%
  filter(is.na(coverage_cov) & !is.na(coverage_original_estimate)) %>%
  select(admin_area_1, indicator_common_id, year, coverage_original_estimate) %>%
  distinct() %>%
  mutate(
    denominator = NA_character_,
    coverage_avgsurveyprojection = NA_real_,
    coverage_cov = NA_real_
  ) %>%
  select(indicator_common_id, year,
         coverage_original_estimate, coverage_avgsurveyprojection, coverage_cov)

combined_national_export <- bind_rows(
  main_export %>%
    mutate(coverage_cov = if_else(abs(coverage_cov) < 1e-8, NA_real_, coverage_cov)) %>%
    select(indicator_common_id,
           year,
           coverage_original_estimate, 
           coverage_avgsurveyprojection, 
           coverage_cov),
  early_survey %>%
    mutate(coverage_cov = if_else(abs(coverage_cov) < 1e-8, NA_real_, coverage_cov))
)

combined_national_export_fixed <- combined_national_export %>%
  arrange(indicator_common_id, year) %>%
  group_by(indicator_common_id, year) %>%
  summarise(
    coverage_original_estimate = first(coverage_original_estimate),
    coverage_cov = first(coverage_cov),
    .groups = "drop"
  ) %>%
  group_by(indicator_common_id) %>%
  group_modify(~ {
    df <- .x
    df <- df %>% arrange(year)
    
    # Anchor year: latest year with a non-missing original estimate
    anchor_year <- max(df$year[!is.na(df$coverage_original_estimate)], na.rm = TRUE)
    anchor_idx <- which(df$year == anchor_year)[1]
    
    # Fill forward coverage_cov from anchor year to next available non-NA
    if (is.na(df$coverage_cov[anchor_idx])) {
      next_cov_idx <- which(!is.na(df$coverage_cov) & df$year > anchor_year)
      if (length(next_cov_idx) > 0) {
        fill_value <- df$coverage_cov[next_cov_idx[1]]
        df$coverage_cov[anchor_idx:next_cov_idx[1]] <- fill_value
      }
    }
    
    # Recalculate delta
    df <- df %>%
      mutate(
        cov_lag = lag(coverage_cov),
        delta = coverage_cov - cov_lag
      )
    
    # Initialize projection
    df$avgsurveyprojection <- df$coverage_original_estimate
    
    if (is.na(df$avgsurveyprojection[anchor_idx]) && !is.na(df$coverage_cov[anchor_idx])) {
      df$avgsurveyprojection[anchor_idx] <- df$coverage_cov[anchor_idx]
    }
    
    # Project forward
    if (!is.na(df$avgsurveyprojection[anchor_idx])) {
      for (i in (anchor_idx + 1):nrow(df)) {
        prev <- i - 1
        if (!is.na(df$avgsurveyprojection[prev]) && !is.na(df$delta[i])) {
          df$avgsurveyprojection[i] <- df$avgsurveyprojection[prev] + df$delta[i]
        }
      }
    }
    
    return(df)
  }) %>%
  ungroup() %>%
  select(
    indicator_common_id,
    year,
    coverage_original_estimate,
    coverage_avgsurveyprojection = avgsurveyprojection,
    coverage_cov
  )

# Clean projections - keep only years >= last survey
combined_national_export_fixed <- combined_national_export_fixed %>%
  group_by(indicator_common_id) %>%
  mutate(
    max_survey_year = ifelse(any(!is.na(coverage_original_estimate)),
                             max(year[!is.na(coverage_original_estimate)], na.rm = TRUE),
                             -Inf),
    
    coverage_avgsurveyprojection = ifelse(
      year < max_survey_year,
      NA_real_,
      coverage_avgsurveyprojection
    )
  ) %>%
  select(-max_survey_year) %>%
  ungroup() 

best_denom_summary <- best_denom_per_indicator %>%
  distinct(indicator_common_id, denominator) %>%
  arrange(indicator_common_id)


# ------------------------------ Subnational Analysis -------------------------
# Run separate analyses for admin_area_2 and admin_area_3 to get distinct output files
if (!is.null(hmis_data_subnational) && !is.null(survey_data_subnational)) {
  
  message("\n=== RUNNING SUBNATIONAL ANALYSIS ===")
  
  # Get admin_area_1 value for consistency
  admin_area_1_value <- adjusted_volume_data %>% distinct(admin_area_1) %>% pull(admin_area_1)
  hmis_data_subnational <- hmis_data_subnational %>% mutate(admin_area_1 = admin_area_1_value)
  
  # Initialize export variables
  combined_admin2_export <- NULL
  combined_admin3_export <- NULL
  
  # === ADMIN_AREA_2 ANALYSIS ===
  if (ANALYSIS_LEVEL %in% c("NATIONAL_PLUS_AA2", "NATIONAL_PLUS_AA2_AA3")) {
    message("Running admin_area_2 level analysis...")
    
    # Prepare HMIS admin_area_2 data
    hmis_admin2 <- hmis_data_subnational %>% select(-admin_area_3)
    
    # Run pipeline up to coverage evaluation
    hmis_processed_admin2 <- process_hmis_adjusted_volume(hmis_admin2, SELECTED_COUNT_VARIABLE)
    survey_processed_admin2 <- process_survey_data(survey_data_subnational, hmis_processed_admin2$hmis_countries)
    denominators_admin2 <- calculate_denominators(hmis_processed_admin2$annual_hmis, survey_processed_admin2$carried)
    coverage_eval_admin2 <- evaluate_coverage_by_denominator(denominators_admin2)
    
    # SKIP PROJECTION - just use the best coverage estimates directly
    combined_admin2_export <- coverage_eval_admin2$best_only %>%
      select(admin_area_2, indicator_common_id, year, coverage) %>%
      rename(coverage_cov = coverage)
    
    message("✓ Admin_area_2 analysis completed - ", nrow(combined_admin2_export), " result rows")
  }
  
  # === ADMIN_AREA_3 ANALYSIS ===
  if (ANALYSIS_LEVEL == "NATIONAL_PLUS_AA2_AA3") {
    message("Running admin_area_3 level analysis...")
    
    # Check if admin_area_3 data is actually usable
    if ("admin_area_3" %in% names(hmis_data_subnational)) {
      # Prepare HMIS admin_area_3 data (rename admin_area_3 to admin_area_2 for pipeline)
      hmis_admin3 <- hmis_data_subnational %>% 
        filter(!is.na(admin_area_3) & admin_area_3 != "" & admin_area_3 != "ZONE") %>%
        select(-admin_area_2) %>% 
        rename(admin_area_2 = admin_area_3)
      
      if (nrow(hmis_admin3) > 0) {
        # Run pipeline up to coverage evaluation (skip projection)
        hmis_processed_admin3 <- process_hmis_adjusted_volume(hmis_admin3, SELECTED_COUNT_VARIABLE)
        survey_processed_admin3 <- process_survey_data(survey_data_subnational, hmis_processed_admin3$hmis_countries)
        denominators_admin3 <- calculate_denominators(hmis_processed_admin3$annual_hmis, survey_processed_admin3$carried)
        coverage_eval_admin3 <- evaluate_coverage_by_denominator(denominators_admin3)
        
        # SKIP PROJECTION - just use the best coverage estimates directly (rename back to admin_area_3)
        combined_admin3_export <- coverage_eval_admin3$best_only %>%
          select(admin_area_2, indicator_common_id, year, coverage) %>%
          rename(admin_area_3 = admin_area_2, coverage_cov = coverage)
        
        message("✓ Admin_area_3 analysis completed - ", nrow(combined_admin3_export), " result rows")
      } else {
        message("No usable admin_area_3 data found")
      }
    } else {
      message("No admin_area_3 column found in HMIS data")
    }
  }
  
} else {
  message("\n=== SKIPPING SUBNATIONAL ANALYSIS ===")
  combined_admin2_export <- NULL
  combined_admin3_export <- NULL
}
# ------------------------------ Write Output Files -------------------------
# Conditionally rename pnc1_mother back to pnc1 if the original data was pnc1
if (pnc1_renamed_to_mother) {
  if (exists("combined_national_export_fixed")) {
    combined_national_export_fixed$indicator_common_id[combined_national_export_fixed$indicator_common_id == "pnc1_mother"] <- "pnc1"
    message("✓ Renamed 'pnc1_mother' to 'pnc1' in national output")
  }
  
  if (exists("combined_admin2_export")) {
    combined_admin2_export$indicator_common_id[combined_admin2_export$indicator_common_id == "pnc1_mother"] <- "pnc1"
    message("✓ Renamed 'pnc1_mother' to 'pnc1' in admin_area_2 output")
  }
  
  if (exists("combined_admin3_export")) {
    combined_admin3_export$indicator_common_id[combined_admin3_export$indicator_common_id == "pnc1_mother"] <- "pnc1"
    message("✓ Renamed 'pnc1_mother' to 'pnc1' in admin_area_3 output")
  }
  
  if (exists("best_denom_summary")) {
    best_denom_summary$indicator_common_id[best_denom_summary$indicator_common_id == "pnc1_mother"] <- "pnc1"
    message("✓ Renamed 'pnc1_mother' to 'pnc1' in best_denom_summary output")
  }
}


# Write national CSV
write.csv(combined_national_export_fixed, "M4_coverage_estimation.csv", row.names = FALSE, fileEncoding = "UTF-8")
message("✓ Saved national results: M4_coverage_estimation.csv")

# Best denominator summary
write.csv(best_denom_summary, "M4_selected_denominator_per_indicator.csv", row.names = FALSE)
message("✓ Saved denominator summary: M4_selected_denominator_per_indicator.csv")

# Write admin_area_2 CSV 
if (exists("combined_admin2_export") && 
    is.data.frame(combined_admin2_export) && 
    nrow(combined_admin2_export) > 0) {
  write.csv(combined_admin2_export, "M4_coverage_estimation_admin_area_2.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved admin_area_2 results: ", nrow(combined_admin2_export), " rows")
} else {
  # Create empty file
  dummy_data_admin2 <- data.frame(admin_area_2 = character(), indicator_common_id = character(), 
                                  year = numeric(), coverage_cov = numeric())
  write.csv(dummy_data_admin2, "M4_coverage_estimation_admin_area_2.csv", row.names = FALSE)
  message("✓ No admin_area_2 results - saved empty file")
}

# Write admin_area_3 CSV
if (exists("combined_admin3_export") && 
    is.data.frame(combined_admin3_export) && 
    nrow(combined_admin3_export) > 0) {
  write.csv(combined_admin3_export, "M4_coverage_estimation_admin_area_3.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved admin_area_3 results: ", nrow(combined_admin3_export), " rows")
} else {
  # Create empty file
  dummy_data_admin3 <- data.frame(admin_area_3 = character(), indicator_common_id = character(), 
                                  year = numeric(), coverage_cov = numeric())
  write.csv(dummy_data_admin3, "M4_coverage_estimation_admin_area_3.csv", row.names = FALSE)
  message("✓ No admin_area_3 results - saved empty file")
}

message("\n✓ Analysis complete for level: ", ANALYSIS_LEVEL)
if (exists("original_level") && original_level != ANALYSIS_LEVEL) {
  message("  (Originally requested: ", original_level, ", adjusted due to data availability)")
}
message("Output files:")
message("  - M4_coverage_estimation.csv (national)")
message("  - M4_coverage_estimation_admin_area_2.csv (zone/province level)")
message("  - M4_coverage_estimation_admin_area_3.csv (district level)")  
message("  - M4_selected_denominator_per_indicator.csv (denominator summary)")