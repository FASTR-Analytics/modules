SELECTED_COUNT_VARIABLE <- "count_final_both"  # Options: "count_final_none", "count_final_outlier", "count_final_completeness", "count_final_both"

CURRENT_YEAR <- as.numeric(format(Sys.Date(), "%Y"))  # Dynamically get current year
MIN_YEAR <- 2000  # Set a fixed minimum year for filtering

PREGNANCY_LOSS_RATE <- 0.03
TWIN_RATE <- 0.015
STILLBIRTH_RATE <- 0.02
P1_NMR <- 0.041     #Default = 0.03
P2_PNMR <- 0.022
INFANT_MORTALITY_RATE <- 0.063  #Default = 0.05

PROJECT_DATA_COVERAGE <-"survey_data_unified.csv"
PROJECT_DATA_POPULATION <- "population_estimates_only.csv"

ANALYSIS_LEVEL <- "NATIONAL_PLUS_AA2" # Options: "NATIONAL_ONLY", "NATIONAL_PLUS_AA2", "NATIONAL_PLUS_AA2_AA3"

#-------------------------------------------------------------------------------------------------------------
# CB - R code FASTR PROJECT
# Last edit: 2025 Sep 8
# Module: COVERAGE ESTIMATES (PART1 - DENOMINATORS)
#-------------------------------------------------------------------------------------------------------------

# ------------------------------ Load Required Libraries -----------------------------------------------------
library(dplyr)
library(tidyr)
library(zoo)
library(stringr)
library(purrr)

# ------------------------------ Define File Paths -----------------------------------------------------------
# Input Datasets
adjusted_volume_data <- read.csv("M2_adjusted_data_national.csv", fileEncoding = "UTF-8")
adjusted_volume_data_subnational <- read.csv("M2_adjusted_data_admin_area.csv", fileEncoding = "UTF-8")
survey_data_unified <- read.csv(PROJECT_DATA_COVERAGE, fileEncoding = "UTF-8")
population_estimates_only <- read.csv(PROJECT_DATA_POPULATION, fileEncoding = "UTF-8")

# ------------------------------ Prepare Data for Analysis ---------------------------------------------------
# A flag to track if pnc1 was renamed to pnc1_mother
pnc1_renamed_to_mother <- FALSE

# Dynamically detect the country from the HMIS data and filter other datasets
COUNTRY_NAME <- unique(adjusted_volume_data$admin_area_1)

if (length(COUNTRY_NAME) > 1) {
  warning("More than one country detected in adjusted_volume_data. Using the first one found: ", COUNTRY_NAME[1])
  COUNTRY_NAME <- COUNTRY_NAME[1]
}

message("Analyzing data for country: ", COUNTRY_NAME)

survey_data_unified <- survey_data_unified %>% filter(admin_area_1 == COUNTRY_NAME)
population_estimates_only <- population_estimates_only %>% filter(admin_area_1 == COUNTRY_NAME)

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

# ------------------------------ Rename for Test Instance ----------------------------------------------------
adjusted_volume_data <- adjusted_volume_data %>%
  mutate(admin_area_1 = case_when(
    admin_area_1 %in% c("Pays 001", "Country 001") ~ "Federal Govt of Somalia",
    TRUE ~ admin_area_1
  ))

# ------------------------------ Define Parameters -----------------------------------------------------------
# Coverage Estimation Parameters
coverage_params <- list(
  indicators = c(
    "anc1",
    "anc4",
    "delivery",
    "bcg",
    "penta1",
    "penta3",
    "measles1",
    "measles2",
    "rota1",
    "rota2",
    "opv1",
    "opv2",
    "opv3",
    "pnc1_mother",
    "nmr",
    "imr"
  )
)

# List of survey variables to carry forward (for forward-fill and projections)
survey_vars <- c(
  "avgsurvey_anc1",
  "avgsurvey_anc4",
  "avgsurvey_delivery",
  "avgsurvey_bcg",
  "avgsurvey_penta1",
  "avgsurvey_penta3",
  "avgsurvey_measles1",
  "avgsurvey_measles2",
  "avgsurvey_rota1",
  "avgsurvey_rota2",
  "avgsurvey_opv1",
  "avgsurvey_opv2",
  "avgsurvey_opv3",
  "avgsurvey_pnc1_mother",
  "avgsurvey_nmr",
  "avgsurvey_imr",
  "postnmr"
)

# ------------------------------ Define Functions ------------------------------------------------------------
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
  } else if (has_pnc1_mother && has_pnc1) {
    adjusted_volume_data <- adjusted_volume_data %>% filter(indicator_common_id != "pnc1")
    message("Note: Both 'pnc1_mother' and 'pnc1' found. Keeping only 'pnc1_mother'.")
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

# Part 2 - prepare survey data (DHS-preferred, preserves source_detail, robust to missing columns)
process_survey_data <- function(survey_data, hmis_countries,
                                min_year = MIN_YEAR, max_year = CURRENT_YEAR) {
  
  # --- Harmonize, scope, coerce ---
  survey_data <- survey_data %>%
    filter(admin_area_1 %in% hmis_countries) %>%
    mutate(
      indicator_common_id = recode(
        indicator_common_id,
        "polio1" = "opv1", "polio2" = "opv2", "polio3" = "opv3",
        "pnc1"   = "pnc1_mother", .default = indicator_common_id
      ),
      source = tolower(source),
      year   = as.integer(year)
    )
  
  # national vs subnational
  is_national <- all(survey_data$admin_area_2 == "NATIONAL", na.rm = TRUE)
  
  indicators <- c(
    "anc1","anc4","delivery","bcg","penta1","penta3",
    "measles1","measles2","rota1","rota2","opv1","opv2","opv3",
    "pnc1_mother","nmr","imr"
  )
  
  survey_filtered <- if (is_national) {
    survey_data %>% filter(admin_area_2 == "NATIONAL")
  } else {
    survey_data %>% filter(admin_area_2 != "NATIONAL")
  }
  
  # normalize source labels
  survey_filtered <- survey_filtered %>%
    mutate(source = case_when(
      str_detect(source, "dhs")   ~ "dhs",
      str_detect(source, "mics")  ~ "mics",
      str_detect(source, "unwpp") ~ "unwpp",
      TRUE ~ source
    ))
  
  
  # aggregate within (geo,year,indicator,source), keep a representative source_detail
  raw_by_source <- survey_filtered %>%
    filter(source %in% c("dhs","mics"),
                  year >= min_year, year <= max_year) %>%
    group_by(admin_area_1, admin_area_2, year, indicator_common_id, source) %>%
    summarise(
      survey_value   = mean(survey_value, na.rm = TRUE),
      source_detail  = first(source_detail[!is.na(source_detail)], default = NA_character_),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from  = source,
      values_from = c(survey_value, source_detail),
      names_glue  = "{.value}_{source}"
    )
  
  # ensure columns exist even if one source is entirely absent
  raw_by_source <- {
    df <- raw_by_source
    if (!"survey_value_dhs"   %in% names(df)) df$survey_value_dhs   <- NA_real_
    if (!"survey_value_mics"  %in% names(df)) df$survey_value_mics  <- NA_real_
    if (!"source_detail_dhs"  %in% names(df)) df$source_detail_dhs  <- NA_character_
    if (!"source_detail_mics" %in% names(df)) df$source_detail_mics <- NA_character_
    df
  }
  
  # pick DHS if present, else MICS
  raw_pick_long <- raw_by_source %>%
    mutate(
      survey_value  = coalesce(.data$survey_value_dhs, .data$survey_value_mics),
      source        = case_when(!is.na(.data$survey_value_dhs)  ~ "dhs",
                                       !is.na(.data$survey_value_mics) ~ "mics",
                                       TRUE ~ NA_character_),
      source_detail = case_when(source == "dhs"  ~ .data$source_detail_dhs,
                                       source == "mics" ~ .data$source_detail_mics,
                                       TRUE ~ NA_character_)
    ) %>%
    select(admin_area_1, admin_area_2, year, indicator_common_id,
                  source, source_detail, survey_value) %>%
    drop_na(survey_value)
  
  # wide values + sources + details (for convenience)
  raw_vals_wide <- raw_pick_long %>%
    select(admin_area_1, admin_area_2, year, indicator_common_id, survey_value) %>%
    pivot_wider(
      names_from  = indicator_common_id,
      values_from = survey_value,
      names_glue  = "rawsurvey_{indicator_common_id}"
    )
  
  raw_srcs_wide <- raw_pick_long %>%
    select(admin_area_1, admin_area_2, year, indicator_common_id, source) %>%
    pivot_wider(
      names_from  = indicator_common_id,
      values_from = source,
      names_glue  = "rawsource_{indicator_common_id}"
    )
  
  raw_detail_wide <- raw_pick_long %>%
    select(admin_area_1, admin_area_2, year, indicator_common_id, source_detail) %>%
    pivot_wider(
      names_from  = indicator_common_id,
      values_from = source_detail,
      names_glue  = "rawdetail_{indicator_common_id}"
    )
  
  raw_survey_values <- raw_vals_wide %>%
    left_join(raw_srcs_wide,   by = c("admin_area_1","admin_area_2","year")) %>%
    left_join(raw_detail_wide, by = c("admin_area_1","admin_area_2","year"))
  
  
  full_years <- seq(min_year, max_year)
  group_keys <- if (is_national)
    c("admin_area_1","indicator_common_id","source")
  else
    c("admin_area_1","admin_area_2","indicator_common_id","source")
  
  survey_extended <- survey_filtered %>%
    filter(year %in% full_years) %>%
    group_by(across(all_of(group_keys)), .drop = FALSE) %>%
    group_modify(~{
      if (nrow(.x) == 0) return(tibble::tibble())
      .x %>%
        complete(year = full_years) %>%
        arrange(year) %>%
        mutate(survey_value_carry = zoo::na.locf(survey_value, na.rm = FALSE))
    }) %>%
    ungroup()
  
  survey_wide <- survey_extended %>%
    select(all_of(c("admin_area_1", if (!is_national) "admin_area_2")),
                  year, indicator_common_id, source, survey_value_carry) %>%
    pivot_wider(
      names_from = c(source, indicator_common_id),
      values_from = survey_value_carry,
      names_glue  = "{indicator_common_id}_{source}",
      values_fn   = mean
    )
  
  # geo×year grid
  geo_keys <- if (is_national) c("admin_area_1") else c("admin_area_1","admin_area_2")
  geos <- (if (nrow(survey_wide)) survey_wide else survey_filtered) %>%
    distinct(across(all_of(geo_keys)))
  grid <- expand_grid(geos, year = full_years)
  
  build_last_table <- function(sf, src_label, last_name) {
    yrs <- sf %>%
      filter(source == src_label, year %in% full_years) %>%
      distinct(across(all_of(geo_keys)), year)
    grid %>%
      left_join(yrs %>% mutate(obs = 1L), by = c(geo_keys, "year")) %>%
      group_by(across(all_of(geo_keys))) %>%
      arrange(year, .by_group = TRUE) %>%
      mutate(
        tmp = if_else(!is.na(obs) & obs == 1L, year, NA_integer_),
        !!last_name := zoo::na.locf(tmp, na.rm = FALSE)
      ) %>%
      ungroup() %>%
      select(all_of(c(geo_keys, "year", last_name)))
  }
  
  dhs_last  <- build_last_table(survey_filtered, "dhs",  "dhs_lastyear")
  mics_last <- build_last_table(survey_filtered, "mics", "mics_lastyear")
  
  survey_wide <- survey_wide %>%
    right_join(grid, by = c(geo_keys, "year")) %>%
    left_join(dhs_last,  by = c(geo_keys, "year")) %>%
    left_join(mics_last, by = c(geo_keys, "year"))
  
  # choose most-recent source per year (ties → DHS)
  choose_most_recent <- function(df, ind) {
    dhs_col  <- paste0(ind, "_dhs")
    mics_col <- paste0(ind, "_mics")
    avg_col  <- paste0("avgsurvey_", ind)
    if (!(dhs_col %in% names(df)))  df[[dhs_col]]  <- NA_real_
    if (!(mics_col %in% names(df))) df[[mics_col]] <- NA_real_
    if (!("dhs_lastyear"  %in% names(df))) df$dhs_lastyear  <- NA_integer_
    if (!("mics_lastyear" %in% names(df))) df$mics_lastyear <- NA_integer_
    take_mics <- !is.na(df$mics_lastyear) &
      (is.na(df$dhs_lastyear) | df$mics_lastyear > df$dhs_lastyear)
    df[[avg_col]] <- as.numeric(ifelse(take_mics, df[[mics_col]], df[[dhs_col]]))
    df
  }
  for (ind in indicators) survey_wide <- choose_most_recent(survey_wide, ind)
  
  # postnmr
  survey_wide <- survey_wide %>%
    mutate(postnmr = ifelse("avgsurvey_imr" %in% names(.) & "avgsurvey_nmr" %in% names(.),
                                   avgsurvey_imr - avgsurvey_nmr, NA_real_))
  
  # carried panel
  carry_group <- if (is_national) "admin_area_1" else c("admin_area_1","admin_area_2")
  survey_carried <- survey_wide %>%
    group_by(across(all_of(carry_group))) %>%
    complete(year = full_seq(year, 1)) %>%
    arrange(across(all_of(carry_group)), year) %>%
    mutate(across(everything(), ~ zoo::na.locf(.x, na.rm = FALSE))) %>%
    ungroup()
  
  for (ind in c(indicators, "postnmr")) {
    avg_col   <- paste0("avgsurvey_", ind)
    carry_col <- paste0(ind, "carry")
    if (avg_col %in% names(survey_carried)) {
      survey_carried[[carry_col]] <- survey_carried[[avg_col]]
    }
  }
  
  # defaults for subnational when missing
  if (!is_national) {
    for (ind in c("bcg","penta1","penta3")) {
      carry_col <- paste0(ind, "carry")
      if (!(carry_col %in% names(survey_carried))) {
        survey_carried[[carry_col]] <- 0.70
      } else {
        survey_carried[[carry_col]] <- ifelse(is.na(survey_carried[[carry_col]]), 0.70,
                                              survey_carried[[carry_col]])
      }
    }
  }
  
  # tag NATIONAL if needed
  if (is_national) {
    survey_carried    <- survey_carried    %>% mutate(admin_area_2 = "NATIONAL")
    raw_survey_values <- raw_survey_values %>% mutate(admin_area_2 = "NATIONAL")
    raw_pick_long     <- raw_pick_long     %>% mutate(admin_area_2 = "NATIONAL")
  }
  
  list(
    carried  = survey_carried %>% arrange(across(any_of(c("admin_area_1", if (!is_national) "admin_area_2", "year")))),
    raw      = raw_survey_values,   # wide: rawsurvey_* + rawsource_* + rawdetail_* per indicator
    raw_long = raw_pick_long        # long: single row per geo–year–indicator; DHS-preferred; keeps source_detail
  )
}

# Part 2b - prepare unwpp data
process_national_population_data <- function(population_data, hmis_countries,
                                             min_year = MIN_YEAR, max_year = CURRENT_YEAR) {
  
  base <- population_data %>%
    filter(admin_area_2 == "NATIONAL",
                  admin_area_1 %in% hmis_countries) %>%
    mutate(
      source = tolower(source),
      year   = as.integer(year)
    )
  
  #original wide
  wide <- base %>%
    select(admin_area_1, year, indicator_common_id, survey_value, source) %>%
    tidyr::pivot_wider(
      names_from  = c(indicator_common_id, source),
      values_from = survey_value,
      names_glue  = "{indicator_common_id}_{source}",
      values_fn   = mean
    )
  
  #raw
  raw_long <- base %>%
    filter(source == "unwpp",
                  between(year, min_year, max_year)) %>%
    select(admin_area_1, year, indicator_common_id,
                  survey_value, source, any_of("source_detail")) %>%
    mutate(admin_area_2 = "NATIONAL") %>%
    relocate(admin_area_2, .after = admin_area_1)
  
  # ensure source_detail exists if absent
  if (!"source_detail" %in% names(raw_long)) raw_long$source_detail <- NA_character_
  
  # optional: revert pnc1_mother -> pnc1 if you use that convention elsewhere
  if (exists("rename_back_pnc1")) {
    raw_long$indicator_common_id <- rename_back_pnc1(raw_long$indicator_common_id)
  }
  
  raw_long <- raw_long %>%
    select(admin_area_1, admin_area_2, year,
                  indicator_common_id, source, source_detail, survey_value) %>%
    arrange(admin_area_1, admin_area_2, year, indicator_common_id)
  
  list(wide = wide, raw_long = raw_long)
}

#Part 3 - calculate denominators
calculate_denominators <- function(hmis_data, survey_data, population_data = NULL) {
  # nmrcarry is handled by the survey data processing, no need for redundant check here
  
  has_admin_area_2 <- "admin_area_2" %in% names(hmis_data)
  
  if (has_admin_area_2) {
    # Standard join admin_area_2 to admin_area_2 (survey data is restructured)
    data <- hmis_data %>%
      full_join(survey_data, by = c("admin_area_1", "admin_area_2", "year"))
  } else {
    data <- hmis_data %>%
      full_join(survey_data, by = c("admin_area_1", "year")) %>%
      full_join(population_data, by = c("admin_area_1", "year"))
  }
  
  indicator_vars <- list(
    anc1 = c("countanc1", "anc1carry"),
    anc4 = c("countanc4", "anc4carry"),
    delivery = c("countdelivery", "deliverycarry"),
    penta1 = c("countpenta1", "penta1carry"),
    penta2 = c("countpenta2", "penta2carry"),
    penta3 = c("countpenta3", "penta3carry"),
    opv1 = c("countopv1", "opv1carry"),
    opv2 = c("countopv2", "opv2carry"),
    opv3 = c("countopv3", "opv3carry"),
    measles1 = c("countmeasles1", "measles1carry"),
    measles2 = c("countmeasles2", "measles2carry"),
    bcg = c("countbcg", "bcgcarry"),
    livebirth = c("countlivebirth", "livebirthcarry"),
    pnc1_mother = c("countpnc1_mother", "pnc1_mothercarry"),
    nmr = c("countnmr", "nmrcarry")
  )
  
  available_vars <- names(data)
  
  safe_mutate <- function(var_name, formula) {
    required_vars <- indicator_vars[[var_name]]
    if (all(required_vars %in% available_vars)) formula else NA_real_
  }
  
  safe_calc <- function(expr) {
    tryCatch(expr, error = function(e) NA_real_)
  }
  
  # DENOMINATORS FROM LIVE BIRTH DATA
  if ("countlivebirth" %in% available_vars) {
    data <- data %>%
      mutate(
        countlivebirth = ifelse(is.na(countlivebirth), 0, countlivebirth),
        dlivebirths_livebirth = ifelse(countlivebirth > 0, countlivebirth / livebirthcarry, 0),
        dlivebirths_pregnancy = safe_calc(dlivebirths_livebirth * (1 - 0.5 * TWIN_RATE) / ((1 - STILLBIRTH_RATE) * (1 - PREGNANCY_LOSS_RATE))),
        dlivebirths_delivery = safe_calc(dlivebirths_pregnancy * (1 - PREGNANCY_LOSS_RATE)),
        dlivebirths_birth = safe_calc(dlivebirths_livebirth / (1 - STILLBIRTH_RATE)),
        dlivebirths_dpt = safe_calc(dlivebirths_livebirth * (1 - P1_NMR)),
        dlivebirths_measles1 = safe_calc(dlivebirths_dpt * (1 - P2_PNMR)),
        dlivebirths_measles2 = safe_calc(dlivebirths_dpt * (1 - 2 * P2_PNMR))
      )
  }
  
  # DENOMINATORS FROM ANC1 DATA
  if (all(indicator_vars$anc1 %in% available_vars)) {
    data <- data %>% mutate(
      danc1_pregnancy = safe_mutate("anc1", countanc1 / anc1carry),
      danc1_delivery = safe_calc(danc1_pregnancy * (1 - PREGNANCY_LOSS_RATE)),
      danc1_birth = safe_calc(danc1_delivery / (1 - 0.5 * TWIN_RATE)),
      danc1_livebirth = safe_calc(danc1_birth * (1 - STILLBIRTH_RATE)),
      danc1_dpt = safe_calc(danc1_livebirth * (1 - P1_NMR)),
      danc1_measles1 = safe_calc(danc1_dpt * (1 - P2_PNMR)),
      danc1_measles2 = safe_calc(danc1_dpt * (1 - 2 * P2_PNMR))
    )
  }
  
  # DENOMINATORS FROM DELIVERY DATA
  if (all(indicator_vars$delivery %in% available_vars)) {
    data <- data %>% mutate(
      ddelivery_livebirth = safe_mutate("delivery", countdelivery / deliverycarry),
      ddelivery_birth = safe_calc(ddelivery_livebirth / (1 - STILLBIRTH_RATE)),
      ddelivery_pregnancy = safe_calc(ddelivery_birth * (1 - 0.5 * TWIN_RATE) / (1 - PREGNANCY_LOSS_RATE)),
      ddelivery_dpt = safe_calc(ddelivery_livebirth * (1 - P1_NMR)),
      ddelivery_measles1 = safe_calc(ddelivery_dpt * (1 - P2_PNMR)),
      ddelivery_measles2 = safe_calc(ddelivery_dpt * (1 - 2 * P2_PNMR))
    )
  }
  
  # DENOMINATORS FROM PENTA1 DATA
  if (all(indicator_vars$penta1 %in% available_vars)) {
    data <- data %>% mutate(
      dpenta1_dpt = safe_mutate("penta1", countpenta1 / penta1carry),
      dpenta1_measles1 = safe_calc(dpenta1_dpt * (1 - P2_PNMR)),
      dpenta1_measles2 = safe_calc(dpenta1_dpt * (1 - 2 * P2_PNMR))
    )
  }
  
  # DENOMINATORS FROM BCG DATA (NATIONAL ANALYSIS ONLY)
  if (!has_admin_area_2 && all(indicator_vars$bcg %in% available_vars)) {
    data <- data %>% mutate(
      dbcg_pregnancy = safe_mutate("bcg", (countbcg / bcgcarry) / (1 - PREGNANCY_LOSS_RATE) / (1 + TWIN_RATE) / (1 - STILLBIRTH_RATE)),
      dbcg_livebirth = safe_mutate("bcg", countbcg / bcgcarry),
      dbcg_dpt = safe_mutate("bcg", (countbcg / bcgcarry) * (1 - P1_NMR))
    )
  }
  
  # DENOMINATORS FROM UNWPP DATA (NATIONAL ANALYSIS ONLY)
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
        dwpp_dpt = if_else(nummonth < 12, dwpp_dpt * (nummonth / 12), dwpp_dpt),
        dwpp_measles1 = if_else(nummonth < 12, dwpp_measles1 * (nummonth / 12), dwpp_measles1),
        dwpp_measles2 = if_else(nummonth < 12, dwpp_measles2 * (nummonth / 12), dwpp_measles2)
      )
  }
  
  return(data)
}

#Part 4 - prepare summary results
create_denominator_summary <- function(denominators_data, analysis_type = "NATIONAL") {
  
  # Get actual denominator columns (calculated denominators start with d + source type)
  denominator_cols <- names(denominators_data)[grepl("^d(livebirth|anc1|delivery|penta1|bcg|wpp)_", names(denominators_data))]
  
  if (length(denominator_cols) == 0) {
    warning("No denominator columns found")
    return(NULL)
  }
  
  # Determine which admin columns to include
  has_admin2 <- "admin_area_2" %in% names(denominators_data)
  has_admin3 <- "admin_area_3" %in% names(denominators_data)
  
  # Select appropriate columns
  if (has_admin3) {
    select_cols <- c("admin_area_1", "admin_area_3", "year", denominator_cols)
  } else if (has_admin2) {
    select_cols <- c("admin_area_1", "admin_area_2", "year", denominator_cols)
  } else {
    select_cols <- c("admin_area_1", "year", denominator_cols)
  }
  
  # Create simple table with denominator descriptions
  summary_stats <- denominators_data %>%
    select(all_of(select_cols)) %>%
    pivot_longer(
      cols = all_of(denominator_cols),
      names_to = "denominator_type", 
      values_to = "value"
    ) %>%
    filter(!is.na(value)) %>%
    mutate(
      target_population = case_when(
        str_ends(denominator_type, "_livebirth") ~ "Live births",
        str_ends(denominator_type, "_pregnancy") ~ "Pregnancies", 
        str_ends(denominator_type, "_delivery") ~ "Deliveries",
        str_ends(denominator_type, "_birth") ~ "Total births (live + stillbirths)",
        str_ends(denominator_type, "_dpt") ~ "Infants eligible for DPT vaccine",
        str_ends(denominator_type, "_measles1") ~ "Children eligible for measles vaccine dose 1",
        str_ends(denominator_type, "_measles2") ~ "Children eligible for measles vaccine dose 2",
        str_ends(denominator_type, "_mcv") ~ "Children eligible for measles vaccine",
        TRUE ~ "Other population"
      )
    ) %>%
    arrange(year, denominator_type)
  
  cat("\n=== ", analysis_type, " DENOMINATORS ===\n")
  print(summary_stats)
  
  return(summary_stats)
}

# ------------------------------ Main Execution --------------------------------------------------------------
# 1 - prepare the hmis data
hmis_processed <- process_hmis_adjusted_volume(adjusted_volume_data)

# 2 - prepare the survey data
survey_processed_national <- process_survey_data(
  survey_data = survey_data_national,
  hmis_countries = hmis_processed$hmis_countries
)
# 2 - prepare the survey data
national_population_processed <- process_national_population_data(
  population_data = population_estimates_only,
  hmis_countries = hmis_processed$hmis_countries
)

# 3 - calculate the denominators
denominators_national <- calculate_denominators(
  hmis_data = hmis_processed$annual_hmis,
  survey_data = survey_processed_national$carried,
  population_data = national_population_processed$wide
)

# 4 - summary 
national_summary <- create_denominator_summary(denominators_national, "NATIONAL")
# ------------------------------ Subnational Analysis --------------------------------------------------------
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
    
    # Prepare HMIS admin_area_2 data (drop admin_area_3)
    hmis_admin2 <- hmis_data_subnational %>% select(-admin_area_3)
    
    # Run pipeline
    hmis_processed_admin2 <- process_hmis_adjusted_volume(hmis_admin2, SELECTED_COUNT_VARIABLE)
    survey_processed_admin2 <- process_survey_data(survey_data_subnational, hmis_processed_admin2$hmis_countries)
    denominators_admin2 <- calculate_denominators(hmis_processed_admin2$annual_hmis, survey_processed_admin2$carried)
    admin2_summary <- create_denominator_summary(denominators_admin2, "ADMIN2")

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
      
      # Prepare survey data for admin_area_3 analysis
      survey_admin3 <- survey_data_subnational %>%
        select(-admin_area_2) %>%
        rename(admin_area_2 = admin_area_3)
      
      if (nrow(hmis_admin3) > 0) {
        # Run pipeline
        hmis_processed_admin3 <- process_hmis_adjusted_volume(hmis_admin3, SELECTED_COUNT_VARIABLE)
        survey_processed_admin3 <- process_survey_data(survey_admin3, hmis_processed_admin3$hmis_countries)
        denominators_admin3 <- calculate_denominators(hmis_processed_admin3$annual_hmis, survey_processed_admin3$carried)
        admin3_summary <- create_denominator_summary(denominators_admin3, "ADMIN3")
      }
    }
  }
}



# ------------------------------ Write Output Files ----------------------------------------------------------
# Helper: consistent rename back to pnc1 iff original data had pnc1
rename_back_pnc1 <- function(x) {
  if (isTRUE(pnc1_renamed_to_mother)) recode(x, "pnc1_mother" = "pnc1", .default = x) else x
}

# ---------- Results 1: Numerators (HMIS annual counts) ------------------------------------------------------
make_numerators_long <- function(annual_hmis_df) {
  if (is.null(annual_hmis_df) || nrow(annual_hmis_df) == 0) return(NULL)
  annual_hmis_df %>%
    { if (!"admin_area_2" %in% names(.)) mutate(., admin_area_2 = "NATIONAL") else . } %>%
    pivot_longer(
      cols = starts_with("count"),
      names_to = "indicator_common_id",
      values_to = "count",
      names_pattern = "^count(.*)$"
    ) %>%
    mutate(indicator_common_id = rename_back_pnc1(indicator_common_id)) %>%
    select(admin_area_1, admin_area_2, year, indicator_common_id, count) %>%
    arrange(admin_area_1, admin_area_2, year, indicator_common_id)
}

numerators_national_long <- make_numerators_long(hmis_processed$annual_hmis)
numerators_admin2_long   <- if (exists("hmis_processed_admin2")) make_numerators_long(hmis_processed_admin2$annual_hmis) else NULL
numerators_admin3_long   <- if (exists("hmis_processed_admin3")) make_numerators_long(hmis_processed_admin3$annual_hmis) else NULL


# ---------- Results 2: Denominators (from *_summary) --------------------------------------------------------
make_denominators_results <- function(summary_df) {
  if (is.null(summary_df) || nrow(summary_df) == 0) return(NULL)
  
  df <- summary_df
  
  # Normalize admin columns for a consistent output schema
  if ("admin_area_3" %in% names(df) && !("admin_area_2" %in% names(df))) {
    df <- rename(df, admin_area_2 = admin_area_3)
  }
  if (!"admin_area_2" %in% names(df)) {
    df <- mutate(df, admin_area_2 = "NATIONAL")
  }
  
  df %>%
    mutate(
      den_source = str_replace(denominator_type, "^d([^_]+).*", "\\1"),
      den_target = str_replace(denominator_type, ".*_([^_]+)$", "\\1")
    ) %>%
    select(
      admin_area_1, admin_area_2, year,
      denominator = denominator_type, den_source, den_target,
      target_population, value
    )
}

# Build results (only if those summary objects exist)
denominators_national_results <- if (exists("national_summary")) make_denominators_results(national_summary) else NULL
denominators_admin2_results   <- if (exists("admin2_summary"))   make_denominators_results(admin2_summary)   else NULL
denominators_admin3_results   <- if (exists("admin3_summary"))   make_denominators_results(admin3_summary)   else NULL


# ---------- Results 3: Survey RAW (long; DHS-preferred; keep source + detail) -------------------------------
# Replace your function with this two-argument version
make_survey_raw_long <- function(dhs_mics_raw_long, unwpp_raw_long = NULL) {
  
  norm <- function(df) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
    if (!"admin_area_2"  %in% names(df)) df$admin_area_2  <- "NATIONAL"
    if (!"source_detail" %in% names(df)) df$source_detail <- NA_character_
    if (!"source"        %in% names(df)) df$source        <- NA_character_
    if (!"survey_value"  %in% names(df)) df$survey_value  <- NA_real_
    
    df %>%
      mutate(
        admin_area_2 = if_else(is.na(admin_area_2) | admin_area_2 == "", "NATIONAL", admin_area_2),
        year = as.integer(year),
        indicator_common_id = if (exists("rename_back_pnc1"))
          rename_back_pnc1(indicator_common_id) else indicator_common_id
      ) %>%
      select(admin_area_1, admin_area_2, year,
                    indicator_common_id, source, source_detail, survey_value) %>%
      distinct()
  }
  
  parts <- list(norm(dhs_mics_raw_long), norm(unwpp_raw_long))
  parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, parts)
  if (length(parts) == 0) return(NULL)
  
  bind_rows(parts) %>%
    arrange(admin_area_1, admin_area_2, year, indicator_common_id)
}

survey_raw_national_long <- make_survey_raw_long(
  dhs_mics_raw_long = survey_processed_national$raw_long,
  unwpp_raw_long    = national_population_processed$raw_long
)

survey_raw_admin2_long <- if (exists("survey_processed_admin2"))
  make_survey_raw_long(survey_processed_admin2$raw_long) else NULL

survey_raw_admin3_long <- if (exists("survey_processed_admin3"))
  make_survey_raw_long(survey_processed_admin3$raw_long) else NULL



# ==============================================================================
# ============================ SAVE CSV OUTPUTS ================================
# ==============================================================================


# ---- Results 1: Numerators ---------------------------------------------------------------------------------
# National
if (exists("numerators_national_long") && is.data.frame(numerators_national_long) && nrow(numerators_national_long) > 0) {
  write.csv(numerators_national_long, "M4_numerators_national.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved numerators_national: ", nrow(numerators_national_long), " rows")
} else {
  dummy <- data.frame(
    admin_area_1 = character(),
    admin_area_2 = character(),
    year = integer(),
    indicator_common_id = character(),
    count = double()
  )
  write.csv(dummy, "M4_numerators_national.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No numerators_national results - saved empty file")
}

# Admin2
if (exists("numerators_admin2_long") && is.data.frame(numerators_admin2_long) && nrow(numerators_admin2_long) > 0) {
  write.csv(numerators_admin2_long, "M4_numerators_admin2.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved numerators_admin2: ", nrow(numerators_admin2_long), " rows")
} else {
  dummy <- data.frame(
    admin_area_1 = character(),
    admin_area_2 = character(),
    year = integer(),
    indicator_common_id = character(),
    count = double()
  )
  write.csv(dummy, "M4_numerators_admin2.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No numerators_admin2 results - saved empty file")
}

# Admin3
if (exists("numerators_admin3_long") && is.data.frame(numerators_admin3_long) && nrow(numerators_admin3_long) > 0) {
  write.csv(numerators_admin3_long, "M4_numerators_admin3.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved numerators_admin3: ", nrow(numerators_admin3_long), " rows")
} else {
  dummy <- data.frame(
    admin_area_1 = character(),
    admin_area_2 = character(),
    year = integer(),
    indicator_common_id = character(),
    count = double()
  )
  write.csv(dummy, "M4_numerators_admin3.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No numerators_admin3 results - saved empty file")
}

# ---- Results 2: Denominators (from *_summary normalized) ---------------------------------------------------
# National
if (exists("denominators_national_results") && is.data.frame(denominators_national_results) && nrow(denominators_national_results) > 0) {
  write.csv(denominators_national_results, "M4_denominators_national.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved denominators_national: ", nrow(denominators_national_results), " rows")
} else {
  dummy <- data.frame(
    admin_area_1 = character(),
    admin_area_2 = character(),
    year = integer(),
    denominator = character(),
    den_source = character(),
    den_target = character(),
    target_population = character(),
    value = double()
  )
  write.csv(dummy, "M4_denominators_national.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No denominators_national results - saved empty file")
}

# Admin2
if (exists("denominators_admin2_results") && is.data.frame(denominators_admin2_results) && nrow(denominators_admin2_results) > 0) {
  write.csv(denominators_admin2_results, "M4_denominators_admin2.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved denominators_admin2: ", nrow(denominators_admin2_results), " rows")
} else {
  dummy <- data.frame(
    admin_area_1 = character(),
    admin_area_2 = character(),
    year = integer(),
    denominator = character(),
    den_source = character(),
    den_target = character(),
    target_population = character(),
    value = double()
  )
  write.csv(dummy, "M4_denominators_admin2.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No denominators_admin2 results - saved empty file")
}

# Admin3
if (exists("denominators_admin3_results") && is.data.frame(denominators_admin3_results) && nrow(denominators_admin3_results) > 0) {
  write.csv(denominators_admin3_results, "M4_denominators_admin3.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved denominators_admin3: ", nrow(denominators_admin3_results), " rows")
} else {
  dummy <- data.frame(
    admin_area_1 = character(),
    admin_area_2 = character(),
    year = integer(),
    denominator = character(),
    den_source = character(),
    den_target = character(),
    target_population = character(),
    value = double()
  )
  write.csv(dummy, "M4_denominators_admin3.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No denominators_admin3 results - saved empty file")
}

# ---- Results 3: Survey RAW (DHS-preferred long with source + detail) ---------------------------------------
# National
if (exists("survey_raw_national_long") && is.data.frame(survey_raw_national_long) && nrow(survey_raw_national_long) > 0) {
  write.csv(survey_raw_national_long, "M4_survey_raw_national.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved survey_raw_national: ", nrow(survey_raw_national_long), " rows")
} else {
  dummy <- data.frame(
    admin_area_1 = character(),
    admin_area_2 = character(),
    year = integer(),
    indicator_common_id = character(),
    source = character(),
    source_detail = character(),
    survey_value = double()
  )
  write.csv(dummy, "M4_survey_raw_national.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No survey_raw_national results - saved empty file")
}

# Admin2
if (exists("survey_raw_admin2_long") && is.data.frame(survey_raw_admin2_long) && nrow(survey_raw_admin2_long) > 0) {
  write.csv(survey_raw_admin2_long, "M4_survey_raw_admin2.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved survey_raw_admin2: ", nrow(survey_raw_admin2_long), " rows")
} else {
  dummy <- data.frame(
    admin_area_1 = character(),
    admin_area_2 = character(),
    year = integer(),
    indicator_common_id = character(),
    source = character(),
    source_detail = character(),
    survey_value = double()
  )
  write.csv(dummy, "M4_survey_raw_admin2.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No survey_raw_admin2 results - saved empty file")
}

# Admin3
if (exists("survey_raw_admin3_long") && is.data.frame(survey_raw_admin3_long) && nrow(survey_raw_admin3_long) > 0) {
  write.csv(survey_raw_admin3_long, "M4_survey_raw_admin3.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved survey_raw_admin3: ", nrow(survey_raw_admin3_long), " rows")
} else {
  dummy <- data.frame(
    admin_area_1 = character(),
    admin_area_2 = character(),
    year = integer(),
    indicator_common_id = character(),
    source = character(),
    source_detail = character(),
    survey_value = double()
  )
  write.csv(dummy, "M4_survey_raw_admin3.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No survey_raw_admin3 results - saved empty file")
}

# ---- Results 4: Survey EXPANDED (carried/expanded panels) --------------------------------------------------
# National
if (exists("survey_processed_national") &&
    is.list(survey_processed_national) &&
    "carried" %in% names(survey_processed_national) &&
    is.data.frame(survey_processed_national$carried) &&
    nrow(survey_processed_national$carried) > 0) {
  
  write.csv(survey_processed_national$carried,
            "M4_survey_expanded_national.csv",
            row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved survey_expanded_national: ",
          nrow(survey_processed_national$carried), " rows")
} else {
  dummy <- data.frame(
    admin_area_1 = character(),
    admin_area_2 = character(),
    year = integer()
  )
  write.csv(dummy, "M4_survey_expanded_national.csv",
            row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No survey_expanded_national results - saved empty file")
}

# Admin2
if (exists("survey_processed_admin2") &&
    is.list(survey_processed_admin2) &&
    "carried" %in% names(survey_processed_admin2) &&
    is.data.frame(survey_processed_admin2$carried) &&
    nrow(survey_processed_admin2$carried) > 0) {
  
  write.csv(survey_processed_admin2$carried,
            "M4_survey_expanded_admin2.csv",
            row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved survey_expanded_admin2: ",
          nrow(survey_processed_admin2$carried), " rows")
} else {
  dummy <- data.frame(
    admin_area_1 = character(),
    admin_area_2 = character(),
    year = integer()
  )
  write.csv(dummy, "M4_survey_expanded_admin2.csv",
            row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No survey_expanded_admin2 results - saved empty file")
}

# Admin3
if (exists("survey_processed_admin3") &&
    is.list(survey_processed_admin3) &&
    "carried" %in% names(survey_processed_admin3) &&
    is.data.frame(survey_processed_admin3$carried) &&
    nrow(survey_processed_admin3$carried) > 0) {
  
  write.csv(survey_processed_admin3$carried,
            "M4_survey_expanded_admin3.csv",
            row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved survey_expanded_admin3: ",
          nrow(survey_processed_admin3$carried), " rows")
} else {
  dummy <- data.frame(
    admin_area_1 = character(),
    admin_area_2 = character(),
    year = integer()
  )
  write.csv(dummy, "M4_survey_expanded_admin3.csv",
            row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No survey_expanded_admin3 results - saved empty file")
}
