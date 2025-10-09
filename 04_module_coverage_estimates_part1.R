COUNTRY_ISO3 <- "SOM"

SELECTED_COUNT_VARIABLE <- "count_final_both"  # Options: "count_final_none", "count_final_outlier", "count_final_completeness", "count_final_both"

PREGNANCY_LOSS_RATE <- 0.03
TWIN_RATE <- 0.015
STILLBIRTH_RATE <- 0.02
P1_NMR <- 0.041     #Default = 0.03
P2_PNMR <- 0.022
INFANT_MORTALITY_RATE <- 0.063  #Default = 0.05


ANALYSIS_LEVEL <- "NATIONAL_PLUS_AA2" # Options: "NATIONAL_ONLY", "NATIONAL_PLUS_AA2", "NATIONAL_PLUS_AA2_AA3"

#-------------------------------------------------------------------------------------------------------------
# CB - R code FASTR PROJECT
# Last edit: 2025 Oct 9
# Module: COVERAGE ESTIMATES (PART1 - DENOMINATORS)
#-------------------------------------------------------------------------------------------------------------

# ------------------------------ Load Required Libraries -----------------------------------------------------
library(dplyr)
library(tidyr)
library(zoo)
library(stringr)
library(purrr)

# ------------------------------ Define File Paths -----------------------------------------------------------
# Direct read from GitHub
PROJECT_DATA_COVERAGE <- "https://raw.githubusercontent.com/FASTR-Analytics/modules/refs/heads/main/survey_data_unified.csv"
PROJECT_DATA_POPULATION <- "https://raw.githubusercontent.com/FASTR-Analytics/modules/refs/heads/main/population_estimates_only.csv"

CURRENT_YEAR <- as.numeric(format(Sys.Date(), "%Y"))  # Dynamically get current year
MIN_YEAR <- 2000  # Set a fixed minimum year for filtering

message("✓ Step 1/7: Loading input datasets...")

# Input Datasets
message("  → Loading adjusted HMIS data (national)...")
adjusted_volume_data <- read.csv("M2_adjusted_data_national.csv", fileEncoding = "UTF-8") %>%
  mutate(iso3_code = COUNTRY_ISO3)

message("  → Loading adjusted HMIS data (subnational)...")
adjusted_volume_data_subnational <- read.csv("M2_adjusted_data_admin_area.csv", fileEncoding = "UTF-8") %>%
  mutate(iso3_code = COUNTRY_ISO3)

message("  → Loading survey data from GitHub...")
survey_data_unified <- read.csv(PROJECT_DATA_COVERAGE, fileEncoding = "UTF-8")

# Filter by ISO3 code
if ("iso3_code" %in% names(survey_data_unified)) {
  survey_data_unified <- survey_data_unified %>% filter(iso3_code == COUNTRY_ISO3)
  message("    Filtered survey data for ISO3: ", COUNTRY_ISO3)

  # Check if data exists for this country
  if (nrow(survey_data_unified) == 0) {
    warning("WARNING: No survey data found for country ISO3 code '", COUNTRY_ISO3, "'. Analysis will proceed with UNWPP data only where available.")
    # Create empty survey data structure to prevent crashes
    survey_data_unified <- data.frame(
      iso3_code = character(),
      admin_area_1 = character(),
      admin_area_2 = character(),
      year = integer(),
      indicator_common_id = character(),
      survey_value = numeric(),
      source = character(),
      source_detail = character(),
      stringsAsFactors = FALSE
    )
  } else {
    message("    ✓ Found ", nrow(survey_data_unified), " survey records for ", COUNTRY_ISO3)
  }
} else {
  warning("iso3_code column not found in survey data - cannot filter by country")
}

message("  → Loading population estimates from GitHub...")
population_estimates_only <- read.csv(PROJECT_DATA_POPULATION, fileEncoding = "UTF-8")

# Filter by ISO3 code
if ("iso3_code" %in% names(population_estimates_only)) {
  population_estimates_only <- population_estimates_only %>% filter(iso3_code == COUNTRY_ISO3)
  message("    Filtered population data for ISO3: ", COUNTRY_ISO3)

  # Check if data exists for this country
  if (nrow(population_estimates_only) == 0) {
    stop("ERROR: No population data found for country ISO3 code '", COUNTRY_ISO3, "'. Please check the ISO3 code and data availability.")
  }
  message("    ✓ Found ", nrow(population_estimates_only), " population records for ", COUNTRY_ISO3)
} else {
  warning("iso3_code column not found in population data - cannot filter by country")
}

message("✓ Step 1/7 completed: All datasets loaded successfully!")
message("================================================================================")

# ------------------------------ Prepare Data for Analysis ---------------------------------------------------

message("\n✓ Step 2/7: Preparing data for analysis...")

# A flag to track if pnc1 was renamed to pnc1_mother
pnc1_renamed_to_mother <- FALSE

# Extract country name from HMIS data (used for joins and display)
message("  → Extracting country name from HMIS data...")
COUNTRY_NAME <- unique(adjusted_volume_data$admin_area_1)

if (length(COUNTRY_NAME) > 1) {
  warning("More than one country detected in adjusted_volume_data. Using the first one found: ", COUNTRY_NAME[1])
  COUNTRY_NAME <- COUNTRY_NAME[1]
}

message("Analyzing data for country: ", COUNTRY_NAME, " (ISO3: ", COUNTRY_ISO3, ")")

message("Analysis mode: ", ANALYSIS_LEVEL)

# Always run national analysis
survey_data_national <- survey_data_unified %>% filter(admin_area_2 == "NATIONAL")

# Check if we have any national survey data
if (nrow(survey_data_national) == 0) {
  warning("No national-level survey data available for ", COUNTRY_ISO3, ". Analysis will use UNWPP population data only.")
}

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
    message("SAFEGUARD: No subnational survey data found. Falling back to: NATIONAL_ONLY")
    original_level <- ANALYSIS_LEVEL
    ANALYSIS_LEVEL <- "NATIONAL_ONLY"
    hmis_data_subnational <- NULL
    survey_data_subnational <- NULL
  } else {
    message("  ✓ Found ", nrow(survey_data_subnational), " subnational survey records")
    
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
        message("SAFEGUARD: admin_area_3 data not usable. Falling back to: NATIONAL_PLUS_AA2")
        original_level <- ANALYSIS_LEVEL
        ANALYSIS_LEVEL <- "NATIONAL_PLUS_AA2"
      }
    }
    
    # Prepare HMIS data based on final analysis level
    hmis_data_subnational <- adjusted_volume_data_subnational
  }
}

message("Final analysis level: ", ANALYSIS_LEVEL)

message("✓ Step 2/7 completed: Data preparation finished!")
message("================================================================================")

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
  
  nummonth_data <- adjusted_volume %>%
    distinct(across(all_of(c(group_vars, "month")))) %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(nummonth = n_distinct(month, na.rm = TRUE), .groups = "drop")

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
  
  # Extract ISO3 code if available
  hmis_iso3 <- if ("iso3_code" %in% names(adjusted_volume_data)) {
    unique(adjusted_volume_data$iso3_code)
  } else {
    NULL
  }

  list(
    annual_hmis = annual_hmis,
    hmis_countries = hmis_countries,
    hmis_iso3 = hmis_iso3
  )
}

# Part 2 - prepare survey data (DHS-preferred, preserves source_detail, robust to missing columns)
process_survey_data <- function(survey_data, hmis_countries, hmis_iso3 = NULL,
                                min_year = MIN_YEAR, max_year = CURRENT_YEAR) {

  # --- Harmonize, scope, coerce ---
  # For national data, filter by ISO3 if available; otherwise use admin_area_1
  if (!is.null(hmis_iso3) && "iso3_code" %in% names(survey_data)) {
    survey_data <- survey_data %>%
      filter(iso3_code %in% hmis_iso3) %>%
      mutate(
        source = tolower(source),
        year   = as.integer(year)
      )
  } else {
    survey_data <- survey_data %>%
      filter(admin_area_1 %in% hmis_countries) %>%
      mutate(
        source = tolower(source),
        year   = as.integer(year)
      )
  }

  # Intelligent PNC1 matching based on what exists in survey data
  survey_pnc_indicators <- unique(survey_data$indicator_common_id)
  has_survey_pnc1 <- "pnc1" %in% survey_pnc_indicators
  has_survey_pnc1_mother <- "pnc1_mother" %in% survey_pnc_indicators

  # Apply recoding logic
  survey_data <- survey_data %>%
    mutate(
      indicator_common_id = recode(
        indicator_common_id,
        "polio1" = "opv1", "polio2" = "opv2", "polio3" = "opv3",
        .default = indicator_common_id
      )
    )

  # Handle PNC1 mapping intelligently
  if (has_survey_pnc1 && has_survey_pnc1_mother) {
    # Both exist: keep as-is (pnc1 stays pnc1, pnc1_mother stays pnc1_mother)
    message("Survey data has both pnc1 and pnc1_mother - keeping both")
  } else if (has_survey_pnc1 && !has_survey_pnc1_mother) {
    # Only pnc1 exists: decide based on HMIS data
    if (isTRUE(pnc1_renamed_to_mother)) {
      # HMIS originally had pnc1 (renamed to pnc1_mother) -> keep survey pnc1 as pnc1
      message("Survey has pnc1, HMIS originally had pnc1 - keeping survey pnc1 as pnc1")
    } else {
      # HMIS has pnc1_mother -> rename survey pnc1 to pnc1_mother
      survey_data <- survey_data %>%
        mutate(indicator_common_id = recode(indicator_common_id, "pnc1" = "pnc1_mother"))
      message("Survey has pnc1, HMIS has pnc1_mother - renaming survey pnc1 to pnc1_mother")
    }
  } else if (!has_survey_pnc1 && has_survey_pnc1_mother) {
    # Only pnc1_mother exists: decide based on HMIS data
    if (isTRUE(pnc1_renamed_to_mother)) {
      # HMIS originally had pnc1 -> rename survey pnc1_mother to pnc1
      survey_data <- survey_data %>%
        mutate(indicator_common_id = recode(indicator_common_id, "pnc1_mother" = "pnc1"))
      message("Survey has pnc1_mother, HMIS originally had pnc1 - renaming survey pnc1_mother to pnc1")
    } else {
      # HMIS has pnc1_mother -> keep survey pnc1_mother as pnc1_mother
      message("Survey has pnc1_mother, HMIS has pnc1_mother - keeping survey pnc1_mother as pnc1_mother")
    }
  }
  
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
process_national_population_data <- function(population_data, hmis_countries, hmis_iso3 = NULL,
                                             min_year = MIN_YEAR, max_year = CURRENT_YEAR) {

  # For national data, filter by ISO3 if available; otherwise use admin_area_1
  if (!is.null(hmis_iso3) && "iso3_code" %in% names(population_data)) {
    base <- population_data %>%
      filter(admin_area_2 == "NATIONAL",
                    iso3_code %in% hmis_iso3) %>%
      mutate(
        source = tolower(source),
        year   = as.integer(year)
      )
  } else {
    base <- population_data %>%
      filter(admin_area_2 == "NATIONAL",
                    admin_area_1 %in% hmis_countries) %>%
      mutate(
        source = tolower(source),
        year   = as.integer(year)
      )
  }
  
  #original wide
  wide <- base %>%
    select(admin_area_1, year, indicator_common_id, survey_value, source) %>%
    pivot_wider(
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
  
  # optional: revert pnc1_mother -> pnc1
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
    pnc1 = c("countpnc1", "pnc1carry"),
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
  if (all(indicator_vars$livebirth %in% available_vars)) {
    data <- data %>%
      mutate(
        countlivebirth = ifelse(is.na(countlivebirth), 0, countlivebirth),
        dlivebirths_livebirth = safe_mutate("livebirth", countlivebirth / livebirthcarry),
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
  denominator_cols <- names(denominators_data)[grepl("^d(livebirth|anc1|delivery|penta1|bcg|wpp)_", names(denominators_data))]
  if (length(denominator_cols) == 0) {
    warning("No denominator columns found"); return(NULL)
  }
  
  has_admin2 <- "admin_area_2" %in% names(denominators_data)
  has_admin3 <- "admin_area_3" %in% names(denominators_data)
  
  select_cols <- if (has_admin3) {
    c("admin_area_1", "admin_area_3", "year", denominator_cols)
  } else if (has_admin2) {
    c("admin_area_1", "admin_area_2", "year", denominator_cols)
  } else {
    c("admin_area_1", "year", denominator_cols)
  }
  
  summary_stats <- denominators_data %>%
    select(all_of(select_cols)) %>%
    pivot_longer(
      cols = all_of(denominator_cols),
      names_to = "denominator_type",
      values_to = "value"
    ) %>%
    filter(!is.na(value)) %>%
    arrange(year, denominator_type)

  summary_stats
}


# Add denominator_label
add_denominator_labels <- function(df, denom_col = "denominator") {
  stopifnot(is.data.frame(df), denom_col %in% names(df))
  
  df %>%
    mutate(
      .den = .data[[denom_col]],
      den_source_key = str_replace(.den, "^d([^_]+)_.*$", "\\1"),
      den_target_key = str_replace(.den, "^d[^_]+_(.*)$", "\\1"),
      den_target_key = recode(den_target_key,
                              "livebirths" = "livebirth",
                              .default = den_target_key),
      source_phrase = case_when(
        den_source_key == "anc1"       ~ "derived from HMIS data on ANC 1st visits",
        den_source_key == "delivery"   ~ "derived from HMIS data on institutional deliveries",
        den_source_key == "bcg"        ~ "derived from HMIS data on BCG doses",
        den_source_key == "penta1"     ~ "derived from HMIS data on Penta-1 doses",
        den_source_key == "wpp"        ~ "based on UN WPP estimates",
        den_source_key == "livebirths" ~ "derived from HMIS data on live births",
        TRUE ~ "from other sources"
      ),
      target_phrase = case_when(
        den_target_key == "pregnancy" ~ "Estimated number of pregnancies",
        den_target_key == "delivery"  ~ "Estimated number of deliveries",
        den_target_key == "birth"     ~ "Estimated number of total births (live + stillbirths)",
        den_target_key == "livebirth" ~ "Estimated number of live births",
        den_target_key == "dpt"       ~ "Estimated number of infants eligible for DPT1",
        den_target_key == "measles1"  ~ "Estimated number of children eligible for measles dose 1 (MCV1)",
        den_target_key == "measles2"  ~ "Estimated number of children eligible for measles dose 2 (MCV2)",
        TRUE ~ paste("Estimated population for target", den_target_key)
      ),
      denominator_label = paste0(target_phrase, " ", source_phrase, ".")
    ) %>%
    select(-.den, -den_source_key, -den_target_key, -source_phrase, -target_phrase)
}

# Part 5 - Calculate coverage estimates (from Part 2)
calculate_coverage <- function(denominators_data, numerators_data) {

  # Dynamically determine geographic keys based on available columns
  base_geo_keys <- c("admin_area_1", "year")

  # Add admin_area_2 if it exists, otherwise add a default
  if ("admin_area_2" %in% names(denominators_data)) {
    base_geo_keys <- c(base_geo_keys, "admin_area_2")
  } else {
    denominators_data <- denominators_data %>% mutate(admin_area_2 = "NATIONAL")
    numerators_data <- numerators_data %>% mutate(admin_area_2 = "NATIONAL")
    base_geo_keys <- c(base_geo_keys, "admin_area_2")
  }

  # Add admin_area_3 if it exists
  if ("admin_area_3" %in% names(denominators_data)) {
    geo_keys <- c(base_geo_keys, "admin_area_3")
  } else {
    geo_keys <- base_geo_keys
  }

  # Map denominator targets to indicators
  target_indicator_map <- tibble::tribble(
    ~den_target, ~indicators,
    "pregnancy", c("anc1", "anc4"),
    "livebirth", c("delivery", "bcg", "pnc1_mother", "pnc1"),
    "dpt",       c("penta1", "penta2", "penta3", "opv1", "opv2", "opv3",
                   "pcv1", "pcv2", "pcv3", "rota1", "rota2", "ipv1", "ipv2"),
    "measles1",  c("measles1"),
    "measles2",  c("measles2")
  )

  # Expand denominators to match indicators
  denominator_expanded <- denominators_data %>%
    left_join(target_indicator_map, by = "den_target") %>%
    unnest_longer(indicators) %>%
    filter(!is.na(indicators)) %>%
    rename(indicator_common_id = indicators,
           denominator_value = value) %>%
    select(all_of(geo_keys), denominator, den_source, den_target,
           indicator_common_id, denominator_value)

  # Join numerators with denominators and calculate coverage
  coverage_data <- numerators_data %>%
    rename(numerator = count) %>%
    left_join(denominator_expanded, by = c(geo_keys, "indicator_common_id")) %>%
    filter(!is.na(denominator_value), denominator_value > 0) %>%
    mutate(coverage = numerator / denominator_value) %>%
    filter(!is.na(coverage), !is.infinite(coverage))

  return(coverage_data)
}

# Part 6 - Compare coverage vs carried Survey (LONG format only)
compare_coverage_to_survey <- function(coverage_data, survey_expanded_df) {
  stopifnot(is.data.frame(coverage_data), is.data.frame(survey_expanded_df))
  need_long <- c("admin_area_1","year","indicator_common_id","reference_value")
  if (!all(need_long %in% names(survey_expanded_df))) {
    stop("survey_expanded_df must be LONG with columns: ",
         paste(need_long, collapse = ", "),
         " (admin_area_2 and admin_area_3 optional; defaults to 'NATIONAL').")
  }

  # Ensure admin_area_2 exists on both sides
  if (!"admin_area_2" %in% names(coverage_data))     coverage_data     <- coverage_data %>% mutate(admin_area_2 = "NATIONAL")
  if (!"admin_area_2" %in% names(survey_expanded_df)) survey_expanded_df <- survey_expanded_df %>% mutate(admin_area_2 = "NATIONAL")

  # Handle admin_area_3 - ensure both datasets have same admin level structure
  has_admin3_cov <- "admin_area_3" %in% names(coverage_data)
  has_admin3_sur <- "admin_area_3" %in% names(survey_expanded_df)

  if (has_admin3_cov && !has_admin3_sur) {
    survey_expanded_df <- survey_expanded_df %>% mutate(admin_area_3 = "NATIONAL")
  } else if (!has_admin3_cov && has_admin3_sur) {
    coverage_data <- coverage_data %>% mutate(admin_area_3 = "NATIONAL")
  }

  # Determine geo_keys based on available columns after standardization
  geo_keys <- c("admin_area_1", "admin_area_2", "year")
  if ("admin_area_3" %in% names(coverage_data) && "admin_area_3" %in% names(survey_expanded_df)) {
    geo_keys <- c(geo_keys, "admin_area_3")
  }

  # Geographic grouping keys (exclude year for consistent denominator selection)
  geo_only_keys <- geo_keys[geo_keys != "year"]

  # Types & keys
  coverage_data$year        <- as.integer(coverage_data$year)
  survey_expanded_df$year   <- as.integer(survey_expanded_df$year)

  # Classify denominator source type
  dpt_family <- c("penta1","penta2","penta3","opv1","opv2","opv3",
                  "pcv1","pcv2","pcv3","rota1","rota2","ipv1","ipv2")
  classify_source_type <- function(denominator, ind) {
    if (startsWith(denominator, "danc1_")     && ind %in% c("anc1")) return("reference_based")
    if (startsWith(denominator, "ddelivery_") && ind %in% c("delivery"))    return("reference_based")
    if (startsWith(denominator, "dpenta1_")   && ind %in% c("penta1"))       return("reference_based")
    if (startsWith(denominator, "dbcg_")      && ind %in% c("bcg"))         return("reference_based")
    if (startsWith(denominator, "dwpp_"))                                   return("unwpp_based")
    "independent"
  }

  # Join coverage with survey reference values
  coverage_with_reference <- coverage_data %>%
    left_join(
      survey_expanded_df %>%
        select(all_of(geo_keys), indicator_common_id, reference_value),
      by = c(geo_keys, "indicator_common_id")
    ) %>%
    mutate(
      squared_error = (coverage - reference_value)^2,
      source_type   = mapply(classify_source_type, denominator, indicator_common_id)
    )
    # NOTE: Don't filter out NAs here - we need all denominators for ranking

  # Step 1: Select best denominator for each geo × indicator
  # Logic: 1) Prefer independent (non-reference, non-UNWPP) with lowest error
  #        2) Fallback to reference-based only if no independent options exist
  #        3) UNWPP denominators are EXCLUDED from "best" selection
  #        4) Exception: UNWPP used as last resort if no other denominators exist

  # Try to select best denominators excluding UNWPP
  best_denominators_no_unwpp <- coverage_with_reference %>%
    filter(!is.na(squared_error)) %>%  # Only consider denominators with survey comparisons
    filter(source_type != "unwpp_based") %>%  # EXCLUDE UNWPP from best selection
    group_by(across(all_of(geo_only_keys)), indicator_common_id) %>%
    mutate(is_reference_based = source_type == "reference_based") %>%
    arrange(
      is_reference_based,    # FALSE (independent) comes first, TRUE (reference) comes last
      squared_error,         # Among each type, pick lowest error
      .by_group = TRUE
    ) %>%
    slice_head(n = 1) %>%  # Take the top one from each group
    ungroup() %>%
    select(all_of(geo_only_keys), indicator_common_id, best_denominator = denominator)

  # Fallback: if no non-UNWPP denominators exist for some indicators, include UNWPP as last resort
  missing_indicators <- coverage_with_reference %>%
    filter(!is.na(squared_error)) %>%
    distinct(across(all_of(geo_only_keys)), indicator_common_id) %>%
    anti_join(best_denominators_no_unwpp, by = c(geo_only_keys, "indicator_common_id"))

  if (nrow(missing_indicators) > 0) {
    warning("Some indicators only have UNWPP denominators available. Including UNWPP as fallback for these cases.")

    unwpp_fallbacks <- coverage_with_reference %>%
      filter(!is.na(squared_error)) %>%
      semi_join(missing_indicators, by = c(geo_only_keys, "indicator_common_id")) %>%
      filter(source_type == "unwpp_based") %>%
      group_by(across(all_of(geo_only_keys)), indicator_common_id) %>%
      arrange(squared_error, .by_group = TRUE) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      select(all_of(geo_only_keys), indicator_common_id, best_denominator = denominator)

    best_denominators <- bind_rows(best_denominators_no_unwpp, unwpp_fallbacks)
  } else {
    best_denominators <- best_denominators_no_unwpp
  }

  # Step 2: Filter coverage data to use only the best denominator for each geo × indicator
  # and apply ranking within years for the selected denominators
  final_result <- coverage_with_reference %>%
    left_join(best_denominators, by = c(geo_only_keys, "indicator_common_id")) %>%
    filter(denominator == best_denominator) %>%
    select(-best_denominator) %>%
    # Now rank within geo × indicator × year (should be rank 1 for each since we have only one denominator per geo × indicator)
    # But only rank where we have survey reference values
    group_by(across(all_of(geo_keys)), indicator_common_id) %>%
    arrange(squared_error, .by_group = TRUE) %>%
    mutate(rank = ifelse(!is.na(squared_error), row_number(), NA_integer_)) %>%
    ungroup()

  return(final_result)
}

# Part 7 - Create combined coverage and survey result table
create_combined_results_table <- function(coverage_comparison, survey_raw_df, all_coverage_data = NULL) {

  # Get geographic keys based on what's available
  geo_keys <- c("admin_area_1", "year")
  if ("admin_area_2" %in% names(coverage_comparison)) geo_keys <- c(geo_keys, "admin_area_2")
  if ("admin_area_3" %in% names(coverage_comparison)) geo_keys <- c(geo_keys, "admin_area_3")

  # Step 1: Prepare ALL coverage results with original denominator names
  # Use all_coverage_data if provided (to include ALL denominators including UNWPP)
  if (!is.null(all_coverage_data)) {
    coverage_all <- all_coverage_data %>%
      select(
        all_of(geo_keys),
        indicator_common_id,
        denominator_best_or_survey = denominator,
        denominator_label,
        value = coverage
      )
  } else {
    # Fallback to using coverage_comparison data
    coverage_all <- coverage_comparison %>%
      select(
        all_of(geo_keys),
        indicator_common_id,
        denominator_best_or_survey = denominator,
        denominator_label,
        value = coverage,
        rank
      ) %>%
      select(-rank)
  }

  # Step 2: Create separate "best" entries (duplicate the rank=1 values)
  # Only create "best" entries where we actually have a rank (i.e., where survey comparison was possible)
  coverage_best <- coverage_comparison %>%
    filter(!is.na(rank), rank == 1) %>%
    mutate(denominator_best_or_survey = "best") %>%
    # Keep original denominator_label to show which denominator was selected as best
    select(
      all_of(geo_keys),
      indicator_common_id,
      denominator_best_or_survey,
      denominator_label,  # Keep original label
      value = coverage
    )

  # Step 3: Combine all coverage results (original denominators + best)
  coverage_results <- bind_rows(coverage_all, coverage_best)

  # Step 4: Prepare survey raw results
  # Ensure survey_raw_df has the same geographic structure as coverage
  if (!"admin_area_2" %in% names(survey_raw_df) && "admin_area_2" %in% geo_keys) {
    survey_raw_df <- survey_raw_df %>% mutate(admin_area_2 = "NATIONAL")
  }
  if (!"admin_area_3" %in% names(survey_raw_df) && "admin_area_3" %in% geo_keys) {
    survey_raw_df <- survey_raw_df %>% mutate(admin_area_3 = "NATIONAL")
  }

  # Get list of indicators that have coverage estimates
  coverage_indicators <- coverage_results %>%
    distinct(indicator_common_id) %>%
    pull(indicator_common_id)

  survey_results <- survey_raw_df %>%
    filter(!is.na(survey_value)) %>%  # Only actual survey observations
    filter(indicator_common_id %in% coverage_indicators) %>%  # Only indicators with coverage estimates
    mutate(
      denominator_best_or_survey = "survey",
      denominator_label = "Survey estimate"
    ) %>%
    select(
      all_of(geo_keys),
      indicator_common_id,
      denominator_best_or_survey,
      denominator_label,
      value = survey_value
    )

  # Step 5: Combine all results
  combined_results <- bind_rows(coverage_results, survey_results) %>%
    arrange(
      admin_area_1,
      if("admin_area_2" %in% names(.)) admin_area_2 else NULL,
      if("admin_area_3" %in% names(.)) admin_area_3 else NULL,
      indicator_common_id,
      year,
      denominator_best_or_survey
    )

  return(combined_results)
}

# ---- Helpers needed EARLY (moved up from the write-out section) ----

# Consistent rename back to pnc1 if original data had pnc1
rename_back_pnc1 <- function(x) {
  if (isTRUE(pnc1_renamed_to_mother)) recode(x, "pnc1_mother" = "pnc1", .default = x) else x
}

# Results 1: Numerators (HMIS annual counts)
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

# Results 2: Denominators (from *_summary)
make_denominators_results <- function(summary_df) {
  if (is.null(summary_df) || nrow(summary_df) == 0) return(NULL)
  
  df <- summary_df
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
      denominator = denominator_type, den_source, den_target, value
    )
}

# Results 3: Survey RAW (long; DHS-preferred; keep source + detail)
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

# Results 4: Survey REFERENCE (from carried values)
make_survey_reference_long <- function(survey_expanded_df) {
  if (is.null(survey_expanded_df) || nrow(survey_expanded_df) == 0) return(NULL)
  
  if (!"admin_area_2" %in% names(survey_expanded_df)) {
    survey_expanded_df$admin_area_2 <- "NATIONAL"
  }
  
  carry_cols <- grep("carry$", names(survey_expanded_df), value = TRUE)
  if (length(carry_cols) == 0) return(NULL)
  
  survey_expanded_df |>
    select(admin_area_1, admin_area_2, year, all_of(carry_cols)) |>
    pivot_longer(
      cols          = all_of(carry_cols),
      names_to      = "indicator_common_id",
      names_pattern = "(.*)carry$",
      values_to     = "reference_value"
    ) |>
    filter(!is.na(reference_value)) |>
    arrange(admin_area_1, admin_area_2, year, indicator_common_id)
}

# Helper: ensure admin3 outputs have the right column name
normalize_admin3_for_output <- function(df) {
  if (!is.data.frame(df) || nrow(df) == 0) return(df)
  if (!"admin_area_3" %in% names(df) && "admin_area_2" %in% names(df)) {
    df <- rename(df, admin_area_3 = admin_area_2)
  }
  df
}




# ============================== EXECUTION FLOW   ==============================

message("✓ Step 3/7: Processing national data")

# --- NATIONAL PREP ---
message("  → Processing HMIS adjusted volume data...")
hmis_processed <- process_hmis_adjusted_volume(adjusted_volume_data)

message("  → Processing survey data...")
survey_processed_national <- process_survey_data(
  survey_data    = survey_data_national,
  hmis_countries = hmis_processed$hmis_countries,
  hmis_iso3      = hmis_processed$hmis_iso3
)

message("  → Processing population data...")
national_population_processed <- process_national_population_data(
  population_data = population_estimates_only,
  hmis_countries  = hmis_processed$hmis_countries,
  hmis_iso3       = hmis_processed$hmis_iso3
)

message("  → Calculating denominators...")
denominators_national <- calculate_denominators(
  hmis_data      = hmis_processed$annual_hmis,
  survey_data    = survey_processed_national$carried,
  population_data= national_population_processed$wide
)

message("  → Creating denominator summary...")
national_summary <- create_denominator_summary(denominators_national, "NATIONAL")

# --- NATIONAL RESULTS BUILDERS (must come after summary) ---
numerators_national_long <- make_numerators_long(hmis_processed$annual_hmis)

denominators_national_results <- if (exists("national_summary")) make_denominators_results(national_summary) else NULL
if (!is.null(denominators_national_results)) {
  denominators_national_results <- add_denominator_labels(denominators_national_results, "denominator")
}

survey_raw_national_long <- make_survey_raw_long(
  dhs_mics_raw_long = survey_processed_national$raw_long,
  unwpp_raw_long    = national_population_processed$raw_long
)

survey_reference_national <- if (exists("survey_processed_national") &&
                                 is.list(survey_processed_national) &&
                                 "carried" %in% names(survey_processed_national) &&
                                 is.data.frame(survey_processed_national$carried) &&
                                 nrow(survey_processed_national$carried) > 0) {
  make_survey_reference_long(survey_processed_national$carried)
} else NULL

# --- NATIONAL COVERAGE / COMPARISON / COMBINED ---
if (!is.null(denominators_national_results) &&
    !is.null(numerators_national_long) &&
    !is.null(survey_reference_national)) {

  message("  → Calculating national coverage estimates...")
  national_coverage <- calculate_coverage(denominators_national_results, numerators_national_long)
  national_coverage <- add_denominator_labels(national_coverage)

  message("  → Comparing coverage to survey data...")
  national_comparison <- compare_coverage_to_survey(national_coverage, survey_reference_national)

  message("  → Creating combined results table...")
  national_combined_results <- create_combined_results_table(
    coverage_comparison = national_comparison,
    survey_raw_df       = survey_raw_national_long,
    all_coverage_data   = national_coverage
  )
}


message("✓ Step 3/7 completed: National analysis finished!")
message("================================================================================")

# ============================ SUBNATIONAL FLOW (IF APPLICABLE) ============================

if (!is.null(hmis_data_subnational) && !is.null(survey_data_subnational)) {

  message("✓ Step 4/7: Processing subnational data")

  # Ensure admin_area_1 is consistent
  message("  → Ensuring data consistency...")
  admin_area_1_value <- adjusted_volume_data %>% distinct(admin_area_1) %>% pull(admin_area_1)
  hmis_data_subnational <- hmis_data_subnational %>% mutate(admin_area_1 = admin_area_1_value)

  # ----------------- ADMIN_AREA_2 -----------------
  if (ANALYSIS_LEVEL %in% c("NATIONAL_PLUS_AA2", "NATIONAL_PLUS_AA2_AA3")) {

    message("  → Processing admin area 2 data...")
    
    hmis_admin2 <- hmis_data_subnational %>% select(-admin_area_3)
    
    hmis_processed_admin2   <- process_hmis_adjusted_volume(hmis_admin2, SELECTED_COUNT_VARIABLE)
    survey_processed_admin2 <- process_survey_data(survey_data_subnational, hmis_processed_admin2$hmis_countries)
    
    denominators_admin2 <- calculate_denominators(
      hmis_data   = hmis_processed_admin2$annual_hmis,
      survey_data = survey_processed_admin2$carried
    )
    
    admin2_summary <- create_denominator_summary(denominators_admin2, "ADMIN2")
    
    # --- ADMIN2 RESULTS BUILDERS ---
    numerators_admin2_long <- make_numerators_long(hmis_processed_admin2$annual_hmis)
    
    denominators_admin2_results <- if (exists("admin2_summary")) make_denominators_results(admin2_summary) else NULL
    if (!is.null(denominators_admin2_results)) {
      denominators_admin2_results <- add_denominator_labels(denominators_admin2_results, "denominator")
    }
    
    survey_raw_admin2_long <- make_survey_raw_long(
      dhs_mics_raw_long = survey_processed_admin2$raw_long,
      unwpp_raw_long    = NULL
    )
    
    survey_reference_admin2 <- if (exists("survey_processed_admin2") &&
                                   is.list(survey_processed_admin2) &&
                                   "carried" %in% names(survey_processed_admin2) &&
                                   is.data.frame(survey_processed_admin2$carried) &&
                                   nrow(survey_processed_admin2$carried) > 0) {
      make_survey_reference_long(survey_processed_admin2$carried)
    } else NULL
    
    # --- ADMIN2 COVERAGE / COMPARISON / COMBINED ---
    if (!is.null(denominators_admin2_results) &&
        !is.null(numerators_admin2_long) &&
        !is.null(survey_reference_admin2)) {

      message("  → Calculating admin area 2 coverage estimates...")
      admin2_coverage <- calculate_coverage(denominators_admin2_results, numerators_admin2_long)
      admin2_coverage <- add_denominator_labels(admin2_coverage)

      message("  → Comparing admin area 2 coverage to survey data...")
      admin2_comparison <- compare_coverage_to_survey(admin2_coverage, survey_reference_admin2)

      message("  → Creating admin area 2 combined results table...")
      admin2_combined_results <- create_combined_results_table(
        coverage_comparison = admin2_comparison,
        survey_raw_df       = survey_raw_admin2_long,
        all_coverage_data   = admin2_coverage
      )
    }
  }

  # ----------------- ADMIN_AREA_3 -----------------
  if (ANALYSIS_LEVEL == "NATIONAL_PLUS_AA2_AA3" && "admin_area_3" %in% names(hmis_data_subnational)) {

    message("  → Processing admin area 3 data...")

    hmis_admin3 <- hmis_data_subnational %>%
      filter(!is.na(admin_area_3) & admin_area_3 != "" & admin_area_3 != "ZONE") %>%
      select(-admin_area_2) %>%
      rename(admin_area_2 = admin_area_3)
    
    if (nrow(hmis_admin3) > 0) {
      hmis_processed_admin3   <- process_hmis_adjusted_volume(hmis_admin3, SELECTED_COUNT_VARIABLE)
      survey_processed_admin3 <- process_survey_data(survey_data_subnational, hmis_processed_admin3$hmis_countries)
      
      denominators_admin3 <- calculate_denominators(
        hmis_data   = hmis_processed_admin3$annual_hmis,
        survey_data = survey_processed_admin3$carried
      )
      
      admin3_summary <- create_denominator_summary(denominators_admin3, "ADMIN3")
      
      # --- ADMIN3 RESULTS BUILDERS ---
      numerators_admin3_long <- make_numerators_long(hmis_processed_admin3$annual_hmis)
      
      denominators_admin3_results <- if (exists("admin3_summary")) make_denominators_results(admin3_summary) else NULL
      if (!is.null(denominators_admin3_results)) {
        denominators_admin3_results <- add_denominator_labels(denominators_admin3_results, "denominator")
      }
      
      survey_raw_admin3_long <- make_survey_raw_long(
        dhs_mics_raw_long = survey_processed_admin3$raw_long,
        unwpp_raw_long    = NULL
      )
      survey_raw_admin3_long <- normalize_admin3_for_output(survey_raw_admin3_long)
      
      survey_reference_admin3 <- if (exists("survey_processed_admin3") &&
                                     is.list(survey_processed_admin3) &&
                                     "carried" %in% names(survey_processed_admin3) &&
                                     is.data.frame(survey_processed_admin3$carried) &&
                                     nrow(survey_processed_admin3$carried) > 0) {
        make_survey_reference_long(survey_processed_admin3$carried) %>% normalize_admin3_for_output()
      } else NULL
      
      # --- ADMIN3 COVERAGE / COMPARISON / COMBINED ---
      if (!is.null(denominators_admin3_results) &&
          !is.null(numerators_admin3_long) &&
          !is.null(survey_reference_admin3)) {

        message("  → Calculating admin area 3 coverage estimates...")
        admin3_coverage <- calculate_coverage(denominators_admin3_results, numerators_admin3_long)
        admin3_coverage <- add_denominator_labels(admin3_coverage)

        message("  → Comparing admin area 3 coverage to survey data...")
        admin3_comparison <- compare_coverage_to_survey(admin3_coverage, survey_reference_admin3)

        message("  → Creating admin area 3 combined results table...")
        admin3_combined_results <- create_combined_results_table(
          coverage_comparison = admin3_comparison,
          survey_raw_df       = survey_raw_admin3_long,
          all_coverage_data   = admin3_coverage
        )
      }
    }
  }

  
  message("✓ Step 4/7 completed: Subnational analysis finished!")

} else {
  message("✓ Step 4/7 completed: No subnational analysis (national only)!")
}
message("================================================================================")

message("✓ Step 5/7: Data processing completed! Beginning output generation...")
message("================================================================================")

message("✓ Step 6/7: Saving coverage analysis results")

message("  → Saving denominators results...")
# National
if (exists("denominators_national_results") && is.data.frame(denominators_national_results) && nrow(denominators_national_results) > 0) {
  denominators_national_results <- add_denominator_labels(denominators_national_results, "denominator")
  # Remove admin_area_2 for national results and remove label columns
  denominators_national_results <- denominators_national_results %>%
    select(-admin_area_2, -denominator_label, -den_source, -den_target)
  write.csv(denominators_national_results, "M4_denominators_national.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved denominators_national: ", nrow(denominators_national_results), " rows")
} else {
  dummy <- data.frame(
    admin_area_1      = character(),
    year              = integer(),
    denominator       = character(),
    value             = double()
  )
  write.csv(dummy, "M4_denominators_national.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No denominators_national results - saved empty file")
}

# Admin2
if (exists("denominators_admin2_results") && is.data.frame(denominators_admin2_results) && nrow(denominators_admin2_results) > 0) {
  denominators_admin2_results <- add_denominator_labels(denominators_admin2_results, "denominator")
  # Remove label columns
  denominators_admin2_results <- denominators_admin2_results %>%
    select(-denominator_label, -den_source, -den_target)
  write.csv(denominators_admin2_results, "M4_denominators_admin2.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved denominators_admin2: ", nrow(denominators_admin2_results), " rows")
} else {
  dummy <- data.frame(
    admin_area_1      = character(),
    admin_area_2      = character(),
    year              = integer(),
    denominator       = character(),
    value             = double()
  )
  write.csv(dummy, "M4_denominators_admin2.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No denominators_admin2 results - saved empty file")
}

# Admin3
if (exists("denominators_admin3_results") &&
    is.data.frame(denominators_admin3_results) &&
    nrow(denominators_admin3_results) > 0) {

  df <- normalize_admin3_for_output(denominators_admin3_results) %>%
    add_denominator_labels("denominator") %>%
    # Remove label columns
    select(admin_area_1, admin_area_3, year, denominator, value)

  write.csv(df, "M4_denominators_admin3.csv",
            row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved denominators_admin3: ", nrow(df), " rows")

} else {
  dummy <- data.frame(
    admin_area_1      = character(),
    admin_area_3      = character(),
    year              = integer(),
    denominator       = character(),
    value             = double()
  )
  write.csv(dummy, "M4_denominators_admin3.csv",
            row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No denominators_admin3 results - saved empty file")
}

# Combined Coverage and Survey Results

message("  → Saving combined coverage and survey results...")
# National
if (exists("national_combined_results") && is.data.frame(national_combined_results) && nrow(national_combined_results) > 0) {
  write.csv(national_combined_results, "M4_combined_results_national.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved combined_results_national: ", nrow(national_combined_results), " rows")
} else {
  dummy <- data.frame(
    admin_area_1 = character(),
    year = integer(),
    indicator_common_id = character(),
    denominator_best_or_survey = character(),
    denominator_label = character(),
    value = double()
  )
  write.csv(dummy, "M4_combined_results_national.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No combined_results_national results - saved empty file")
}

# Admin2
if (exists("admin2_combined_results") && is.data.frame(admin2_combined_results) && nrow(admin2_combined_results) > 0) {
  write.csv(admin2_combined_results, "M4_combined_results_admin2.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved combined_results_admin2: ", nrow(admin2_combined_results), " rows")
} else {
  dummy <- data.frame(
    admin_area_1 = character(),
    admin_area_2 = character(),
    year = integer(),
    indicator_common_id = character(),
    denominator_best_or_survey = character(),
    denominator_label = character(),
    value = double()
  )
  write.csv(dummy, "M4_combined_results_admin2.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No combined_results_admin2 results - saved empty file")
}

# Admin3
if (exists("admin3_combined_results") && is.data.frame(admin3_combined_results) && nrow(admin3_combined_results) > 0) {
  # Ensure proper column naming for admin3
  if ("admin_area_2" %in% names(admin3_combined_results)) {
    admin3_combined_results <- admin3_combined_results %>%
      rename(admin_area_3 = admin_area_2)
  }
  write.csv(admin3_combined_results, "M4_combined_results_admin3.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved combined_results_admin3: ", nrow(admin3_combined_results), " rows")
} else {
  dummy <- data.frame(
    admin_area_1 = character(),
    admin_area_3 = character(),
    year = integer(),
    indicator_common_id = character(),
    denominator_best_or_survey = character(),
    denominator_label = character(),
    value = double()
  )
  write.csv(dummy, "M4_combined_results_admin3.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No combined_results_admin3 results - saved empty file")
}


message("✓ Step 6/7 completed: All results saved successfully!")
message("================================================================================")

message("✓ Step 7/7: COVERAGE ESTIMATION ANALYSIS COMPLETE!")
message("================================================================================")
message("All output files have been generated and saved to the working directory.")
message("Check the M4_*.csv files for your coverage analysis results.")
message("================================================================================")
