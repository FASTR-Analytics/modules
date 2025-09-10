## GLOBAL PARAMETERS GO HERE

#-------------------------------------------------------------------------------------------------------------
# CB - R code FASTR PROJECT
# Last edit: 2025 Sep 9
# Module: COVERAGE ESTIMATES (PART2 - COVERAGE)
#-------------------------------------------------------------------------------------------------------------

# ------------------------------ Load Required Libraries -----------------------------------------------------
library(dplyr)
library(tidyr)
library(zoo)
library(stringr)
library(purrr)

# ------------------------------ Define Analysis Parameters --------------------------------------------------
# These parameters control the administrative levels (national, admin2, admin3)
# for which the coverage analysis will be performed.
RUN_NATIONAL <- TRUE  # Always run national
RUN_ADMIN2 <- FALSE   # Will be set based on data
RUN_ADMIN3 <- FALSE   # Will be set based on data

#------------------------------- Load the Data ---------------------------------------------------------------
# Load denominators first to check data availability
denominators_national <- read.csv("M4_denominators_national.csv", fileEncoding = "UTF-8")
denominators_admin2 <- read.csv("M4_denominators_admin2.csv", fileEncoding = "UTF-8")
denominators_admin3 <- read.csv("M4_denominators_admin3.csv", fileEncoding = "UTF-8")

# Check which admin levels have data and update global parameters
RUN_ADMIN2 <- nrow(denominators_admin2) > 0
RUN_ADMIN3 <- nrow(denominators_admin3) > 0

# Message about data availability
if (!RUN_ADMIN2) message("No data in admin2 denominators - admin2 coverage analysis will be skipped")
if (!RUN_ADMIN3) message("No data in admin3 denominators - admin3 coverage analysis will be skipped")

# Load numerators
numerators_national <- read.csv("M4_numerators_national.csv", fileEncoding = "UTF-8")
if (RUN_ADMIN2) numerators_admin2 <- read.csv("M4_numerators_admin2.csv", fileEncoding = "UTF-8")
if (RUN_ADMIN3) numerators_admin3 <- read.csv("M4_numerators_admin3.csv", fileEncoding = "UTF-8")


# Load carried survey data
survey_expanded_national <- read.csv("M4_survey_expanded_national.csv", fileEncoding = "UTF-8")
if (RUN_ADMIN2) survey_expanded_admin2 <- read.csv("M4_survey_expanded_admin2.csv", fileEncoding = "UTF-8")
if (RUN_ADMIN3) survey_expanded_admin3 <- read.csv("M4_survey_expanded_admin3.csv", fileEncoding = "UTF-8")

# Load raw survey data
survey_raw_national <- read.csv("M4_survey_raw_national.csv", fileEncoding = "UTF-8")
if (RUN_ADMIN2) survey_raw_admin2 <- read.csv("M4_survey_raw_admin2.csv", fileEncoding = "UTF-8")
if (RUN_ADMIN3) survey_raw_admin3 <- read.csv("M4_survey_raw_admin3.csv", fileEncoding = "UTF-8")


# ------------------------------ Define Functions ------------------------------------------------------------
# Part 1 - calculate coverage
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
    "livebirth", c("delivery", "bcg", "pnc1_mother"),
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

# Part 2 - compare coverage vs carried Survey (LONG format only)
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
  
  # Types & keys
  coverage_data$year        <- as.integer(coverage_data$year)
  survey_expanded_df$year   <- as.integer(survey_expanded_df$year)
  
  # Classify denominator source type
  dpt_family <- c("penta1","penta2","penta3","opv1","opv2","opv3",
                  "pcv1","pcv2","pcv3","rota1","rota2","ipv1","ipv2")
  classify_source_type <- function(denominator, ind) {
    if (startsWith(denominator, "danc1_")     && ind %in% c("anc1","anc4")) return("reference_based")
    if (startsWith(denominator, "ddelivery_") && ind %in% c("delivery"))    return("reference_based")
    if (startsWith(denominator, "dpenta1_")   && ind %in% dpt_family)       return("reference_based")
    if (startsWith(denominator, "dbcg_")      && ind %in% c("bcg"))         return("reference_based")
    if (startsWith(denominator, "dwpp_"))                                   return("unwpp_based")
    "independent"
  }
  
  # Join, compute error, rank within geo × indicator
  coverage_data %>%
    left_join(
      survey_expanded_df %>%
        select(all_of(geo_keys), indicator_common_id, reference_value),
      by = c(geo_keys, "indicator_common_id")
    ) %>%
    mutate(
      squared_error = (coverage - reference_value)^2,
      source_type   = mapply(classify_source_type, denominator, indicator_common_id)
    ) %>%
    filter(!is.na(squared_error)) %>%
    group_by(across(all_of(geo_keys)), indicator_common_id) %>%
    arrange(squared_error, .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup()
}

# Part 3 — calculate delta per indicator × denominator × geo
coverage_deltas <- function(coverage_df,
                            lag_n = 1,
                            complete_years = TRUE) {
  stopifnot(is.data.frame(coverage_df))
  if (!"coverage" %in% names(coverage_df)) stop("`coverage_df` must contain `coverage`.")
  if (!"denominator" %in% names(coverage_df)) stop("`coverage_df` must contain `denominator`.")
  if (!"admin_area_2" %in% names(coverage_df)) {
    coverage_df <- mutate(coverage_df, admin_area_2 = "NATIONAL")
  }
  
  # Determine group keys based on available columns
  group_keys <- c("admin_area_1", "admin_area_2", "indicator_common_id", "denominator")
  if ("admin_area_3" %in% names(coverage_df)) {
    group_keys <- c(group_keys, "admin_area_3")
  }
  
  coverage_df %>%
    mutate(year = as.integer(year)) %>%
    group_by(across(all_of(group_keys))) %>%
    { if (complete_years) complete(., year = full_seq(year, 1)) else . } %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(delta = coverage - lag(coverage, n = lag_n)) %>%
    ungroup()
}

# Part 4 — project survey values using coverage deltas
project_survey_from_deltas <- function(deltas_df, survey_raw_long) {
  stopifnot(is.data.frame(deltas_df), is.data.frame(survey_raw_long))
  
  # Add admin_area_2 if missing from survey_raw_long
  if (!"admin_area_2" %in% names(survey_raw_long)) {
    survey_raw_long <- survey_raw_long %>% mutate(admin_area_2 = "NATIONAL")
  }
  
  # Determine required columns based on available admin levels
  need_d <- c("admin_area_1","admin_area_2","year","indicator_common_id","denominator","coverage")
  need_s <- c("admin_area_1","admin_area_2","year","indicator_common_id","survey_value")
  
  if ("admin_area_3" %in% names(deltas_df)) {
    need_d <- c(need_d, "admin_area_3")
  }
  if ("admin_area_3" %in% names(survey_raw_long)) {
    need_s <- c(need_s, "admin_area_3")
  }
  
  # Check required columns exist
  missing_d <- setdiff(need_d, names(deltas_df))
  missing_s <- setdiff(need_s, names(survey_raw_long))
  
  if (length(missing_d) > 0) {
    stop(paste("Missing columns in deltas_df:", paste(missing_d, collapse = ", ")))
  }
  if (length(missing_s) > 0) {
    stop(paste("Missing columns in survey_raw_long:", paste(missing_s, collapse = ", ")))
  }
  
  # Determine grouping keys
  group_keys <- c("admin_area_1", "admin_area_2", "indicator_common_id")
  if ("admin_area_3" %in% names(survey_raw_long)) {
    group_keys <- c(group_keys, "admin_area_3")
  }
  
  # last observed survey per (geo, indicator)
  baseline <- survey_raw_long %>%
    group_by(across(all_of(group_keys))) %>%
    filter(year == max(year, na.rm = TRUE)) %>%
    slice_tail(n = 1) %>%   # tie-break safety
    ungroup() %>%
    transmute(
      across(all_of(group_keys)),
      baseline_year = as.integer(year),
      baseline_value = as.numeric(survey_value)
    )
  
  # Determine delta grouping keys
  delta_group_keys <- c(group_keys, "denominator")
  
  # compute year-on-year deltas (ensure strictly by denominator)
  deltas <- deltas_df %>%
    group_by(across(all_of(delta_group_keys))) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(delta = coverage - lag(coverage)) %>%
    ungroup() %>%
    filter(!is.na(delta)) # only years with a defined lag
  
  # attach baseline to every denom path for that (geo, indicator)
  seeds <- deltas %>%
    distinct(across(all_of(delta_group_keys))) %>%
    left_join(baseline, by = group_keys) %>%
    filter(!is.na(baseline_year), !is.na(baseline_value))
  
  # build projections
  proj <- deltas %>%
    inner_join(seeds, by = delta_group_keys) %>%
    group_by(across(all_of(delta_group_keys))) %>%
    # ensure we start at baseline_year (copy forward baseline to first delta year)
    arrange(year, .by_group = TRUE) %>%
    mutate(
      # cumulative sum of deltas AFTER the baseline year
      cum_delta = cumsum(if_else(year > baseline_year, delta, 0)),
      projected = baseline_value + cum_delta
    ) %>%
    ungroup() %>%
    select(
      all_of(group_keys), year, indicator_common_id, denominator,
      baseline_year, projected
    )
  
  # also include an explicit baseline row for traceability (optional)
  baseline_rows <- seeds %>%
    transmute(
      across(all_of(group_keys)),
      year = baseline_year,
      indicator_common_id, denominator,
      baseline_year, projected = baseline_value
    )
  
  bind_rows(proj, baseline_rows) %>%
    distinct() %>%
    arrange(across(all_of(c(group_keys, "indicator_common_id", "denominator", "year"))))
}

# Part 5 - prepare result tables
build_final_results <- function(coverage_df, proj_df, survey_raw_df = NULL) {
  
  # Ensure admin_area_2 exists
  if (!"admin_area_2" %in% names(coverage_df)) coverage_df <- mutate(coverage_df, admin_area_2 = "NATIONAL")
  if (!"admin_area_2" %in% names(proj_df))     proj_df     <- mutate(proj_df,     admin_area_2 = "NATIONAL")
  if (!is.null(survey_raw_df) && !"admin_area_2" %in% names(survey_raw_df)) {
    survey_raw_df <- mutate(survey_raw_df, admin_area_2 = "NATIONAL")
  }
  
  # Required columns
  need_cov  <- c("admin_area_1","admin_area_2","year","indicator_common_id","denominator","coverage")
  need_proj <- c("admin_area_1","admin_area_2","year","indicator_common_id","denominator","projected")
  if ("admin_area_3" %in% names(coverage_df)) need_cov  <- c(need_cov, "admin_area_3")
  if ("admin_area_3" %in% names(proj_df))     need_proj <- c(need_proj, "admin_area_3")
  stopifnot(all(need_cov  %in% names(coverage_df)))
  stopifnot(all(need_proj %in% names(proj_df)))
  
  # Add labels if missing
  if (!"denominator_label" %in% names(coverage_df)) {
    coverage_df <- add_denominator_labels(coverage_df)
  }
  
  # Keys
  base_keys <- c("admin_area_1", "admin_area_2")
  if ("admin_area_3" %in% names(coverage_df)) base_keys <- c(base_keys, "admin_area_3")
  join_keys <- c(base_keys, "year", "indicator_common_id", "denominator")
  
  # 1) HMIS coverage
  cov_base <- coverage_df %>%
    select(
      all_of(base_keys), year,
      indicator_common_id, denominator, denominator_label,
      coverage_cov = coverage
    )
  
  # 2) Projections
  cov_proj <- cov_base %>%
    left_join(
      proj_df %>%
        select(
          all_of(base_keys), year,
          indicator_common_id, denominator,
          coverage_avgsurveyprojection = projected
        ),
      by = join_keys
    )
  
  # If no survey: return HMIS + projections only
  if (is.null(survey_raw_df)) {
    return(
      cov_proj %>%
        mutate(
          coverage_original_estimate = NA_real_,
          survey_raw_source = NA_character_,
          survey_raw_source_detail = NA_character_
        ) %>%
        distinct() %>%
        arrange(across(all_of(c(base_keys, "indicator_common_id", "denominator", "year"))))
    )
  }
  
  # 3) Collapse survey RAW
  survey_group_keys <- base_keys
  if ("admin_area_3" %in% names(survey_raw_df)) survey_group_keys <- c(survey_group_keys, "admin_area_3")
  survey_group_keys <- c(survey_group_keys, "year", "indicator_common_id")
  
  survey_slim <- survey_raw_df %>%
    group_by(across(all_of(survey_group_keys))) %>%
    summarise(
      coverage_original_estimate = mean(survey_value, na.rm = TRUE),
      survey_raw_source          = paste(sort(unique(stats::na.omit(source))), collapse = "; "),
      survey_raw_source_detail   = paste(sort(unique(stats::na.omit(source_detail))), collapse = "; "),
      .groups = "drop"
    )
  
  # 4) Denominator universe
  denom_index_keys <- c(base_keys, "indicator_common_id")
  denom_index <- coverage_df %>%
    distinct(across(all_of(c(denom_index_keys, "denominator", "denominator_label"))))
  
  # 5) Expand survey across ALL denominators
  survey_expanded <- denom_index %>%
    inner_join(survey_slim, by = denom_index_keys)
  
  # 6) Union HMIS+proj with survey-expanded
  final_join_keys <- c(base_keys, "year", "indicator_common_id", "denominator", "denominator_label")
  
  final <- cov_proj %>%
    full_join(survey_expanded, by = final_join_keys) %>%
    distinct() %>%
    arrange(across(all_of(c(base_keys, "indicator_common_id", "denominator", "year")))) %>%
    
    # For each (geo, indicator, denominator), find the last year with a non-NA survey value,
    # and set projections to NA for years BEFORE that last survey year.
    group_by(across(all_of(c(base_keys, "indicator_common_id", "denominator")))) %>%
    mutate(
      .last_svy_year = suppressWarnings(max(year[!is.na(coverage_original_estimate)], na.rm = TRUE)),
      .last_svy_year = ifelse(is.infinite(.last_svy_year), NA_real_, .last_svy_year),
      coverage_avgsurveyprojection = if_else(
        !is.na(.last_svy_year) & year < .last_svy_year,
        NA_real_,
        coverage_avgsurveyprojection
      )
    ) %>%
    ungroup() %>%
    select(-.last_svy_year)
  
  final
}

# ------------------------------ Main Execution ------------------------------

# ===== NATIONAL (always) =====
message("Step 1 (NATIONAL): Calculating coverage...")
coverage_national <- calculate_coverage(denominators_national, numerators_national)
coverage_national <- add_denominator_labels(coverage_national)
message("✓ Coverage complete")

message("Step 2 (NATIONAL): Comparing to carried survey values...")
comparison_national <- compare_coverage_to_survey(coverage_national, survey_expanded_national)
message("✓ Comparison complete")

message("Step 3 (NATIONAL): Computing deltas...")
coverage_delta_national <- coverage_deltas(coverage_national)
message("✓ Deltas complete")

message("Step 4 (NATIONAL): Projecting survey from deltas...")
proj_survey_national <- project_survey_from_deltas(
  deltas_df       = coverage_delta_national,
  survey_raw_long = survey_raw_national
)
message("✓ Projection complete")

message("Step 5 (NATIONAL): Preparing final results...")
final_national <- build_final_results(
  coverage_df   = coverage_national,
  proj_df       = proj_survey_national,
  survey_raw_df = survey_raw_national
)
message("✓ Final results (NATIONAL) ready")

# ===== ADMIN2 (conditional) =====
if (RUN_ADMIN2) {
  message("Step 1 (ADMIN2): Calculating coverage...")
  coverage_admin2 <- calculate_coverage(denominators_admin2, numerators_admin2)
  coverage_admin2 <- add_denominator_labels(coverage_admin2)
  message("✓ Coverage complete")
  
  message("Step 2 (ADMIN2): Comparing to carried survey values...")
  comparison_admin2 <- compare_coverage_to_survey(coverage_admin2, survey_expanded_admin2)
  message("✓ Comparison complete")
  
  message("Step 3 (ADMIN2): Computing deltas...")
  coverage_delta_admin2 <- coverage_deltas(coverage_admin2)
  message("✓ Deltas complete")
  
  message("Step 4 (ADMIN2): Projecting survey from deltas...")
  proj_survey_admin2 <- project_survey_from_deltas(
    deltas_df       = coverage_delta_admin2,
    survey_raw_long = survey_raw_admin2
  )
  message("✓ Projection complete")
  
  message("Step 5 (ADMIN2): Preparing final results...")
  final_admin2 <- build_final_results(
    coverage_df   = coverage_admin2,
    proj_df       = proj_survey_admin2,
    survey_raw_df = survey_raw_admin2
  )
  message("✓ Final results (ADMIN2) ready")
} else {
  message("Admin2 disabled or no data; skipping ADMIN2 block.")
}

# ===== ADMIN3 (conditional) =====
if (RUN_ADMIN3) {
  message("Step 1 (ADMIN3): Calculating coverage...")
  coverage_admin3 <- calculate_coverage(denominators_admin3, numerators_admin3)
  coverage_admin3 <- add_denominator_labels(coverage_admin3)
  message("✓ Coverage complete")
  
  message("Step 2 (ADMIN3): Comparing to carried survey values...")
  comparison_admin3 <- compare_coverage_to_survey(coverage_admin3, survey_expanded_admin3)
  message("✓ Comparison complete")
  
  message("Step 3 (ADMIN3): Computing deltas...")
  coverage_delta_admin3 <- coverage_deltas(coverage_admin3)
  message("✓ Deltas complete")
  
  message("Step 4 (ADMIN3): Projecting survey from deltas...")
  proj_survey_admin3 <- project_survey_from_deltas(
    deltas_df       = coverage_delta_admin3,
    survey_raw_long = survey_raw_admin3
  )
  message("✓ Projection complete")
  
  message("Step 5 (ADMIN3): Preparing final results...")
  final_admin3 <- build_final_results(
    coverage_df   = coverage_admin3,
    proj_df       = proj_survey_admin3,
    survey_raw_df = survey_raw_admin3
  )
  message("✓ Final results (ADMIN3) ready")
} else {
  message("Admin3 disabled or no data; skipping ADMIN3 block.")
}

# ==============================================================================
# ============================ SAVE CSV OUTPUTS ================================
# ==============================================================================
message("Saving CSVs...")

# ---- Required fields ----
nat_required_cols <- c(
  "admin_area_1", 
  "year", 
  "indicator_common_id",
  "denominator", 
  "denominator_label",
  "coverage_original_estimate",
  "coverage_avgsurveyprojection",
  "coverage_cov",
  "survey_raw_source",
  "survey_raw_source_detail"
)

admin2_required_cols <- c(
  "admin_area_1",
  "admin_area_2",
  "year",
  "indicator_common_id",
  "denominator",
  "denominator_label",
  "coverage_original_estimate",
  "coverage_avgsurveyprojection",
  "coverage_cov",
  "survey_raw_source",
  "survey_raw_source_detail"
)

admin3_required_cols <- c(
  "admin_area_1",
  "admin_area_3",  # Changed from admin_area_2 to admin_area_3
  "year",
  "indicator_common_id",
  "denominator",
  "denominator_label",
  "coverage_original_estimate",
  "coverage_avgsurveyprojection",
  "coverage_cov",
  "survey_raw_source",
  "survey_raw_source_detail"
)

# ---------------- NATIONAL (no admin_area_2) ----------------
if (exists("final_national") && is.data.frame(final_national) && nrow(final_national) > 0) {
  # drop admin_area_2 if it exists
  if ("admin_area_2" %in% names(final_national)) final_national$admin_area_2 <- NULL
  # drop admin_area_3 if it exists (shouldn't be in national)
  if ("admin_area_3" %in% names(final_national)) final_national$admin_area_3 <- NULL
  # add any missing cols as NA, and order
  for (cn in setdiff(nat_required_cols, names(final_national))) final_national[[cn]] <- NA
  final_national <- final_national[, nat_required_cols]
  write.csv(final_national, "M5_final_national.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved M5_final_national.csv: ", nrow(final_national), " rows")
} else {
  dummy_nat <- data.frame(
    admin_area_1 = character(),
    year = integer(),
    indicator_common_id = character(),
    denominator = character(),
    denominator_label = character(),
    coverage_original_estimate = double(),
    coverage_avgsurveyprojection = double(),
    coverage_cov = double(),
    survey_raw_source = character(),
    survey_raw_source_detail = character(),
    stringsAsFactors = FALSE
  )
  write.csv(dummy_nat, "M5_final_national.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No NATIONAL final results - saved empty file")
}

# ---------------- ADMIN2 (keeps admin_area_2) ----------------
if (exists("final_admin2") && is.data.frame(final_admin2) && nrow(final_admin2) > 0) {
  # drop admin_area_3 if it exists (shouldn't be in admin2)
  if ("admin_area_3" %in% names(final_admin2)) final_admin2$admin_area_3 <- NULL
  for (cn in setdiff(admin2_required_cols, names(final_admin2))) final_admin2[[cn]] <- NA
  final_admin2 <- final_admin2[, admin2_required_cols]
  write.csv(final_admin2, "M5_final_admin2.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved M5_final_admin2.csv: ", nrow(final_admin2), " rows")
} else {
  dummy_a2 <- data.frame(
    admin_area_1 = character(),
    admin_area_2 = character(),
    year = integer(),
    indicator_common_id = character(),
    denominator = character(),
    denominator_label = character(),
    coverage_original_estimate = double(),
    coverage_avgsurveyprojection = double(),
    coverage_cov = double(),
    survey_raw_source = character(),
    survey_raw_source_detail = character(),
    stringsAsFactors = FALSE
  )
  write.csv(dummy_a2, "M5_final_admin2.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No ADMIN2 final results - saved empty file (or ADMIN2 skipped)")
}

# ---------------- ADMIN3 (keeps admin_area_3, removes admin_area_2) ----------------
if (exists("final_admin3") && is.data.frame(final_admin3) && nrow(final_admin3) > 0) {
  # Remove admin_area_2 if it exists (it was added by our functions but shouldn't be in final output)
  if ("admin_area_2" %in% names(final_admin3)) final_admin3$admin_area_2 <- NULL
  # Add any missing columns as NA
  for (cn in setdiff(admin3_required_cols, names(final_admin3))) final_admin3[[cn]] <- NA
  # Reorder columns to match schema
  final_admin3 <- final_admin3[, admin3_required_cols]
  write.csv(final_admin3, "M5_final_admin3.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved M5_final_admin3.csv: ", nrow(final_admin3), " rows")
} else {
  dummy_a3 <- data.frame(
    admin_area_1 = character(),
    admin_area_3 = character(),  # Changed from admin_area_2 to admin_area_3
    year = integer(),
    indicator_common_id = character(),
    denominator = character(),
    denominator_label = character(),
    coverage_original_estimate = double(),
    coverage_avgsurveyprojection = double(),
    coverage_cov = double(),
    survey_raw_source = character(),
    survey_raw_source_detail = character(),
    stringsAsFactors = FALSE
  )
  write.csv(dummy_a3, "M5_final_admin3.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No ADMIN3 final results - saved empty file (or ADMIN3 skipped)")
}

message("✓ All done.")