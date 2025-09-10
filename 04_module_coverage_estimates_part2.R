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
  
  geo_keys <- c("admin_area_1", "admin_area_2", "year")
  
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

# Part 2 - compare coverage vs carried Survey (returns ONE tibble)
compare_coverage_to_survey <- function(coverage_data, survey_expanded_df) {
  stopifnot(is.data.frame(coverage_data), is.data.frame(survey_expanded_df))
  
  # Geo keys: NATIONAL has admin_area_2 == "NATIONAL" in your exports
  geo_keys <- c("admin_area_1", "admin_area_2", "year")
  if (!"admin_area_2" %in% names(coverage_data)) {
    coverage_data <- coverage_data %>% mutate(admin_area_2 = "NATIONAL")
  }
  if (!"admin_area_2" %in% names(survey_expanded_df)) {
    survey_expanded_df <- survey_expanded_df %>% mutate(admin_area_2 = "NATIONAL")
  }
  
  # Build reference (carried) values from survey_expanded_df
  carry_cols <- grep("carry$", names(survey_expanded_df), value = TRUE)
  if (length(carry_cols) == 0) {
    warning("No *carry columns found in survey_expanded_df; returning empty comparison.")
    return(
      coverage_data %>%
        mutate(
          reference_value = NA_real_,
          squared_error   = NA_real_,
          source_type     = NA_character_,
          rank            = NA_integer_
        ) %>%
        slice(0)
    )
  }
  
  carry_values <- survey_expanded_df %>%
    select(any_of(geo_keys), all_of(carry_cols)) %>%
    tidyr::pivot_longer(
      cols = all_of(carry_cols),
      names_to   = "indicator_common_id",
      names_pattern = "(.*)carry$",
      values_to  = "reference_value"
    ) %>%
    filter(!is.na(reference_value)) %>%
    group_by(across(all_of(c(geo_keys, "indicator_common_id")))) %>%
    summarise(reference_value = mean(reference_value, na.rm = TRUE), .groups = "drop")
  
  # Classify denominator "source type"
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
  
  # Join, compute error, rank
  coverage_with_error <- coverage_data %>%
    left_join(carry_values, by = c(geo_keys, "indicator_common_id")) %>%
    mutate(
      squared_error = (coverage - reference_value)^2,
      source_type   = mapply(classify_source_type, denominator, indicator_common_id)
    )
  
  ranked <- coverage_with_error %>%
    filter(!is.na(squared_error)) %>%
    group_by(across(all_of(geo_keys)), indicator_common_id) %>%
    arrange(squared_error, .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup()
  
  ranked
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
  
  group_keys <- c("admin_area_1", "admin_area_2",
                  "indicator_common_id", "denominator")
  
  coverage_df %>%
    mutate(year = as.integer(year)) %>%
    group_by(across(all_of(group_keys))) %>%
    { if (complete_years) tidyr::complete(., year = tidyr::full_seq(year, 1)) else . } %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(delta = coverage - lag(coverage, n = lag_n)) %>%
    ungroup()
}

# Part 4 — project survey values using coverage deltas
project_survey_from_deltas <- function(deltas_df, survey_raw_long) {
  stopifnot(is.data.frame(deltas_df), is.data.frame(survey_raw_long))
  need_d <- c("admin_area_1","admin_area_2","year","indicator_common_id","denominator","coverage")
  need_s <- c("admin_area_1","admin_area_2","year","indicator_common_id","survey_value")
  stopifnot(all(need_d %in% names(deltas_df)))
  stopifnot(all(need_s %in% names(survey_raw_long)))
  
  # last observed survey per (geo, indicator)
  baseline <- survey_raw_long %>%
    group_by(admin_area_1, admin_area_2, indicator_common_id) %>%
    filter(year == max(year, na.rm = TRUE)) %>%
    slice_tail(n = 1) %>%   # tie-break safety
    ungroup() %>%
    transmute(
      admin_area_1, admin_area_2, indicator_common_id,
      baseline_year = as.integer(year),
      baseline_value = as.numeric(survey_value)
    )
  
  # compute year-on-year deltas (ensure strictly by denominator)
  deltas <- deltas_df %>%
    group_by(admin_area_1, admin_area_2, indicator_common_id, denominator) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(delta = coverage - lag(coverage)) %>%
    ungroup() %>%
    filter(!is.na(delta)) # only years with a defined lag
  
  # attach baseline to every denom path for that (geo, indicator)
  seeds <- deltas %>%
    distinct(admin_area_1, admin_area_2, indicator_common_id, denominator) %>%
    left_join(baseline, by = c("admin_area_1","admin_area_2","indicator_common_id")) %>%
    filter(!is.na(baseline_year), !is.na(baseline_value))
  
  # build projections
  proj <- deltas %>%
    inner_join(seeds,
                      by = c("admin_area_1","admin_area_2","indicator_common_id","denominator")) %>%
    group_by(admin_area_1, admin_area_2, indicator_common_id, denominator) %>%
    # ensure we start at baseline_year (copy forward baseline to first delta year)
    arrange(year, .by_group = TRUE) %>%
    mutate(
      # cumulative sum of deltas AFTER the baseline year
      cum_delta = cumsum(if_else(year > baseline_year, delta, 0)),
      projected = baseline_value + cum_delta
    ) %>%
    ungroup() %>%
    select(
      admin_area_1, admin_area_2, year, indicator_common_id, denominator,
      baseline_year, projected
    )
  
  # also include an explicit baseline row for traceability (optional)
  baseline_rows <- seeds %>%
    transmute(
      admin_area_1, admin_area_2,
      year = baseline_year,
      indicator_common_id, denominator,
      baseline_year, projected = baseline_value
    )
  
  bind_rows(proj, baseline_rows) %>%
    distinct() %>%
    arrange(admin_area_1, admin_area_2, indicator_common_id, denominator, year)
}

# Part 5 - prepare result tables
build_final_results <- function(coverage_df, proj_df, survey_raw_df = NULL) {

  # Ensure admin_area_2 exists
  if (!"admin_area_2" %in% names(coverage_df)) coverage_df <- mutate(coverage_df, admin_area_2 = "NATIONAL")
  if (!"admin_area_2" %in% names(proj_df))     proj_df     <- mutate(proj_df,     admin_area_2 = "NATIONAL")
  if (!is.null(survey_raw_df) && !"admin_area_2" %in% names(survey_raw_df)) {
    survey_raw_df <- mutate(survey_raw_df, admin_area_2 = "NATIONAL")
  }
  
  # Required cols
  need_cov  <- c("admin_area_1","admin_area_2","year","indicator_common_id","denominator","coverage")
  need_proj <- c("admin_area_1","admin_area_2","year","indicator_common_id","denominator","projected")
  stopifnot(all(need_cov  %in% names(coverage_df)))
  stopifnot(all(need_proj %in% names(proj_df)))
  
  # Add labels if missing (does not change your existing wording)
  if (!"denominator_label" %in% names(coverage_df)) {
    coverage_df <- add_denominator_labels(coverage_df)
  }
  
  # 1) HMIS coverage (coverage_cov)
  cov_base <- coverage_df %>%
    select(
      admin_area_1, admin_area_2, year,
      indicator_common_id, denominator, denominator_label,
      coverage_cov = coverage
    )
  
  # 2) Projections (coverage_avgsurveyprojection)
  cov_proj <- cov_base %>%
    left_join(
      proj_df %>%
        select(
          admin_area_1, admin_area_2, year,
          indicator_common_id, denominator,
          coverage_avgsurveyprojection = projected
        ),
      by = c("admin_area_1","admin_area_2","year","indicator_common_id","denominator")
    )
  
  # If no survey provided, return HMIS+proj only
  if (is.null(survey_raw_df)) {
    return(
      cov_proj %>%
        mutate(
          coverage_original_estimate = NA_real_,
          survey_raw_source = NA_character_,
          survey_raw_source_detail = NA_character_
        ) %>%
        distinct() %>%
        arrange(admin_area_1, admin_area_2, indicator_common_id, denominator, year)
    )
  }
  
  # 3) Collapse survey RAW to one row per geo-year-indicator
  survey_slim <- survey_raw_df %>%
    group_by(admin_area_1, admin_area_2, year, indicator_common_id) %>%
    summarise(
      coverage_original_estimate = mean(survey_value, na.rm = TRUE),
      survey_raw_source          = paste(sort(unique(stats::na.omit(source))), collapse = "; "),
      survey_raw_source_detail   = paste(sort(unique(stats::na.omit(source_detail))), collapse = "; "),
      .groups = "drop"
    )
  
  # 4) Build the denominator universe per (geo, indicator)
  denom_index <- coverage_df %>%
    distinct(admin_area_1, admin_area_2, indicator_common_id, denominator, denominator_label)
  
  # 5) Expand survey years across ALL denominators for that (geo, indicator)
  #    (this is what brings in earlier survey years, replicated per denominator)
  survey_expanded <- denom_index %>%
    inner_join(
      survey_slim,
      by = c("admin_area_1","admin_area_2","indicator_common_id")
    )
  # Note: joined *without year* initially to replicate across denoms,
  # then the year comes from survey_slim; no many-to-many warning because
  # we intend the cartesian expansion at this stage.
  
  # 6) Union HMIS+proj with survey-expanded (includes early years)
  final <- cov_proj %>%
    full_join(
      survey_expanded,
      by = c("admin_area_1","admin_area_2","year","indicator_common_id","denominator","denominator_label")
    ) %>%
    distinct() %>%
    arrange(admin_area_1, admin_area_2, indicator_common_id, denominator, year)
  
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

# ---- Required schemas ----
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

subnat_required_cols <- c(
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

# ---------------- NATIONAL (no admin_area_2) ----------------
if (exists("final_national") && is.data.frame(final_national) && nrow(final_national) > 0) {
  # drop admin_area_2 if it exists
  if ("admin_area_2" %in% names(final_national)) final_national$admin_area_2 <- NULL
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
    coverage_original_estimate_source = character(),
    coverage_original_estimate = double(),
    coverage_avgsurveyprojection = double(),
    coverage_cov = double(),
    survey_raw_value = double(),
    survey_raw_source = character(),
    survey_raw_source_detail = character(),
    stringsAsFactors = FALSE
  )
  write.csv(dummy_nat, "M5_final_national.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No NATIONAL final results - saved empty file")
}

# ---------------- ADMIN2 (keeps admin_area_2) ----------------
if (exists("final_admin2") && is.data.frame(final_admin2) && nrow(final_admin2) > 0) {
  for (cn in setdiff(subnat_required_cols, names(final_admin2))) final_admin2[[cn]] <- NA
  final_admin2 <- final_admin2[, subnat_required_cols]
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
    coverage_original_estimate_source = character(),
    coverage_original_estimate = double(),
    coverage_avgsurveyprojection = double(),
    coverage_cov = double(),
    survey_raw_value = double(),
    survey_raw_source = character(),
    survey_raw_source_detail = character(),
    stringsAsFactors = FALSE
  )
  write.csv(dummy_a2, "M5_final_admin2.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No ADMIN2 final results - saved empty file (or ADMIN2 skipped)")
}

# ---------------- ADMIN3 (keeps admin_area_2) ----------------
if (exists("final_admin3") && is.data.frame(final_admin3) && nrow(final_admin3) > 0) {
  for (cn in setdiff(subnat_required_cols, names(final_admin3))) final_admin3[[cn]] <- NA
  final_admin3 <- final_admin3[, subnat_required_cols]
  write.csv(final_admin3, "M5_final_admin3.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ Saved M5_final_admin3.csv: ", nrow(final_admin3), " rows")
} else {
  dummy_a3 <- data.frame(
    admin_area_1 = character(),
    admin_area_2 = character(),
    year = integer(),
    indicator_common_id = character(),
    denominator = character(),
    denominator_label = character(),
    coverage_original_estimate_source = character(),
    coverage_original_estimate = double(),
    coverage_avgsurveyprojection = double(),
    coverage_cov = double(),
    survey_raw_value = double(),
    survey_raw_source = character(),
    survey_raw_source_detail = character(),
    stringsAsFactors = FALSE
  )
  write.csv(dummy_a3, "M5_final_admin3.csv", row.names = FALSE, fileEncoding = "UTF-8")
  message("✓ No ADMIN3 final results - saved empty file (or ADMIN3 skipped)")
}

message("✓ All done.")
