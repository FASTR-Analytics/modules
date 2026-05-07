SELECTED_COUNT_VARIABLE <- "count_final_none"
INTERPOLATE_POPULATION <- FALSE

#-------------------------------------------------------------------------------------------------------------
# M8: Catalog-Driven Scorecard Module
#
# Computes scorecard indicators at AA3 × month grain from the calculated indicators catalog.
# Outputs numerator and denominator columns; SQL aggregates to any time period.
#
# OUTPUT:
#   M8_output_scorecard.csv — One row per admin area × period × indicator, with num/denom columns
#
#-------------------------------------------------------------------------------------------------------------

message("Starting M8 Scorecard module...")

library(dplyr)
library(readr)
library(tidyr)

# Load Data ------------------------------------------------------------------------------------------------
ADJUSTED_DATA_FILE <- "M2_adjusted_data.csv"
POPULATION_FILE <- "population.csv"

message("Loading adjusted data from M2...")
adjusted_data <- read_csv(ADJUSTED_DATA_FILE, show_col_types = FALSE) %>%
  mutate(
    year = as.integer(substr(as.character(period_id), 1, 4)),
    month = as.integer(substr(as.character(period_id), 5, 6))
  )

message("Loading population data...")
population_raw <- read_csv(POPULATION_FILE, show_col_types = FALSE)
if ("admin_area_1" %in% names(population_raw)) {
  population_raw <- population_raw %>% select(-admin_area_1)
}
# Pivot so each population_type becomes a column (total_population, u5, wra, etc.)
population <- population_raw %>%
  pivot_wider(names_from = population_type, values_from = count)

# Derive nhmis_timely_and_data: on-time reports that also contain data
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

# Detect available periods (monthly) -----------------------------------------------------------------------
available_periods <- adjusted_data %>%
  distinct(period_id, year) %>%
  arrange(period_id)

pop_years <- unique(as.integer(substr(as.character(population$period_id), 1, 4)))

# Filter to periods where population is available (same year or previous year as fallback)
available_periods <- available_periods %>%
  filter(year %in% pop_years | (year - 1) %in% pop_years)

if (nrow(available_periods) == 0) {
  stop("ERROR: No periods found with both indicator data and population data!")
}

message(sprintf("Found %d period(s) with indicator and population data", nrow(available_periods)))

# Helper functions -----------------------------------------------------------------------------------------
# Determine geo columns dynamically based on what exists in BOTH datasets
all_possible_geo_cols <- c("admin_area_2", "admin_area_3", "admin_area_4")
geo_cols_in_adjusted <- intersect(all_possible_geo_cols, names(adjusted_data))
geo_cols_in_population <- intersect(all_possible_geo_cols, names(population))

message(sprintf("Adjusted data columns: %s", paste(names(adjusted_data), collapse=", ")))
message(sprintf("Population columns: %s", paste(names(population), collapse=", ")))
message(sprintf("Geo cols in adjusted data: %s", paste(geo_cols_in_adjusted, collapse=", ")))
message(sprintf("Geo cols in population: %s", paste(geo_cols_in_population, collapse=", ")))

# Use the common geo columns (population may be at coarser level than adjusted data)
geo_cols <- intersect(geo_cols_in_adjusted, geo_cols_in_population)
if (length(geo_cols) == 0) {
  stop("ERROR: No common geographic columns between adjusted data and population data!")
}
message(sprintf("Using common geo_cols for aggregation: %s", paste(geo_cols, collapse=", ")))

aggregate_to_period <- function(data, target_period_id) {
  data %>%
    filter(period_id == target_period_id) %>%
    group_by(across(all_of(geo_cols)), indicator_common_id) %>%
    summarise(count = sum(.data[[SELECTED_COUNT_VARIABLE]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = indicator_common_id, values_from = count)
}

finest_geo_col <- geo_cols[length(geo_cols)]
message(sprintf("Finest common geo column: %s", finest_geo_col))

get_population_for_period <- function(pop_data, target_period_id, area_data) {
  area_names <- unique(area_data[[finest_geo_col]])
  message(sprintf("  Looking for %d unique areas in %s", length(area_names), finest_geo_col))
  target_year <- as.integer(substr(as.character(target_period_id), 1, 4))

  # Population data is already wide-format with columns like total_population, u5, wra, etc.
  # Find matching period (try current year, then fallback to previous year)
  current_year_periods <- sprintf("%04d%02d", target_year, 12:1)
  pop_result <- NULL

  for (period in current_year_periods) {
    pop_result <- pop_data %>%
      filter(period_id == period, .data[[finest_geo_col]] %in% area_names)
    if (nrow(pop_result) > 0) {
      break
    }
  }

  # Fallback to previous year if needed
  if (is.null(pop_result) || nrow(pop_result) == 0) {
    prev_periods <- sprintf("%04d%02d", target_year - 1, 12:1)
    for (period in prev_periods) {
      pop_result <- pop_data %>%
        filter(period_id == period, .data[[finest_geo_col]] %in% area_names)
      if (nrow(pop_result) > 0) {
        break
      }
    }
  }

  if (is.null(pop_result) || nrow(pop_result) == 0) {
    # Check if it's a geographic mismatch
    pop_areas <- unique(pop_data[[finest_geo_col]])
    matching <- intersect(area_names, pop_areas)
    if (length(matching) == 0) {
      stop(sprintf(
        "ERROR: No matching areas between indicator data and population data!\n  Indicator data %s values (sample): %s\n  Population data %s values (sample): %s\n  This usually means the datasets are from different countries or have mismatched admin area names.",
        finest_geo_col, paste(head(area_names, 3), collapse=", "),
        finest_geo_col, paste(head(pop_areas, 3), collapse=", ")
      ))
    }
    stop(sprintf("ERROR: No population data found for period %s or fallback (year %d or %d)!",
                 target_period_id, target_year, target_year - 1))
  }

  # Return only finest_geo_col + population type columns (avoid duplicate geo cols in join)
  pop_type_cols <- setdiff(names(pop_result), c(geo_cols, "period_id"))
  pop_result %>% select(all_of(c(finest_geo_col, pop_type_cols)))
}

# Build scorecard rows for one period ----------------------------------------------------------------------
build_scorecard_for_period <- function(data, target_period_id) {
  # Aggregate to finest geo level × period and pivot wide
  agg <- aggregate_to_period(data, target_period_id)

  if (nrow(agg) == 0) {
    message(sprintf("  No data for period %s, skipping", target_period_id))
    return(tibble())
  }

  # Merge population and add period_id back (lost during pivot)
  pop <- get_population_for_period(population, target_period_id, agg)
  data <- agg %>%
    left_join(pop, by = finest_geo_col) %>%
    mutate(period_id = target_period_id)

  # Apply calculated indicator blocks (generated by codegen)
  __CALCULATED_INDICATOR_BLOCKS__
}

# Main execution -------------------------------------------------------------------------------------------
message(sprintf("\nCalculating M8 Scorecard for %d period(s) using %s...",
                nrow(available_periods), SELECTED_COUNT_VARIABLE))

all_results <- list()

for (i in seq_len(nrow(available_periods))) {
  target_period <- available_periods$period_id[i]
  message(sprintf("  Processing period %s...", target_period))

  result <- build_scorecard_for_period(adjusted_data, target_period)
  if (nrow(result) > 0) {
    all_results[[length(all_results) + 1]] <- result
  }
}

if (length(all_results) == 0) {
  message("WARNING: No results produced!")
  # Create empty tibble with dynamic geo columns
  empty_cols <- setNames(
    lapply(c(geo_cols, "period_id", "indicator_common_id", "numerator", "denominator"),
           function(x) if (x %in% c("numerator", "denominator")) numeric() else if (x == "period_id") integer() else character()),
    c(geo_cols, "period_id", "indicator_common_id", "numerator", "denominator")
  )
  final_output <- as_tibble(empty_cols)
} else {
  final_output <- bind_rows(all_results) %>%
    arrange(across(all_of(c(finest_geo_col, "period_id", "indicator_common_id"))))
}

message("Writing output file...")
write_csv(final_output, "M8_output_scorecard.csv")

n_indicators <- length(unique(final_output$indicator_common_id))
n_areas <- length(unique(final_output[[finest_geo_col]]))
message(sprintf("\nM8 Scorecard complete!"))
message(sprintf("  Output: M8_output_scorecard.csv"))
message(sprintf("  Rows: %d", nrow(final_output)))
message(sprintf("  Indicators: %d", n_indicators))
message(sprintf("  Areas: %d", n_areas))
message(sprintf("  Periods: %d", length(all_results)))
