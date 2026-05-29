SELECTED_COUNT_VARIABLE <- "count_final_none"
__SKIP_MISSING_INDICATORS_VALUE__ <- TRUE
__NEEDS_POPULATION_VALUE__ <- FALSE

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
SKIP_MISSING_INDICATORS <- __SKIP_MISSING_INDICATORS_VALUE__

message("Loading adjusted data from M2...")
adjusted_data <- read_csv(ADJUSTED_DATA_FILE, show_col_types = FALSE) %>%
  mutate(
    year = as.integer(substr(as.character(period_id), 1, 4)),
    month = as.integer(substr(as.character(period_id), 5, 6))
  )

message("Loading population data...")
population_raw <- read_csv(POPULATION_FILE, show_col_types = FALSE)
NEEDS_POPULATION <- __NEEDS_POPULATION_VALUE__
if (!NEEDS_POPULATION) {
  message("  No calculated indicators use a population denominator - population.csv will be ignored")
}
HAS_POPULATION <- NEEDS_POPULATION && nrow(population_raw) > 0

if (HAS_POPULATION) {
  if ("admin_area_1" %in% names(population_raw)) {
    population_raw <- population_raw %>% select(-admin_area_1)
  }
  if ("period_id" %in% names(population_raw) && !"year" %in% names(population_raw)) {
    population_raw <- population_raw %>%
      mutate(year = as.integer(substr(as.character(period_id), 1, 4))) %>%
      select(-period_id)
  }
  population <- population_raw %>%
    pivot_wider(names_from = population_type, values_from = count)
  message(sprintf("  Loaded %d population records", nrow(population)))
  pop_years <- sort(unique(population$year))
  message(sprintf("  Population years available: %s", paste(pop_years, collapse=", ")))
} else {
  if (NEEDS_POPULATION) {
    message("  No population data - population-based denominators will not be available")
  }
  population <- NULL
  pop_years <- integer(0)
}

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
  distinct(period_id, year, month) %>%
  arrange(period_id)

if (HAS_POPULATION) {
  min_pop_year <- min(pop_years)
  max_pop_year <- max(pop_years)
  available_periods <- available_periods %>%
    filter(year >= min_pop_year - 1, year <= max_pop_year + 1)
  if (nrow(available_periods) == 0) {
    stop("ERROR: No periods found within range of population data!")
  }
  message(sprintf("Found %d period(s) within population year range (%d-%d +/- 1 year extrapolation)",
                  nrow(available_periods), min_pop_year, max_pop_year))
} else {
  message(sprintf("Found %d period(s) with indicator data", nrow(available_periods)))
}

# Helper functions -----------------------------------------------------------------------------------------
all_possible_geo_cols <- c("admin_area_2", "admin_area_3", "admin_area_4")
geo_cols_in_adjusted <- intersect(all_possible_geo_cols, names(adjusted_data))

message(sprintf("Adjusted data columns: %s", paste(names(adjusted_data), collapse=", ")))
message(sprintf("Geo cols in adjusted data: %s", paste(geo_cols_in_adjusted, collapse=", ")))

if (HAS_POPULATION) {
  geo_cols_in_population <- intersect(all_possible_geo_cols, names(population))
  message(sprintf("Population columns: %s", paste(names(population), collapse=", ")))
  message(sprintf("Geo cols in population: %s", paste(geo_cols_in_population, collapse=", ")))
  geo_cols <- intersect(geo_cols_in_adjusted, geo_cols_in_population)
  if (length(geo_cols) == 0) {
    stop("ERROR: No common geographic columns between adjusted data and population data!")
  }
  message(sprintf("Using common geo_cols for aggregation: %s", paste(geo_cols, collapse=", ")))
} else {
  geo_cols <- geo_cols_in_adjusted
  message(sprintf("Using geo_cols from adjusted data: %s", paste(geo_cols, collapse=", ")))
}

aggregate_to_period <- function(data, target_period_id) {
  data %>%
    filter(period_id == target_period_id) %>%
    group_by(across(all_of(geo_cols)), indicator_common_id) %>%
    summarise(count = sum(.data[[SELECTED_COUNT_VARIABLE]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = indicator_common_id, values_from = count)
}

finest_geo_col <- geo_cols[length(geo_cols)]
message(sprintf("Finest common geo column: %s", finest_geo_col))

PERIOD_FRACTION <- 1/12

if (HAS_POPULATION) {
  pop_type_cols <- setdiff(names(population), c(geo_cols, "year"))

  safe_max <- function(x) if (length(x) == 0) NA_integer_ else max(x)
  safe_min <- function(x) if (length(x) == 0) NA_integer_ else min(x)

  get_interpolated_population <- function(pop_data, target_year, target_month, area_names) {
    month_fraction <- (target_month - 1) / 12

    years_available <- sort(unique(pop_data$year))
    if (length(years_available) == 0) {
      stop("ERROR: No population years available!")
    }

    year_before <- safe_max(years_available[years_available <= target_year])
    year_after <- safe_min(years_available[years_available > target_year])

    if (is.na(year_before) && is.na(year_after)) {
      stop(sprintf("ERROR: No population data available for interpolation (target year: %d)", target_year))
    }

    pop_for_areas <- pop_data %>%
      filter(.data[[finest_geo_col]] %in% area_names)

    if (is.na(year_before)) {
      y1 <- year_after
      y2 <- safe_min(years_available[years_available > y1])
      if (is.na(y2)) {
        return(pop_for_areas %>% filter(year == y1) %>% select(all_of(c(finest_geo_col, pop_type_cols))))
      }
      pop_y1 <- pop_for_areas %>% filter(year == y1)
      pop_y2 <- pop_for_areas %>% filter(year == y2)
      years_back <- y1 - target_year + month_fraction

      result <- pop_y1 %>%
        left_join(pop_y2, by = finest_geo_col, suffix = c("_y1", "_y2"))
      for (col in pop_type_cols) {
        col_y1 <- paste0(col, "_y1")
        col_y2 <- paste0(col, "_y2")
        if (col_y1 %in% names(result) && col_y2 %in% names(result)) {
          annual_growth <- (result[[col_y2]] / result[[col_y1]]) ^ (1 / (y2 - y1))
          result[[col]] <- result[[col_y1]] / (annual_growth ^ years_back)
        }
      }
      return(result %>% select(all_of(c(finest_geo_col, pop_type_cols))))
    }

    if (is.na(year_after)) {
      y1 <- safe_max(years_available[years_available < year_before])
      y2 <- year_before
      if (is.na(y1)) {
        return(pop_for_areas %>% filter(year == y2) %>% select(all_of(c(finest_geo_col, pop_type_cols))))
      }
      pop_y1 <- pop_for_areas %>% filter(year == y1)
      pop_y2 <- pop_for_areas %>% filter(year == y2)
      years_forward <- target_year - y2 + month_fraction

      result <- pop_y1 %>%
        left_join(pop_y2, by = finest_geo_col, suffix = c("_y1", "_y2"))
      for (col in pop_type_cols) {
        col_y1 <- paste0(col, "_y1")
        col_y2 <- paste0(col, "_y2")
        if (col_y1 %in% names(result) && col_y2 %in% names(result)) {
          annual_growth <- (result[[col_y2]] / result[[col_y1]]) ^ (1 / (y2 - y1))
          result[[col]] <- result[[col_y2]] * (annual_growth ^ years_forward)
        }
      }
      return(result %>% select(all_of(c(finest_geo_col, pop_type_cols))))
    }

    pop_before <- pop_for_areas %>% filter(year == year_before)
    pop_after <- pop_for_areas %>% filter(year == year_after)

    total_years <- year_after - year_before
    years_from_before <- target_year - year_before + month_fraction
    weight_after <- years_from_before / total_years
    weight_before <- 1 - weight_after

    result <- pop_before %>%
      left_join(pop_after, by = finest_geo_col, suffix = c("_before", "_after"))

    for (col in pop_type_cols) {
      col_before <- paste0(col, "_before")
      col_after <- paste0(col, "_after")
      if (col_before %in% names(result) && col_after %in% names(result)) {
        result[[col]] <- result[[col_before]] * weight_before + result[[col_after]] * weight_after
      }
    }

    result %>% select(all_of(c(finest_geo_col, pop_type_cols)))
  }

  get_population_for_period <- function(pop_data, target_period_id, area_data) {
    area_names <- unique(area_data[[finest_geo_col]])
    target_year <- as.integer(substr(as.character(target_period_id), 1, 4))
    target_month <- as.integer(substr(as.character(target_period_id), 5, 6))

    pop_areas <- unique(pop_data[[finest_geo_col]])
    matching <- intersect(area_names, pop_areas)
    if (length(matching) == 0) {
      stop(sprintf(
        "ERROR: No matching areas between indicator data and population data!\n  Indicator data %s values (sample): %s\n  Population data %s values (sample): %s\n  This usually means the datasets are from different countries or have mismatched admin area names.",
        finest_geo_col, paste(head(area_names, 3), collapse=", "),
        finest_geo_col, paste(head(pop_areas, 3), collapse=", ")
      ))
    }

    get_interpolated_population(pop_data, target_year, target_month, area_names)
  }
}

# Build scorecard rows for one period ----------------------------------------------------------------------
build_scorecard_for_period <- function(data, target_period_id) {
  agg <- aggregate_to_period(data, target_period_id)

  if (nrow(agg) == 0) {
    message(sprintf("  No data for period %s, skipping", target_period_id))
    return(tibble())
  }

  if (HAS_POPULATION) {
    pop <- get_population_for_period(population, target_period_id, agg)
    data <- agg %>%
      left_join(pop, by = finest_geo_col) %>%
      mutate(period_id = target_period_id)
  } else {
    data <- agg %>%
      mutate(period_id = target_period_id)
  }

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
