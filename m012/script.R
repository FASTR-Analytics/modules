SELECTED_COUNT_VARIABLE <- SELECTEDCOUNT
ADJUSTED_DATA_FILE <- "M2_adjusted_data.csv"

#-------------------------------------------------------------------------------------------------------------
# M12: Indicator values
#
# Materialises the additive INGREDIENTS of every common indicator at
# admin area x month grain. It does NOT compute any indicator's value: the
# formula is catalog data and is applied after aggregation, so that a chart at
# any grouping re-sums the ingredients and evaluates the formula once, exactly.
#
# The ingredient table says which base indicator fills which slot for which
# indicator. This script never parses a formula - it only sums the columns the
# table names. The table is DATA the app substitutes in below (the
# INDICATOR_INGREDIENTS token); the logic here is the same in every country.
#
# INPUTS:
#   M2_adjusted_data.csv        - facility x month x indicator, four count variants
#   INDICATOR_INGREDIENTS       - substituted tribble of
#                                 indicator_common_id, slot, ingredient_common_id
#
# OUTPUT:
#   M12_indicator_values.csv    - indicator x month x area, ing1..ing8
#-------------------------------------------------------------------------------------------------------------

message("Starting M12 indicator values module...")

library(dplyr)
library(readr)
library(tidyr)

SLOTS <- paste0("ing", 1:8)

message("Loading adjusted data from M2...")
adjusted_data <- read_csv(ADJUSTED_DATA_FILE, show_col_types = FALSE)

if (!SELECTED_COUNT_VARIABLE %in% names(adjusted_data)) {
  stop(sprintf(
    "ERROR: count variable '%s' is not a column of %s",
    SELECTED_COUNT_VARIABLE, ADJUSTED_DATA_FILE
  ))
}

# Substituted by the app from the run's resolved indicator catalog. An empty
# table (nothing mapped) is a valid 0-row tibble and yields a header-only
# output rather than an error.
ingredients <- INDICATOR_INGREDIENTS
message(sprintf("Ingredient table: %d row(s)", nrow(ingredients)))

bad_slots <- setdiff(unique(ingredients$slot), SLOTS)
if (length(bad_slots) > 0) {
  stop(sprintf(
    "ERROR: ingredient table names unknown slot(s): %s",
    paste(bad_slots, collapse = ", ")
  ))
}

# Admin columns present in the upstream output, finest last.
all_geo_cols <- c("admin_area_2", "admin_area_3", "admin_area_4")
geo_cols <- intersect(all_geo_cols, names(adjusted_data))
if (length(geo_cols) == 0) {
  stop("ERROR: no admin area columns in the adjusted data")
}
message(sprintf("Aggregating to: %s x period_id", paste(geo_cols, collapse = " x ")))

# Step 1: facilities summed away. This is the ONLY aggregation the module does;
# every later grouping re-sums these same additive numbers.
area_month <- adjusted_data %>%
  group_by(across(all_of(geo_cols)), period_id, indicator_common_id) %>%
  summarise(count = sum(.data[[SELECTED_COUNT_VARIABLE]], na.rm = TRUE), .groups = "drop")

# An ingredient with no rows in this dataset is NOT an error (PLAN_1a §1.5):
# the join below simply produces no row for it and the pivot leaves NA, which
# is the correct answer and what the app's evaluator expects. Failing here
# would abort generation on every instance that does not collect one of the
# seeded default indicators.
missing <- setdiff(unique(ingredients$ingredient_common_id),
                   unique(area_month$indicator_common_id))
if (length(missing) > 0) {
  message(sprintf(
    "Note: no data this run for ingredient indicator(s): %s",
    paste(missing, collapse = ", ")
  ))
}

# Step 2: one row per (indicator, area, month), ingredients in their slots.
# The join fans each base indicator's rows out to every indicator that uses it,
# relabelled by the slot it fills there.
message(sprintf("Building ingredient columns for %d indicator(s)...",
                length(unique(ingredients$indicator_common_id))))

# Renamed BEFORE the join so the two indicator columns never collide: the
# ingredient table's own key becomes target_indicator_id, and its
# ingredient_common_id becomes the join key.
ingredient_map <- ingredients %>%
  select(indicator_common_id, slot, ingredient_common_id) %>%
  rename(
    target_indicator_id = indicator_common_id,
    indicator_common_id = ingredient_common_id
  )

# Many-to-many by construction and by intent: one base indicator feeds several
# indicators, and each indicator has many area x month rows.
output <- area_month %>%
  inner_join(
    ingredient_map,
    by = "indicator_common_id",
    relationship = "many-to-many"
  ) %>%
  select(-indicator_common_id) %>%
  rename(indicator_common_id = target_indicator_id) %>%
  group_by(across(all_of(geo_cols)), period_id, indicator_common_id, slot) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = slot, values_from = count)

# Every declared slot column must exist, even when no indicator uses it, so the
# output schema is the same shape in every country.
for (slot in SLOTS) {
  if (!slot %in% names(output)) {
    output[[slot]] <- NA_real_
  }
}

output <- output %>%
  select(indicator_common_id, period_id, all_of(geo_cols), all_of(SLOTS))
# Base-R ordering: arrange(across(...)) is deprecated, and the sort columns are
# only known at runtime.
output <- output[
  do.call(order, output[c("indicator_common_id", "period_id", geo_cols)]),
]

message(sprintf("Writing %d row(s) to M12_indicator_values.csv", nrow(output)))
write_csv(output, "M12_indicator_values.csv", na = "NA")

message("M12 complete.")
