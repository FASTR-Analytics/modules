SELECTEDCOUNT <- "count_final_outliers"  # local development only: the app strips everything above the #--- marker and substitutes the tokens in the body
POPULATION_PERSON_YEARS <- "population.csv"
POPULATION_ACTIVE <- TRUE
#-------------------------------------------------------------------------------------------------------------
# M12: Indicator values
#
# Materialises the additive INGREDIENTS of every common indicator at
# admin area x month grain, and decides WHICH ROWS EXIST. When no indicator
# formula names a population the grain is the data's own admin level. When
# one does, the grain is the instance's POPULATION LEVEL, the admin columns
# of the person-years file: facilities and any finer admin level in the data
# are summed away.
#
# It does NOT write any indicator's value: the formula is applied after
# aggregation by the app, so that a chart at any grouping re-sums the
# ingredients and evaluates the formula once, exactly. It DOES evaluate the
# formula per row, for one purpose:
#
#   THE RULE: a row (indicator x area x month) is in the output only if the
#   indicator's formula over that row's ingredients produces a number.
#
# So an ingredient no facility in the area ever reports, a zero denominator,
# a nullif that fires, or a month or area the population store does not
# cover all mean the row is absent, and an indicator with no such row is
# absent from the package altogether. What the app cannot compute, nobody
# sees: not in a chart, not in a filter list. Keeping such a row instead
# would let a coarser grouping sum the numerator over cells the
# denominator never covers. A `coalesce` in the formula is honoured, because
# the formula decides, not the presence of a slot.
#
# Two tables are DATA the app substitutes in below; the logic here is the
# same in every country:
#   INDICATOR_INGREDIENTS  - which base indicator (or population) fills which
#                            slot column of which indicator
#   INDICATOR_EXPRESSIONS  - each indicator's formula over its slot columns
#                            (ing1..ing8), in the app's expression language
#
# The expression language's syntax is a subset of R's: decimals, + - * /,
# unary minus, parentheses, and calls to abs, coalesce, nullif. The text is
# evaluated as R inside an environment that binds `/`, coalesce and nullif to
# the app evaluator's semantics (lib/indicator_expression/evaluate.ts):
#   - NA in, NA out, through every operator and abs
#   - division by zero is NA, not Inf/NaN
#   - coalesce(a, b, ...) is the first non-NA argument
#   - nullif(a, b) is NA where a equals b, else a
# A row is kept when the result is finite. The app's parity test runs this
# script against the TypeScript evaluator over the same inputs.
#
# Population rates: the app expands the instance's annual population counts
# into monthly PERSON-YEARS (population / 12) per area, which sum like any
# count. They enter here as one more ingredient under the population type id
# (population_total, ...) - the same id the ingredient table names in
# whichever slot the term was assigned (slots follow order of appearance in
# the formula) - so the join below treats them exactly like a base indicator.
# This script never tells a population ingredient from a base one.
# The file holds rows only for the cells the stored population covers (its
# anchored years plus one year of extrapolation either side, per area), so
# coverage is partial by design and the app records it in the run manifest.
#
# INPUTS:
#   M2_adjusted_data.csv        - facility x month x indicator, four count variants
#   POPULATION_ACTIVE           - substituted TRUE/FALSE: whether any indicator
#                                 formula names a population
#   POPULATION_PERSON_YEARS     - area x month x population_type, person_years
#                                 at the population level, only the covered
#                                 cells (header-only when not active)
#   INDICATOR_INGREDIENTS       - substituted tribble of
#                                 indicator_common_id, slot, ingredient_common_id
#   INDICATOR_EXPRESSIONS       - substituted tribble of
#                                 indicator_common_id, expression
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
adjusted_data <- read_csv("M2_adjusted_data.csv", show_col_types = FALSE)

if (!SELECTEDCOUNT %in% names(adjusted_data)) {
  stop(sprintf(
    "ERROR: count variable '%s' is not a column of the adjusted data",
    SELECTEDCOUNT
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

expressions <- INDICATOR_EXPRESSIONS
message(sprintf("Expression table: %d row(s)", nrow(expressions)))

# The two tables describe the same indicators: one is the catalog's slot maps,
# the other its expressions, both written only for indicators with data.
unmatched <- c(
  setdiff(unique(ingredients$indicator_common_id), expressions$indicator_common_id),
  setdiff(expressions$indicator_common_id, unique(ingredients$indicator_common_id))
)
if (length(unmatched) > 0) {
  stop(sprintf(
    "ERROR: ingredient and expression tables disagree on indicator(s): %s",
    paste(unmatched, collapse = ", ")
  ))
}

# The evaluator's semantics, bound for eval() below. The helpers are closures
# over this script's environment, not formula_env, so arithmetic inside them
# is base R's; `/` says so explicitly because it is the one being rebound.
formula_env <- new.env(parent = baseenv())
formula_env[["/"]] <- function(x, y) {
  out <- base::`/`(x, y)
  out[!is.na(y) & y == 0] <- NA_real_
  out
}
formula_env[["coalesce"]] <- function(...) {
  args <- list(...)
  out <- args[[1]]
  for (a in args[-1]) {
    out <- ifelse(is.na(out), a, out)
  }
  out
}
formula_env[["nullif"]] <- function(x, y) {
  n <- max(length(x), length(y))
  x <- rep_len(x, n)
  y <- rep_len(y, n)
  x[!is.na(x) & !is.na(y) & x == y] <- NA_real_
  x
}
compiled <- setNames(
  lapply(expressions$expression, function(e) parse(text = e, keep.source = FALSE)[[1]]),
  expressions$indicator_common_id
)

all_geo_cols <- c("admin_area_2", "admin_area_3", "admin_area_4")
data_geo_cols <- intersect(all_geo_cols, names(adjusted_data))
if (length(data_geo_cols) == 0) {
  stop("ERROR: no admin area columns in the adjusted data")
}

population_active <- POPULATION_ACTIVE

# When a formula names a population, the person-years file's admin columns
# set this module's grain: the app writes its header at the instance's
# population level. The data may be finer (it is summed up) but never
# coarser: the app refuses such a run before this script runs, so the check
# below is defensive. Otherwise the data keeps its own admin level.
if (population_active) {
  message("Loading population person-years...")
  population <- read_csv(POPULATION_PERSON_YEARS, show_col_types = FALSE,
                         col_types = cols(.default = col_character(),
                                          period_id = col_integer(),
                                          person_years = col_double()))
  geo_cols <- intersect(all_geo_cols, names(population))
  if (length(geo_cols) == 0) {
    stop("ERROR: no admin area columns in the population file")
  }
  deeper_than_data <- setdiff(geo_cols, data_geo_cols)
  if (length(deeper_than_data) > 0) {
    stop(sprintf(
      "ERROR: the population level is deeper than the data: the population file has %s, which the adjusted data does not",
      paste(deeper_than_data, collapse = ", ")
    ))
  }
} else {
  message("No indicator formula names a population: keeping the data's own admin level")
  geo_cols <- data_geo_cols
}
message(sprintf("Aggregating to: %s x period_id", paste(geo_cols, collapse = " x ")))

# Step 1: facilities (and any admin level below the population level) summed
# away. This is the ONLY aggregation the module does; every later grouping
# re-sums these same additive numbers.
area_month <- adjusted_data %>%
  group_by(across(all_of(geo_cols)), period_id, indicator_common_id) %>%
  summarise(count = sum(.data[[SELECTEDCOUNT]], na.rm = TRUE), .groups = "drop")

# Step 1b: person-years join the area x month table as pseudo-indicator rows,
# already at geo_cols by construction. One log line per population type the
# file covers; a referenced type the file lacks shows up in the missing note
# below like any other ingredient with no rows.
if (population_active) {
  for (pop_type in unique(population$population_type)) {
    pop_rows <- population[population$population_type == pop_type, ]
    message(sprintf("  %s: %d person-year row(s), %d area(s), %d to %d",
                    pop_type, nrow(pop_rows), nrow(unique(pop_rows[geo_cols])),
                    min(pop_rows$period_id), max(pop_rows$period_id)))
  }
  area_month <- bind_rows(
    area_month,
    population %>%
      transmute(
        across(all_of(geo_cols)),
        period_id,
        indicator_common_id = population_type,
        count = person_years
      )
  )
}

# An ingredient with no rows in this dataset is NOT an error (PLAN_1a §1.5):
# the join below simply produces no row for it and the pivot leaves NA, which
# is the correct answer and what the app's evaluator expects. Failing here
# would abort generation on every instance that does not collect one of the
# seeded default indicators.
missing <- setdiff(
  unique(ingredients$ingredient_common_id),
  unique(area_month$indicator_common_id)
)
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

# Step 3: THE RULE (header). Evaluate each indicator's formula over its slot
# columns and keep the rows where the result is a number. The value itself is
# not written: the app re-evaluates over sums at whatever grouping a chart
# asks for. An indicator with no surviving row leaves the output entirely.
rows_before <- nrow(output)
if (rows_before > 0) {
  kept <- lapply(split(output, output$indicator_common_id), function(rows) {
    id <- rows$indicator_common_id[1]
    value <- eval(compiled[[id]], envir = rows[SLOTS], enclos = formula_env)
    value <- rep_len(value, nrow(rows))
    dropped <- sum(!is.finite(value))
    if (dropped > 0) {
      message(sprintf("  %s: %d of %d row(s) produce no value, dropped%s",
                      id, dropped, nrow(rows),
                      if (dropped == nrow(rows)) " (indicator absent from output)" else ""))
    }
    rows[is.finite(value), ]
  })
  output <- bind_rows(kept)
}
message(sprintf("Dropped %d row(s) whose formula produces no value", rows_before - nrow(output)))

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
