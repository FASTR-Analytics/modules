#-------------------------------------------------------------------------------------------------------------
# M10: Health Facility Assessment Module
#
# Computes HFA indicators from facility-level survey data. Indicator R code,
# warnings, and metadata are dynamically generated based on the configured
# HFA indicator catalog.
#
# OUTPUT:
#   M10_hfa_results.csv — One row per facility × indicator × time_point
#   M10_hfa_response_status.csv — Response-status companion (don't-know /
#     missing / not-applicable / answered), one row per facility × indicator
#     × time_point. Empty (header only) when the server does not emit
#     __status columns.
#
#-------------------------------------------------------------------------------------------------------------

library(dplyr)
library(tidyr)

print("Starting HFA script...")

__WARNING_PRINTS__

# Read and pivot data to wide format
data <- read.csv(PROJECT_DATA_HFA)

# Sampling weight column (added by newer exports; tolerate its absence)
if (!"weight" %in% names(data)) data$weight <- NA_real_

data_wide <- data %>%
  pivot_wider(names_from = var_name, values_from = value)

# Detect facility columns dynamically
facility_cols <- names(data_wide)[grepl("^(facility_|admin_area_|time_point)", names(data_wide))]

# Convert pivoted variable columns to numeric (they may be character after pivot)
data_wide <- data_wide %>%
  mutate(across(-all_of(facility_cols), as.numeric))

# Resolve sampling weights: weights off, or a facility with no weight row,
# means weight 1 — the weighted formulas then reduce exactly to unweighted ones
if (USE_SAMPLE_WEIGHTS) {
  n_no_weight <- sum(is.na(data_wide$weight))
  print(paste0("Facility/time-point rows without a sampling weight (fallback to weight 1): ", n_no_weight))
}
data_wide <- data_wide %>%
  mutate(weight_final = if (USE_SAMPLE_WEIGHTS) coalesce(weight, 1) else 1)

# Calculate indicators
results <- data_wide %>%
__INDICATOR_MUTATES__

# Select only indicator columns
indicator_cols <- c(__INDICATOR_COLS__)
results_final <- results %>%
  select(all_of(indicator_cols))

# Create indicator metadata mapping (category, type, aggregation)
indicator_metadata <- data.frame(
__INDICATOR_METADATA__
)

# Pivot back to long format and add metadata
facility_info <- data_wide %>%
  select(all_of(facility_cols))

results_long <- facility_info %>%
  bind_cols(data_wide %>% select(weight_final)) %>%
  bind_cols(results_final) %>%
  pivot_longer(
    cols = all_of(indicator_cols),
    names_to = "hfa_indicator",
    values_to = "raw_value"
  ) %>%
  left_join(indicator_metadata, by = "hfa_indicator") %>%
  filter(!is.na(raw_value)) %>%
  mutate(
    sum_val    = ifelse(ind_aggregation == "sum", raw_value * weight_final, NA_real_),
    avg_num    = ifelse(ind_aggregation == "avg", raw_value * weight_final, NA_real_),
    avg_weight = ifelse(ind_aggregation == "avg", weight_final, NA_real_)
  ) %>%
  select(-raw_value, -weight_final, -hfa_short_label, -ind_type, -ind_aggregation)

if (nrow(results_long) == 0) {
  stop("No results generated - all indicator values are NA. Check that HFA indicators have been configured with R code.")
}

# Write output
write.csv(results_long, "M10_hfa_results.csv", row.names = FALSE)

# Response-status companion output: share of facilities answering don't-know /
# missing / not-applicable per indicator. The __status columns are emitted by
# the server-side script generator; without them (older server) write a
# header-only file so ingestion still succeeds.
status_cols <- names(results)[grepl("__status$", names(results))]
if (length(status_cols) > 0) {
  status_long <- facility_info %>%
    bind_cols(data_wide %>% select(weight_final)) %>%
    bind_cols(results %>% select(all_of(status_cols))) %>%
    pivot_longer(
      cols = all_of(status_cols),
      names_to = "hfa_indicator",
      values_to = "response_status"
    ) %>%
    mutate(hfa_indicator = sub("__status$", "", hfa_indicator)) %>%
    left_join(indicator_metadata, by = "hfa_indicator") %>%
    filter(!is.na(response_status)) %>%
    mutate(
      dk_num       = (response_status == "dont_know") * weight_final,
      missing_num  = (response_status == "missing") * weight_final,
      answered_num = (response_status == "answered") * weight_final,
      na_num       = (response_status == "not_applicable") * weight_final,
      resp_weight  = (response_status != "not_applicable") * weight_final,
      total_weight = weight_final
    ) %>%
    select(-response_status, -weight_final, -hfa_short_label, -ind_type, -ind_aggregation)
  write.csv(status_long, "M10_hfa_response_status.csv", row.names = FALSE)
} else {
  writeLines(
    "facility_id,admin_area_4,admin_area_3,admin_area_2,admin_area_1,hfa_indicator,hfa_category,hfa_sub_category,hfa_service_category,time_point,facility_ownership,facility_type,facility_custom_1,facility_custom_2,facility_custom_3,facility_custom_4,facility_custom_5,dk_num,missing_num,answered_num,na_num,resp_weight,total_weight",
    "M10_hfa_response_status.csv"
  )
}

print("HFA script completed successfully!")
