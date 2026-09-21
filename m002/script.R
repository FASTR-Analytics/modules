COUNTRY_ISO3 <- "ZMB"
PROJECT_DATA_HMIS <- "hmis_ZMB.csv"

#-------------------------------------------------------------------------------------------------------------
# CB - R code FASTR PROJECT
# Module: DATA QUALITY ADJUSTMENT
# Last edit: 2026 Sep 21
#-------------------------------------------------------------------------------------------------------------

# -------------------------- KEY OUTPUT ----------------------------------------------------------------------
# FILE: M2_adjusted_data.csv              # Facility-level adjusted volumes (all scenarios), period_id only
# FILE: M2_adjusted_data_admin_area.csv   # Admin-level adjusted volumes (all scenarios), period_id only
# FILE: M2_adjusted_data_national.csv     # National-level adjusted volumes (all scenarios), period_id only

# Libraries --------------------------------------------------------------------------------------------------
library(data.table)
library(zoo)
library(lubridate)

EXCLUDED_PATTERN <- "death|still_birth"
WINDOW_MONTHS    <- 7   # rolling window width; odd so the centred window is symmetric
MIN_VALID_MONTHS <- 3   # minimum valid months required inside a window for it to be used

# Memory notes: the facility-level table can exceed 80 million rows. Only the columns that are used are read,
# the four scenarios are computed on one table (the outlier step is shared by the "outliers" and "both"
# scenarios), temporary columns are created one at a time and dropped as soon as they are used, and the
# same-month-last-year lookup is built only for the series that need it.

# Load (only the columns used)
hmis_cols <- names(fread(PROJECT_DATA_HMIS, nrows = 0))
geo_cols  <- grep("^admin_area_[0-9]+$", hmis_cols, value = TRUE)
raw_data         <- fread(PROJECT_DATA_HMIS,
                          select = c("facility_id", geo_cols, "period_id", "indicator_common_id", "count"))
outlier_data     <- fread("M1_output_outliers.csv",
                          select = c("facility_id", "indicator_common_id", "period_id", "outlier_flag"))
completeness_data<- fread("M1_output_completeness.csv",
                          select = c("facility_id", "indicator_common_id", "period_id", "completeness_flag"))

setDT(raw_data); setDT(outlier_data); setDT(completeness_data)

# Identify low-volume indicators (no observations with count >= 100) - excluded from all adjustments
low_volume_check <- raw_data[, .(has_volume = any(count >= 100, na.rm = TRUE)), by = indicator_common_id]
low_volume_check[, low_volume_exclude := !has_volume]
LOW_VOLUME_INDICATORS <- low_volume_check[has_volume == FALSE, indicator_common_id]

message("Indicators excluded from adjustment (volume < 100): ",
        if (length(LOW_VOLUME_INDICATORS) > 0) paste(LOW_VOLUME_INDICATORS, collapse = ", ") else "None")

# Facility -> admin area lookup, then the geo columns are dropped from the raw table
geo_lookup <- unique(raw_data[, .SD, .SDcols = c("facility_id", geo_cols)])
raw_data[, (geo_cols) := NULL]
if (anyDuplicated(geo_lookup, by = "facility_id")) {
  message("NOTE: some facilities carry more than one admin area combination; the first is used")
  geo_lookup <- unique(geo_lookup, by = "facility_id")
}

# ----------------------------- Adjustment core --------------------------------------------------------------

# Rolling mean over WINDOW_MONTHS valid months for one alignment, kept only where at least MIN_VALID_MONTHS
# of the window are valid.
roll_fill <- function(dt, valid_col, out_col, align) {
  dt[, (out_col) := frollmean(get(valid_col), WINDOW_MONTHS, na.rm = TRUE, align = align),
     by = .(facility_id, indicator_common_id)]
  dt[, n_valid := frollsum(as.integer(!is.na(get(valid_col))), WINDOW_MONTHS, align = align),
     by = .(facility_id, indicator_common_id)]
  dt[is.na(n_valid) | n_valid < MIN_VALID_MONTHS, (out_col) := NA_real_]
  dt[, n_valid := NULL]
  invisible(dt)
}

# Same-month-last-year value for the rows in `target` (facility_id, indicator_common_id, mm, yy): the raw,
# non-outlier count of the same month one year earlier, when exactly one such record exists.
smly_lookup <- function(dt, target) {
  if (nrow(target) == 0L) return(NULL)
  series <- unique(target[, .(facility_id, indicator_common_id)])
  src <- dt[series, on = .(facility_id, indicator_common_id), nomatch = NULL,
            .(facility_id, indicator_common_id, mm, yy = yy + 1L, count, outlier_flag)]
  src <- src[outlier_flag == 0L & !is.na(count)]
  src <- src[, .(smly = if (.N == 1L) count else NA_real_), by = .(facility_id, indicator_common_id, mm, yy)]
  hit <- src[target, on = .(facility_id, indicator_common_id, mm, yy), nomatch = NULL]
  hit[!is.na(smly)]
}

# Outlier adjustment: replaces flagged outliers in count_working (in place)
adjust_outliers_step <- function(dt) {
  message(" -> Adjusting outliers...")
  dt[, adj_method := NA_character_]
  dt[, valid_count := fifelse(outlier_flag == 0L & !is.na(count), count, NA_real_)]
  roll_fill(dt, "valid_count", "roll6", "center")
  roll_fill(dt, "valid_count", "fwd6",  "left")
  roll_fill(dt, "valid_count", "bwd6",  "right")

  dt[outlier_flag == 1L & !is.na(roll6),                        `:=`(count_working = roll6, adj_method = "roll6")]
  dt[outlier_flag == 1L & is.na(roll6) & !is.na(fwd6),          `:=`(count_working = fwd6, adj_method = "forward")]
  dt[outlier_flag == 1L & is.na(roll6) & is.na(fwd6) & !is.na(bwd6),
     `:=`(count_working = bwd6, adj_method = "backward")]
  dt[, c("roll6", "fwd6", "bwd6") := NULL]

  # same-month last year fallback
  hit <- smly_lookup(dt, dt[outlier_flag == 1L & is.na(adj_method), .(facility_id, indicator_common_id, mm, yy)])
  if (!is.null(hit) && nrow(hit) > 0L) {
    dt[hit, on = .(facility_id, indicator_common_id, mm, yy),
       `:=`(count_working = i.smly, adj_method = "same_month_last_year")]
  }

  dt[, fallback := mean(valid_count, na.rm = TRUE), by = .(facility_id, indicator_common_id)]
  dt[outlier_flag == 1L & is.na(adj_method), `:=`(count_working = fallback, adj_method = "fallback")]
  dt[, c("fallback", "valid_count") := NULL]

  message("     Roll6 adjusted: ", sum(dt$adj_method == "roll6", na.rm = TRUE))
  message("     Forward-filled: ", sum(dt$adj_method == "forward", na.rm = TRUE))
  message("     Backward-filled:", sum(dt$adj_method == "backward", na.rm = TRUE))
  message("     Same-month LY:  ", sum(dt$adj_method == "same_month_last_year", na.rm = TRUE))
  message("     Fallback mean:  ", sum(dt$adj_method == "fallback", na.rm = TRUE))
  dt[, adj_method := NULL]
  invisible(dt)
}

# Completeness adjustment: fills missing count_working (in place)
adjust_completeness_step <- function(dt) {
  message(" -> Adjusting for completeness...")
  dt[, valid_count := fifelse(!is.na(count_working) & outlier_flag == 0L, count_working, NA_real_)]
  roll_fill(dt, "valid_count", "roll6", "center")
  roll_fill(dt, "valid_count", "fwd6",  "left")
  roll_fill(dt, "valid_count", "bwd6",  "right")

  dt[, adj_source := NA_character_]
  dt[is.na(count_working) & !is.na(roll6),                        `:=`(count_working = roll6, adj_source = "roll6")]
  dt[is.na(count_working) & is.na(roll6) & !is.na(fwd6),          `:=`(count_working = fwd6, adj_source = "forward")]
  dt[is.na(count_working) & is.na(roll6) & is.na(fwd6) & !is.na(bwd6),
     `:=`(count_working = bwd6, adj_source = "backward")]
  dt[, c("roll6", "fwd6", "bwd6") := NULL]

  hit <- smly_lookup(dt, dt[is.na(count_working), .(facility_id, indicator_common_id, mm, yy)])
  if (!is.null(hit) && nrow(hit) > 0L) {
    dt[hit, on = .(facility_id, indicator_common_id, mm, yy),
       `:=`(count_working = i.smly, adj_source = "same_month_last_year")]
  }

  dt[, fallback := mean(valid_count, na.rm = TRUE), by = .(facility_id, indicator_common_id)]
  dt[is.na(count_working), `:=`(count_working = fallback, adj_source = "fallback")]

  message("     Roll6 filled:    ", sum(dt$adj_source == "roll6",   na.rm = TRUE))
  message("     Forward-filled:  ", sum(dt$adj_source == "forward", na.rm = TRUE))
  message("     Backward-filled: ", sum(dt$adj_source == "backward",na.rm = TRUE))
  message("     Same-month LY:   ", sum(dt$adj_source == "same_month_last_year", na.rm = TRUE))
  message("     Fallback mean:   ", sum(dt$adj_source == "fallback",na.rm = TRUE))
  dt[, c("valid_count", "fallback", "adj_source") := NULL]
  invisible(dt)
}

# ----------------------------- Scenarios --------------------------------------------------------------------
# All four scenarios on one table:
#   none         = raw count
#   outliers     = outlier step
#   both         = outlier step, then completeness step (continues from the outlier-adjusted values)
#   completeness = completeness step on the raw counts
# Excluded indicators (EXCLUDED_PATTERN, low volume) keep the raw count in every scenario.
apply_adjustments_scenarios <- function(raw_data, completeness_data, outlier_data) {
  message("Applying adjustments across scenarios...")

  # Row universe = the completeness table (one row per facility, indicator, period); outlier flags and raw
  # counts are attached in place, so no copy of the full table is made.
  keys <- c("facility_id", "indicator_common_id", "period_id")
  data_adj <- completeness_data[, .(facility_id, indicator_common_id, period_id)]
  data_adj[outlier_data, outlier_flag := i.outlier_flag, on = keys]
  data_adj[, outlier_flag := fifelse(is.na(outlier_flag), 0L, outlier_flag)]
  data_adj[raw_data, count := i.count, on = keys]

  # period_id (YYYYMM) orders the same way as a date; month and year are read from it directly
  data_adj[, period_id := as.integer(period_id)]
  setorder(data_adj, facility_id, indicator_common_id, period_id)
  data_adj[, `:=`(mm = period_id %% 100L, yy = period_id %/% 100L)]
  data_adj[, excluded := grepl(EXCLUDED_PATTERN, indicator_common_id, ignore.case = TRUE) |
                         indicator_common_id %in% LOW_VOLUME_INDICATORS]

  message(" -> Scenario: none")
  data_adj[, count_final_none := as.numeric(count)]

  message(" -> Scenario: outliers")
  data_adj[, count_working := as.numeric(count)]
  adjust_outliers_step(data_adj)
  data_adj[, count_final_outliers := fifelse(excluded, as.numeric(count), count_working)]

  message(" -> Scenario: both")
  adjust_completeness_step(data_adj)
  data_adj[, count_final_both := fifelse(excluded, as.numeric(count), count_working)]

  message(" -> Scenario: completeness")
  data_adj[, count_working := as.numeric(count)]
  adjust_completeness_step(data_adj)
  data_adj[, count_final_completeness := fifelse(excluded, as.numeric(count), count_working)]

  data_adj[, c("count_working", "count", "outlier_flag", "mm", "yy", "excluded") := NULL]
  setcolorder(data_adj, c("facility_id", "indicator_common_id", "period_id",
                          "count_final_none", "count_final_outliers", "count_final_completeness", "count_final_both"))
  data_adj[]
}

# ----------------------------- Main -------------------------------------------------------------------------
message("Running adjustments analysis...")

adjusted_data_export <- apply_adjustments_scenarios(
  raw_data = raw_data,
  completeness_data = completeness_data,
  outlier_data = outlier_data
)
rm(completeness_data, outlier_data, raw_data); invisible(gc())

# Attach admin areas to the facility-level adjusted data (in place)
adjusted_data_export[geo_lookup, on = "facility_id", (geo_cols) := mget(paste0("i.", geo_cols))]

# Geo sets
geo_admin_area_sub <- setdiff(geo_cols, "admin_area_1")

message("Detected admin area columns: ", paste(geo_cols, collapse = ", "))
message("Using for subnational aggregation: ", paste(geo_admin_area_sub, collapse = ", "))

# Order columns for export (no year/quarter)
setcolorder(adjusted_data_export, c(
  "facility_id",
  geo_admin_area_sub,
  "period_id",
  "indicator_common_id"
))

# --------------------------- Subnational Output (period_id only) --------------------------------------------
adjusted_data_admin_area_final <- adjusted_data_export[
  ,
  .(
    count_final_none         = sum(count_final_none,         na.rm = TRUE),
    count_final_outliers     = sum(count_final_outliers,     na.rm = TRUE),
    count_final_completeness = sum(count_final_completeness, na.rm = TRUE),
    count_final_both         = sum(count_final_both,         na.rm = TRUE)
  ),
  by = c(geo_admin_area_sub, "indicator_common_id", "period_id")
]

# --------------------------- National Output (period_id only) -----------------------------------------------
adjusted_data_national_final <- adjusted_data_export[
  ,
  .(
    count_final_none         = sum(count_final_none,         na.rm = TRUE),
    count_final_outliers     = sum(count_final_outliers,     na.rm = TRUE),
    count_final_completeness = sum(count_final_completeness, na.rm = TRUE),
    count_final_both         = sum(count_final_both,         na.rm = TRUE)
  ),
  by = .(admin_area_1, indicator_common_id, period_id)
]

# --------------------------- Save Outputs -------------------------------------------------------------------
fwrite(adjusted_data_admin_area_final, "M2_adjusted_data_admin_area.csv", na = "NA")
fwrite(adjusted_data_national_final,   "M2_adjusted_data_national.csv",   na = "NA")
fwrite(low_volume_check[, .(indicator_common_id, low_volume_exclude)], "M2_low_volume_exclusions.csv", na = "NA")

# Facility-level file without admin_area_1 (dropped in place, no copy of the table)
adjusted_data_export[, admin_area_1 := NULL]
fwrite(adjusted_data_export, "M2_adjusted_data.csv", na = "NA")

message("Adjustments completed and all outputs saved.")
