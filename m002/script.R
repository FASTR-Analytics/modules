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

# Load
raw_data         <- fread(PROJECT_DATA_HMIS)
outlier_data     <- fread("M1_output_outliers.csv")
completeness_data<- fread("M1_output_completeness.csv")

setDT(raw_data); setDT(outlier_data); setDT(completeness_data)

# Identify low-volume indicators (no observations with count >= 100) - excluded from all adjustments
low_volume_check <- raw_data[, .(has_volume = any(count >= 100, na.rm = TRUE)), by = indicator_common_id]
low_volume_check[, low_volume_exclude := !has_volume]
LOW_VOLUME_INDICATORS <- low_volume_check[has_volume == FALSE, indicator_common_id]

message("Indicators excluded from adjustment (volume < 100): ",
        if (length(LOW_VOLUME_INDICATORS) > 0) paste(LOW_VOLUME_INDICATORS, collapse = ", ") else "None")

# Geo columns
geo_cols <- grep("^admin_area_[0-9]+$", names(raw_data), value = TRUE)

# ----------------------------- Adjustment core --------------------------------------------------------------
# Memory notes: the facility-level table can exceed 80 million rows, so temporary columns are created one at a
# time and dropped as soon as they are used, and the same-month-last-year lookup is built only for the series
# that need it. Results are identical to computing everything at once.

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

apply_adjustments <- function(raw_data, completeness_data, outlier_data,
                              adjust_outliers = FALSE, adjust_completeness = FALSE) {
  message("Running adjustments...")
  
  # Merge inputs (period_id only)
  data_adj <- merge(
    completeness_data[, .(facility_id, indicator_common_id, period_id, completeness_flag)],
    outlier_data[, .(facility_id, indicator_common_id, period_id, outlier_flag)],
    by = c("facility_id", "indicator_common_id", "period_id"),
    all.x = TRUE
  )
  data_adj[, outlier_flag := fifelse(is.na(outlier_flag), 0L, outlier_flag)]
  
  data_adj <- merge(
    data_adj,
    raw_data[, .(facility_id, indicator_common_id, period_id, count)],
    by = c("facility_id", "indicator_common_id", "period_id"),
    all.x = TRUE
  )
  
  # period_id (YYYYMM) orders the same way as a date; month and year are read from it directly
  data_adj[, period_id := as.integer(period_id)]
  setorder(data_adj, facility_id, indicator_common_id, period_id)
  data_adj[, `:=`(mm = period_id %% 100L, yy = period_id %/% 100L)]
  
  data_adj[, `:=`(count_working = as.numeric(count),
                  adj_method = NA_character_, adjust_note = NA_character_)]
  
  # -------- Outlier adjustment --------
  if (adjust_outliers) {
    message(" -> Adjusting outliers...")
    data_adj[, valid_count := fifelse(outlier_flag == 0L & !is.na(count), count, NA_real_)]
    roll_fill(data_adj, "valid_count", "roll6", "center")
    roll_fill(data_adj, "valid_count", "fwd6",  "left")
    roll_fill(data_adj, "valid_count", "bwd6",  "right")

    data_adj[outlier_flag == 1L & !is.na(roll6),                        `:=`(count_working = roll6, adj_method = "roll6")]
    data_adj[outlier_flag == 1L & is.na(roll6) & !is.na(fwd6),          `:=`(count_working = fwd6, adj_method = "forward")]
    data_adj[outlier_flag == 1L & is.na(roll6) & is.na(fwd6) & !is.na(bwd6),
             `:=`(count_working = bwd6, adj_method = "backward")]
    data_adj[, c("roll6", "fwd6", "bwd6") := NULL]

    # same-month last year fallback
    hit <- smly_lookup(data_adj, data_adj[outlier_flag == 1L & is.na(adj_method),
                                          .(facility_id, indicator_common_id, mm, yy)])
    if (!is.null(hit) && nrow(hit) > 0L) {
      data_adj[hit, on = .(facility_id, indicator_common_id, mm, yy),
               `:=`(count_working = i.smly,
                    adj_method    = "same_month_last_year",
                    adjust_note   = sprintf("%04d-%02d", i.yy - 1L, i.mm))]
    }

    data_adj[, fallback := mean(valid_count, na.rm = TRUE), by = .(facility_id, indicator_common_id)]
    data_adj[outlier_flag == 1L & is.na(adj_method), `:=`(count_working = fallback, adj_method = "fallback")]
    data_adj[, c("fallback", "valid_count") := NULL]
    
    message("     Roll6 adjusted: ", sum(data_adj$adj_method == "roll6", na.rm = TRUE))
    message("     Forward-filled: ", sum(data_adj$adj_method == "forward", na.rm = TRUE))
    message("     Backward-filled:", sum(data_adj$adj_method == "backward", na.rm = TRUE))
    message("     Same-month LY:  ", sum(data_adj$adj_method == "same_month_last_year", na.rm = TRUE))
    message("     Fallback mean:  ", sum(data_adj$adj_method == "fallback", na.rm = TRUE))
  }
  
  # -------- Completeness adjustment --------
  if (adjust_completeness) {
    message(" -> Adjusting for completeness...")
    data_adj[, valid_count := fifelse(!is.na(count_working) & outlier_flag == 0L, count_working, NA_real_)]
    roll_fill(data_adj, "valid_count", "roll6", "center")
    roll_fill(data_adj, "valid_count", "fwd6",  "left")
    roll_fill(data_adj, "valid_count", "bwd6",  "right")

    data_adj[, adj_source := NA_character_]
    data_adj[is.na(count_working) & !is.na(roll6),                        `:=`(count_working = roll6, adj_source = "roll6")]
    data_adj[is.na(count_working) & is.na(roll6) & !is.na(fwd6),          `:=`(count_working = fwd6, adj_source = "forward")]
    data_adj[is.na(count_working) & is.na(roll6) & is.na(fwd6) & !is.na(bwd6),
             `:=`(count_working = bwd6, adj_source = "backward")]
    data_adj[, c("roll6", "fwd6", "bwd6") := NULL]

    hit <- smly_lookup(data_adj, data_adj[is.na(count_working), .(facility_id, indicator_common_id, mm, yy)])
    if (!is.null(hit) && nrow(hit) > 0L) {
      data_adj[hit, on = .(facility_id, indicator_common_id, mm, yy),
               `:=`(count_working = i.smly, adj_source = "same_month_last_year")]
    }

    data_adj[, fallback := mean(valid_count, na.rm = TRUE), by = .(facility_id, indicator_common_id)]
    data_adj[is.na(count_working),                                       `:=`(count_working = fallback, adj_source = "fallback")]

    message("     Roll6 filled:    ", sum(data_adj$adj_source == "roll6",   na.rm = TRUE))
    message("     Forward-filled:  ", sum(data_adj$adj_source == "forward", na.rm = TRUE))
    message("     Backward-filled: ", sum(data_adj$adj_source == "backward",na.rm = TRUE))
    message("     Same-month LY:   ", sum(data_adj$adj_source == "same_month_last_year", na.rm = TRUE))
    message("     Fallback mean:   ", sum(data_adj$adj_source == "fallback",na.rm = TRUE))

    data_adj[, c("valid_count", "fallback", "adj_source") := NULL]
  }
  
  data_adj[, c("mm", "yy") := NULL]
  return(data_adj[])
}

# ----------------------------- Scenarios wrapper ------------------------------------------------------------
apply_adjustments_scenarios <- function(raw_data, completeness_data, outlier_data) {
  message("Applying adjustments across scenarios...")
  join_cols <- c("facility_id","indicator_common_id","period_id")
  
  scenarios <- list(
    none          = list(adjust_outliers = FALSE, adjust_completeness = FALSE),
    outliers      = list(adjust_outliers = TRUE,  adjust_completeness = FALSE),
    completeness  = list(adjust_outliers = FALSE, adjust_completeness = TRUE),
    both          = list(adjust_outliers = TRUE,  adjust_completeness = TRUE)
  )
  
  results <- lapply(names(scenarios), function(scn) {
    message(" -> Scenario: ", scn)
    opts <- scenarios[[scn]]
    dat  <- apply_adjustments(raw_data, completeness_data, outlier_data,
                              adjust_outliers = opts$adjust_outliers,
                              adjust_completeness = opts$adjust_completeness)
    dat[grepl(EXCLUDED_PATTERN, indicator_common_id, ignore.case = TRUE) |
        indicator_common_id %in% LOW_VOLUME_INDICATORS, count_working := count]
    dat <- dat[, .(facility_id, indicator_common_id, period_id,
                   count_final = count_working)]
    setnames(dat, "count_final", paste0("count_final_", scn))
    dat
  })
  names(results) <- names(scenarios)
  
  Reduce(function(x, y) merge(x, y, by = join_cols, all = TRUE), results)
}

# ----------------------------- Main -------------------------------------------------------------------------
message("Running adjustments analysis...")

adjusted_data_final <- apply_adjustments_scenarios(
  raw_data = raw_data,
  completeness_data = completeness_data,
  outlier_data = outlier_data
)

# Metadata lookups
geo_lookup <- unique(raw_data[, .SD, .SDcols = c("facility_id", geo_cols)])

# Merge metadata into facility-level adjusted data
setDT(adjusted_data_final)
setkey(adjusted_data_final, facility_id)
setkey(geo_lookup, facility_id)

adjusted_data_export <- merge(adjusted_data_final, geo_lookup, by = "facility_id", all.x = TRUE)

# Geo sets
geo_cols <- grep("^admin_area_[0-9]+$", names(adjusted_data_export), value = TRUE)
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
# Drop admin_area_1 from facility-level file for cleanliness (unchanged)
adjusted_data_export_clean <- adjusted_data_export[, !"admin_area_1"]

fwrite(adjusted_data_export_clean,     "M2_adjusted_data.csv",            na = "NA")
fwrite(adjusted_data_admin_area_final, "M2_adjusted_data_admin_area.csv", na = "NA")
fwrite(adjusted_data_national_final,   "M2_adjusted_data_national.csv",   na = "NA")
fwrite(low_volume_check[, .(indicator_common_id, low_volume_exclude)], "M2_low_volume_exclusions.csv", na = "NA")

message("Adjustments completed and all outputs saved.")
