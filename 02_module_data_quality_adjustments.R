
PROJECT_DATA_HMIS <- "hmis_nigeria_q2.csv"

#-------------------------------------------------------------------------------------------------------------
# CB - R code FASTR PROJECT
# Module: DATA QUALITY ADJUSTMENT
# Last edit: 2025 Aug 5

# This script dynamically adjusts raw data for:
#   1. Outliers: Replaces flagged outliers with 12-month rolling averages (excluding outliers).
#   2. Completeness: Replaces missing count with 12-month rolling averages (excluding outliers).

# Ce script ajuste dynamiquement les données brutes pour :
#   1. Les valeurs aberrantes : Remplace les valeurs identifiées comme aberrantes par une moyenne mobile sur 12 mois (hors valeurs aberrantes).
#   2. L'exhaustivité : Remplace les valeurs manquantes par une moyenne mobile sur 12 mois (hors valeurs aberrantes).

# -------------------------- KEY OUTPUT ----------------------------------------------------------------------
# FILE: M2_adjusted_data.csv              # Dataset including facility-level adjusted volumes for all adjustment scenarios.
# FILE: M2_adjusted_data_admin_area.csv   # Dataset including admin-level adjusted volumes for all adjustment scenarios.
# FILE: M2_adjusted_data_national.csv     # Dataset including national-level adjusted volumes for all adjustment scenarios.


# Load Required Libraries -----------------------------------------------------------------------------------
library(data.table)
library(zoo)
library(lubridate)

EXCLUDED_FROM_ADJUSTMENT <- c("u5_deaths", "maternal_deaths")

# Load CSVs and convert to data.table
raw_data <- fread(PROJECT_DATA_HMIS)
outlier_data <- fread("M1_output_outliers.csv")
completeness_data <- fread("M1_output_completeness.csv")

setDT(raw_data)
setDT(outlier_data)
setDT(completeness_data)


# Define Functions ------------------------------------------------------------------------------------------
geo_cols <- colnames(raw_data)[grepl("^admin_area_[0-9]+$", colnames(raw_data))]

# Function to Apply Adjustments -----------------------------------------------------------------------------
apply_adjustments <- function(raw_data, completeness_data, outlier_data,
                              adjust_outliers = FALSE, adjust_completeness = FALSE) {
  message("Running adjustments...")
  
  # Merge inputs
  data_adj <- merge(completeness_data,
                    outlier_data[, .(facility_id, indicator_common_id, period_id, outlier_flag)],
                    by = c("facility_id", "indicator_common_id", "period_id"),
                    all.x = TRUE)
  data_adj[, outlier_flag := fifelse(is.na(outlier_flag), 0L, outlier_flag)]
  
  data_adj <- merge(data_adj,
                    raw_data[, .(facility_id, indicator_common_id, period_id, count)],
                    by = c("facility_id", "indicator_common_id", "period_id"),
                    all.x = TRUE)
  
  data_adj[, count_working := as.numeric(count)]
  
  # Create date column once at the beginning
  data_adj[, date := as.Date(paste0(substr(period_id, 1, 4), "-", substr(period_id, 6, 7), "-01"))]
  setorder(data_adj, facility_id, indicator_common_id, date)
  data_adj[, `:=`(adj_method = NA_character_, adjust_note = NA_character_)]
  
  # ---------------- Outlier Adjustment ----------------
  if (adjust_outliers) {
    message(" -> Adjusting outliers...")
    
    data_adj[, valid_count := fifelse(outlier_flag == 0L & !is.na(count), count, NA_real_)]
    data_adj[, `:=`(
      roll6 = frollmean(valid_count, 6, na.rm = TRUE, align = "center"),
      fwd6 = frollmean(valid_count, 6, na.rm = TRUE, align = "left"),
      bwd6 = frollmean(valid_count, 6, na.rm = TRUE, align = "right"),
      fallback = median(valid_count, na.rm = TRUE)
    ), by = .(facility_id, indicator_common_id)]
    
    # Fixed logic: Check each method sequentially for outliers
    data_adj[outlier_flag == 1L & !is.na(roll6), `:=`(count_working = roll6, adj_method = "roll6")]
    data_adj[outlier_flag == 1L & is.na(roll6) & !is.na(fwd6), `:=`(count_working = fwd6, adj_method = "forward")]
    data_adj[outlier_flag == 1L & is.na(roll6) & is.na(fwd6) & !is.na(bwd6), `:=`(count_working = bwd6, adj_method = "backward")]
    
    data_adj[, `:=`(month = month(date), year = year(date))]
    data_adj <- data_adj[, {
      for (i in which(outlier_flag == 1L & is.na(roll6) & is.na(fwd6) & is.na(bwd6))) {
        current_month <- month[i]
        current_year <- year[i]
        j <- which(month == current_month & year == (current_year - 1) & outlier_flag == 0L & !is.na(count))
        if (length(j) == 1L) {
          count_working[i] <- count[j]
          adj_method[i] <- "same_month_last_year"
          adjust_note[i] <- format(date[j], "%b-%Y")
        }
      }
      .SD
    }, by = .(facility_id, indicator_common_id)]
    
    data_adj[outlier_flag == 1L & is.na(adj_method), `:=`(count_working = fallback, adj_method = "fallback")]
    
    message("     Roll6 adjusted: ", sum(data_adj$adj_method == "roll6", na.rm = TRUE))
    message("     Forward-filled: ", sum(data_adj$adj_method == "forward", na.rm = TRUE))
    message("     Backward-filled: ", sum(data_adj$adj_method == "backward", na.rm = TRUE))
    message("     Same-month fallback: ", sum(data_adj$adj_method == "same_month_last_year", na.rm = TRUE))
    message("     Fallback (median): ", sum(data_adj$adj_method == "fallback", na.rm = TRUE))
    message("   -> Outlier adjustment complete")
    
    data_adj[, c("roll6", "fwd6", "bwd6", "fallback", "valid_count", "month", "year") := NULL]
  }
  
  # ---------------- Completeness Adjustment ----------------
  if (adjust_completeness) {
    message(" -> Adjusting for completeness...")
    
    # Date column already exists, no need to recreate it
    # setorder already done above
    
    data_adj[, valid_count := fifelse(!is.na(count_working) & count_working > 0, count_working, NA_real_)]
    data_adj[, `:=`(
      roll6 = frollmean(valid_count, 6, na.rm = TRUE, align = "center"),
      fwd6 = frollmean(valid_count, 6, na.rm = TRUE, align = "left"),
      bwd6 = frollmean(valid_count, 6, na.rm = TRUE, align = "right"),
      fallback = mean(valid_count, na.rm = TRUE)
    ), by = .(facility_id, indicator_common_id)]
    
    data_adj[, adj_source := NA_character_]
    # Fixed logic: Check each method sequentially for missing values
    data_adj[is.na(count_working) & !is.na(roll6), `:=`(count_working = roll6, adj_source = "roll6")]
    data_adj[is.na(count_working) & is.na(roll6) & !is.na(fwd6), `:=`(count_working = fwd6, adj_source = "forward")]
    data_adj[is.na(count_working) & is.na(roll6) & is.na(fwd6) & !is.na(bwd6), `:=`(count_working = bwd6, adj_source = "backward")]
    data_adj[is.na(count_working), `:=`(count_working = fallback, adj_source = "fallback")]
    
    message("     Roll6 adjusted: ", sum(data_adj$adj_source == "roll6", na.rm = TRUE))
    message("     Forward-filled: ", sum(data_adj$adj_source == "forward", na.rm = TRUE))
    message("     Backward-filled: ", sum(data_adj$adj_source == "backward", na.rm = TRUE))
    message("     Fallback (mean): ", sum(data_adj$adj_source == "fallback", na.rm = TRUE))
    
    data_adj[, c("valid_count", "roll6", "fwd6", "bwd6", "fallback", "adj_source") := NULL]
  }
  
  return(data_adj)
}

# Function to Apply Adjustments Across Scenarios ------------------------------------------------------------
apply_adjustments_scenarios <- function(raw_data, completeness_data, outlier_data) {
  message("Applying adjustments across scenarios...")
  
  join_cols <- c("facility_id", "indicator_common_id", "period_id")
  
  scenarios <- list(
    none = list(adjust_outliers = FALSE, adjust_completeness = FALSE),
    outliers = list(adjust_outliers = TRUE, adjust_completeness = FALSE),
    completeness = list(adjust_outliers = FALSE, adjust_completeness = TRUE),
    both = list(adjust_outliers = TRUE, adjust_completeness = TRUE)
  )
  
  results <- vector("list", length(scenarios))
  names(results) <- names(scenarios)
  
  for (scenario in names(scenarios)) {
    message(" -> Scenario: ", scenario)
    opts <- scenarios[[scenario]]
    
    dat <- apply_adjustments(raw_data, completeness_data, outlier_data,
                             adjust_outliers = opts$adjust_outliers,
                             adjust_completeness = opts$adjust_completeness)
    
    dat[indicator_common_id %in% EXCLUDED_FROM_ADJUSTMENT, count_working := count]
    dat <- dat[, .(facility_id, indicator_common_id, period_id, count_final = count_working)]
    setnames(dat, "count_final", paste0("count_final_", scenario))
    
    results[[scenario]] <- dat
    gc()
  }
  
  combined <- Reduce(function(x, y) merge(x, y, by = join_cols, all = TRUE), results)
  return(combined)
}

# ------------------- Main Execution ------------------------------------------------------------------------

print("Running adjustments analysis...")

# Apply adjustment scenarios
adjusted_data_final <- apply_adjustments_scenarios(
  raw_data = raw_data,
  completeness_data = completeness_data,
  outlier_data = outlier_data
)

# Create metadata lookups
geo_lookup <- unique(raw_data[, .SD, .SDcols = c("facility_id", grep("^admin_area_", names(raw_data), value = TRUE))])
period_lookup <- unique(completeness_data[, .(period_id, quarter_id, year)])

# Merge metadata into facility-level adjusted data
setDT(adjusted_data_final)
setkey(adjusted_data_final, facility_id)
setkey(geo_lookup, facility_id)
setkey(period_lookup, period_id)

adjusted_data_export <- merge(adjusted_data_final, geo_lookup, by = "facility_id", all.x = TRUE)
adjusted_data_export <- merge(adjusted_data_export, period_lookup, by = "period_id", all.x = TRUE)

# Detect geo columns
geo_cols <- grep("^admin_area_[0-9]+$", names(adjusted_data_export), value = TRUE)
geo_admin_area_sub <- setdiff(geo_cols, "admin_area_1")

message("Detected admin area columns: ", paste(geo_cols, collapse = ", "))
message("Using for subnational aggregation: ", paste(geo_admin_area_sub, collapse = ", "))

# Reorder columns for export
setcolorder(adjusted_data_export, c(
  "facility_id", 
  geo_admin_area_sub, 
  "period_id", 
  "quarter_id", 
  "year", 
  "indicator_common_id"
))

# --------------------------- Subnational Output ---------------------------
adjusted_data_admin_area_final <- adjusted_data_export[
  ,
  .(
    count_final_none        = sum(count_final_none, na.rm = TRUE),
    count_final_outliers    = sum(count_final_outliers, na.rm = TRUE),
    count_final_completeness= sum(count_final_completeness, na.rm = TRUE),
    count_final_both        = sum(count_final_both, na.rm = TRUE)
  ),
  by = c(geo_admin_area_sub, "indicator_common_id", "period_id")
]

adjusted_data_admin_area_final <- merge(
  adjusted_data_admin_area_final,
  period_lookup,
  by = "period_id",
  all.x = TRUE
)

setcolorder(adjusted_data_admin_area_final, c(
  geo_admin_area_sub, 
  "period_id", 
  "quarter_id", 
  "year", 
  "indicator_common_id"
))

# --------------------------- National Output ---------------------------
adjusted_data_national_final <- adjusted_data_export[
  ,
  .(
    count_final_none        = sum(count_final_none, na.rm = TRUE),
    count_final_outliers    = sum(count_final_outliers, na.rm = TRUE),
    count_final_completeness= sum(count_final_completeness, na.rm = TRUE),
    count_final_both        = sum(count_final_both, na.rm = TRUE)
  ),
  by = .(admin_area_1, indicator_common_id, period_id)
]

adjusted_data_national_final <- merge(
  adjusted_data_national_final,
  period_lookup,
  by = "period_id",
  all.x = TRUE
)

setcolorder(adjusted_data_national_final, c(
  "admin_area_1", 
  "period_id", 
  "quarter_id", 
  "year", 
  "indicator_common_id"
))

# --------------------------- Save Outputs ---------------------------
adjusted_data_export_clean <- adjusted_data_export[, !"admin_area_1"]

fwrite(adjusted_data_export_clean,     "M2_adjusted_data.csv",            na = "NA")
fwrite(adjusted_data_admin_area_final, "M2_adjusted_data_admin_area.csv", na = "NA")
fwrite(adjusted_data_national_final,   "M2_adjusted_data_national.csv",   na = "NA")

print("Adjustments completed and all outputs saved.")
