COUNTRY_ISO3 <- "NGA"

SELECTEDCOUNT <- "count_final_both"  #use count_final_none or count_final_completeness
VISUALIZATIONCOUNT <- "count_final_outliers" 

SMOOTH_K <- 7                          # Window size (in months) for rolling median smoothing of predicted counts.
                                       # Used in the control chart to reduce noise in trend estimation. MUST BE ODD

MADS_THRESHOLD <- 1.5                 # Threshold (in MAD units) for detecting sharp deviations in robust control chart.
                                       # If residual/MAD > THRESHOLD, the month is flagged as a sharp disruption.

DIP_THRESHOLD <- 0.90                  # Threshold for dips: a month is flagged if actual count falls below
                                       # 90% of the smoothed expected volume (i.e., a ≥10% drop).
                                       # Set to 0.80 for a more conservative detection -> to flag big drops.

DIFFPERCENT <- 10                      # Difference threshold (in percent): if the actual volume differs from the predicted
                                       # volume by more than ±10%, use the predicted value in plotting disruptions.

RUN_DISTRICT_MODEL <- TRUE             # Set to TRUE to run regressions at the lowest geographic level (admin_area_3).
                                       # Set to FALSE for faster runtime.


PROJECT_DATA_HMIS <- "hmis_NGA.csv"
#-------------------------------------------------------------------------------------------------------------
# CB - R code FASTR PROJECT
# Last edit: 2026 Sep 21
# Module: SERVICE UTILIZATION


# This script analyzes disruptions in essential health services using adjusted HMIS data (M2_adjusted_data.csv). 
# It has two main components:
#   1. Control Chart Analysis: Identifies whether deviations in service volumes are part of normal fluctuations 
#      or indicate significant disruptions.
#   2. Disruption Analysis: Quantifies the impact of these disruptions by measuring how service volumes changed 
#      during flagged periods.


# Ce code analyse les perturbations des services de santé essentiels à partir des données HMIS ajustées.
# Il comprend deux parties principales :
#   1. Analyse des cartes de contrôle : Détermine si les écarts dans les volumes de services relèvent de fluctuations normales 
#      ou signalent des perturbations importantes.
#   2. Analyse des perturbations : Quantifie l'impact des perturbations en mesurant les variations des volumes 
#      de services pendant les périodes signalées.

# ------------------------------------- KEY OUTPUTS ----------------------------------------------------------
# FILE: control_chart_results.csv       # Facility-level control chart analysis results with flags for anomalies.
# FILE: indicator_results.csv           # Indicator-level trends with control limits.
# FILE: M3_chartout.csv                 # Filtered dataset of flagged disruptions only.
# FILE: M3_disruptions_analysis.csv     # Outputs from the disruption analysis.

#-------------------------------------------------------------------------------------------------------------
# Script-internal setting (not a platform parameter; must stay below the header so the platform keeps it).
# The admin_area_4 analysis is always off. Empty AA4 compatibility files are still written further down.
RUN_ADMIN_AREA_4_ANALYSIS <- FALSE     # Set to TRUE to run finest-level analysis (admin_area_4)
                                       # Warning: This can be very slow for large datasets

# Load required libraries
library(data.table)  # For memory-efficient data operations
library(lubridate)
library(zoo)
library(MASS)    # For the rlm() function >> robust regression
library(fixest)  # For panel regressions (alternative to 'xtreg' in Stata)
library(stringr)
library(dplyr)
library(tidyr)

# Ensure dplyr::select is used (MASS::select masks it)
select <- dplyr::select

# Memory notes: the facility-level table can exceed 80 million rows. Only the columns that are used are read,
# joins are done in place, dates are integers, and the regression loops keep (unit, period) totals only:
# facility-level predictions are never stored, so memory does not grow with the number of facilities.

# Memory tracking helper - logs to console AND file
mem_usage <- function(msg) {
  # R memory usage
  gc_result <- gc()
  mem_mb <- sum(gc_result[,2])

  # System memory (if available)
  sys_mem <- tryCatch({
    if (Sys.info()["sysname"] == "Linux") {
      meminfo <- readLines("/proc/meminfo")
      total <- as.numeric(gsub("[^0-9]", "", grep("MemTotal", meminfo, value = TRUE)))
      available <- as.numeric(gsub("[^0-9]", "", grep("MemAvailable", meminfo, value = TRUE)))
      used_pct <- round((1 - available/total) * 100, 1)
      sprintf("System: %s%%", used_pct)
    } else {
      "N/A"
    }
  }, error = function(e) "N/A")

  # Log message
  log_msg <- sprintf("[MEM] %s: R=%.0f MB | %s | %s\n",
                     msg, mem_mb, sys_mem, format(Sys.time(), "%H:%M:%S"))

  # Print to console
  cat(log_msg)

  # Append to log file (survives crashes)
  cat(log_msg, file = "M3_memory_log.txt", append = TRUE)

  # Flush to ensure it's written immediately
  flush(stdout())

  invisible(mem_mb)
}

#-------------------------------------------------------------------------------------------------------------
# SAFETY: Clean up any temporary files from previous runs
#-------------------------------------------------------------------------------------------------------------
print("Checking for temporary files from previous runs...")

# Clear old memory log
if (file.exists("M3_memory_log.txt")) {
  file.remove("M3_memory_log.txt")
}
cat(sprintf("=== Memory Log Started: %s ===\n", Sys.time()), file = "M3_memory_log.txt")

old_temp_files <- list.files(pattern = "^M3_temp_.*\\.csv$")
if (length(old_temp_files) > 0) {
  print(paste("  WARNING: Found", length(old_temp_files), "old temporary files from a previous version - deleting"))
  file.remove(old_temp_files)
}

# Set CONTROL_CHART_LEVEL conditionally based on analysis flags
if (RUN_ADMIN_AREA_4_ANALYSIS) {
  CONTROL_CHART_LEVEL <- "admin_area_4"
} else if (RUN_DISTRICT_MODEL) {
  CONTROL_CHART_LEVEL <- "admin_area_3"
} else {
  CONTROL_CHART_LEVEL <- "admin_area_2"  # Default
}

#-------------------------------------------------------------------------------------------------------------
# STEP 1: CONTROL CHART ANALYSIS
#-------------------------------------------------------------------------------------------------------------

RISE_THRESHOLD <- 1 / DIP_THRESHOLD  # Threshold for rises: a month is flagged if actual count exceeds
                                     # ~111% of expected volume (i.e., a ≥10% rise). Mirrors the dip logic.

print("Loading data for control chart analysis...")
# Only the columns that are used are read
admin_area_1_lookup <- unique(fread(PROJECT_DATA_HMIS, select = c("facility_id", "admin_area_1")))
if (anyDuplicated(admin_area_1_lookup, by = "facility_id")) {
  admin_area_1_lookup <- unique(admin_area_1_lookup, by = "facility_id")
}
outlier_data <- fread("M1_output_outliers.csv",
                      select = c("facility_id", "indicator_common_id", "period_id", "outlier_flag"))
m2_cols  <- names(fread("M2_adjusted_data.csv", nrows = 0))
geo_cols <- grep("^admin_area_[0-9]+$", m2_cols, value = TRUE)
data <- fread("M2_adjusted_data.csv",
              select = unique(c("facility_id", geo_cols, "indicator_common_id", "period_id",
                                SELECTEDCOUNT, VISUALIZATIONCOUNT)))
mem_usage("After loading data")

print("Preparing data for the control chart analysis...")

# Join outlier data (in place)
data[outlier_data, outlier_flag := i.outlier_flag, on = .(facility_id, indicator_common_id, period_id)]
rm(outlier_data)
data[, outlier_flag := fifelse(is.na(outlier_flag), 0L, outlier_flag)]

# Filter outliers
data <- data[outlier_flag != 1L]
data[, outlier_flag := NULL]
invisible(gc())

# Calendar columns: period_id is YYYYMM; the date is kept as an integer (days since 1970-01-01, the same
# number a Date carries) so the models see exactly the values they used to
data[, period_id := as.integer(period_id)]
data[, `:=`(year = period_id %/% 100L, month = period_id %% 100L)]
period_dates <- unique(data[, .(period_id)])
period_dates[, date := as.integer(as.Date(sprintf("%04d-%02d-01", period_id %/% 100L, period_id %% 100L)))]
data[period_dates, date := i.date, on = "period_id"]
data[, count_model := get(SELECTEDCOUNT)]

print(paste("Aggregating data to", CONTROL_CHART_LEVEL, "level..."))
province_data <- data[, .(count_original = sum(count_model, na.rm = TRUE)),
                      by = c("indicator_common_id", CONTROL_CHART_LEVEL, "date")]
province_data[, date := as.Date(date, origin = "1970-01-01")]
data[, count_model := NULL]
province_data <- as_tibble(province_data) %>%
  arrange(indicator_common_id, !!sym(CONTROL_CHART_LEVEL), date) %>%
  group_by(indicator_common_id, !!sym(CONTROL_CHART_LEVEL)) %>%
  mutate(panelvar = cur_group_id()) %>%
  ungroup()

print("Filling missing months and metadata...")
province_data <- province_data %>%
  group_by(panelvar) %>%
  complete(date = seq(min(date, na.rm = TRUE), max(date, na.rm = TRUE), by = "month")) %>%
  fill(indicator_common_id, !!sym(CONTROL_CHART_LEVEL), .direction = "downup") %>%
  ungroup()

print("Removing months with extremely low counts...")
province_data <- province_data %>%
  group_by(panelvar) %>%
  mutate(
    globalmean = mean(count_original, na.rm = TRUE),
    count = ifelse(count_original / globalmean < 0.5, NA_real_, count_original)
  ) %>%
  ungroup()

print("Interpolating missing/removed values for modeling...")
province_data <- province_data %>%
  group_by(panelvar) %>%
  arrange(date) %>%
  mutate(
    count = zoo::na.approx(count, na.rm = FALSE, maxgap = Inf, rule = 2)
  ) %>%
  ungroup()

mem_usage("After data preparation")

print("Running robust control chart analysis for each panel...")

# Function control chart -------------------------------------------------------
robust_control_chart <- function(panel_data, selected_count) {
  panel_data <- panel_data %>%
    mutate(month_factor = factor(month(date)))
  
  # Count non-missing obs and unique dates
  n_obs <- sum(!is.na(panel_data[[selected_count]]))
  n_dates <- length(unique(panel_data$date[!is.na(panel_data[[selected_count]])]))
  
  # Model fallback logic
  if (n_obs >= 12 && n_dates > 12) {
    # Safe to use full model
    mod <- tryCatch({
      rlm(as.formula(paste(selected_count, "~ month_factor + as.numeric(date)")),
          data = panel_data, maxit = 100)
    }, error = function(e) {
      warning(paste("Full model failed, fallback to trend-only. Error:", e$message))
      NULL
    })

    # Check convergence
    if (!is.null(mod) && !mod$converged) {
      panel_id <- unique(panel_data[[CONTROL_CHART_LEVEL]])
      indicator_id <- unique(panel_data$indicator_common_id)
      print(paste("WARNING: Full model failed to converge for", CONTROL_CHART_LEVEL, "=", panel_id,
                  "| Indicator =", indicator_id))
    }

  } else if (n_obs >= 12) {
    # Use simpler model (trend only)
    mod <- tryCatch({
      rlm(as.formula(paste(selected_count, "~ as.numeric(date)")),
          data = panel_data, maxit = 100)
    }, error = function(e) {
      warning(paste("Trend-only model failed. Error:", e$message))
      NULL
    })

    # Check convergence
    if (!is.null(mod) && !mod$converged) {
      panel_id <- unique(panel_data[[CONTROL_CHART_LEVEL]])
      indicator_id <- unique(panel_data$indicator_common_id)
      print(paste("WARNING: Trend-only model failed to converge for", CONTROL_CHART_LEVEL, "=", panel_id,
                  "| Indicator =", indicator_id))
    }

  } else {
    mod <- NULL
  }
  
  # Predict or fallback to median
  # Use pmax to ensure predictions are never negative (count data cannot be negative)
  panel_data <- panel_data %>%
    mutate(count_predict = if (!is.null(mod)) {
      pmax(0, predict(mod, newdata = panel_data))
    } else {
      median(panel_data[[selected_count]], na.rm = TRUE)
    })
  
  # Smoothing
  panel_data <- panel_data %>%
    arrange(date) %>%
    mutate(
      count_smooth = zoo::rollmedian(count_predict, k = SMOOTH_K, fill = NA, align = "center"),
      count_smooth = ifelse(is.na(count_smooth), count_predict, count_smooth)
    )
  
  # Residuals and MAD-based control limits
  panel_data <- panel_data %>%
    mutate(
      residual = count_original - count_smooth,
      robust_control = residual / (mad(residual, constant = 1, na.rm = TRUE) + 1e-6),
      tag_sharp = ifelse(!is.na(robust_control) & abs(robust_control) >= MADS_THRESHOLD, 1, 0),
      mild_flag = ifelse(!is.na(robust_control) & abs(robust_control) >= 1 & abs(robust_control) < MADS_THRESHOLD, 1, 0),
      mild_cumulative = zoo::rollapply(mild_flag, width = 3, align = "right", fill = NA, FUN = sum, na.rm = TRUE),
      tag_sustained = ifelse(mild_cumulative >= 3 & abs(robust_control) >= 1.5, 1, 0),
      dip_flag = ifelse(is.na(count_original) | count_original < DIP_THRESHOLD * count_smooth, 1, 0)
    )
  
  # Dips
  dip_rle <- rle(panel_data$dip_flag)
  panel_data$tag_sustained_dip <- inverse.rle(with(dip_rle, list(
    lengths = lengths,
    values = ifelse(values == 1 & lengths >= 3, 1, 0)
  )))
  
  # Missing and rise tagging
  panel_data <- panel_data %>%
    mutate(
      is_missing = is.na(count_original) | count_original == 0,
      missing_roll = zoo::rollapply(is_missing, width = 3, align = "right", fill = NA, FUN = sum, na.rm = TRUE),
      tag_missing = ifelse(missing_roll >= 2, 1, 0),
      rise_flag = ifelse(!is.na(count_original) & count_original > RISE_THRESHOLD * count_smooth, 1, 0)
    )
  
  rise_rle <- rle(panel_data$rise_flag)
  panel_data$tag_sustained_rise <- inverse.rle(with(rise_rle, list(
    lengths = lengths,
    values = ifelse(values == 1 & lengths >= 3, 1, 0)
  )))
  
  # Final tagging
  panel_data <- panel_data %>%
    mutate(
      tagged = case_when(
        tag_sharp == 1 |
          tag_sustained == 1 |
          tag_sustained_dip == 1 |
          tag_sustained_rise == 1 |
          tag_missing == 1 ~ 1,
        TRUE ~ 0
      ),
      tagged = replace_na(tagged, 0)
    ) %>%
    group_by(!!sym(CONTROL_CHART_LEVEL)) %>%
    mutate(
      last_6_months = ifelse(date >= max(date) - months(6), 1, 0),
      tagged = ifelse(last_6_months == 1, 1, tagged)
    ) %>%
    ungroup()
  
  return(panel_data)
}

# Run for all panels (the panel table is small: one row per unit x indicator x month) -----------------------
panel_list <- unique(province_data$panelvar)
results_list <- vector("list", length(panel_list))
for (i in seq_along(panel_list)) {
  panel <- panel_list[i]
  if (i %% 100 == 0) print(paste0("Processing panel ", i, " of ", length(panel_list)))
  panel_data <- province_data %>% filter(panelvar == panel)
  results_list[[i]] <- robust_control_chart(panel_data, "count")
}
M3_chartout <- rbindlist(results_list)
rm(results_list, province_data)
invisible(gc())
print("Control chart analysis complete")
mem_usage("After control chart analysis")

#-------------------------------------------------------------------------------------------------------------
# STEP 2: DISRUPTION REGRESSION ANALYSIS
#-------------------------------------------------------------------------------------------------------------
print("Loading and preparing data for disruption analysis...")
M3_chartout_selected <- M3_chartout[, c("date", "indicator_common_id", CONTROL_CHART_LEVEL, "tagged"), with = FALSE]
rm(M3_chartout)

# Join control chart results and admin_area_1 (both in place)
chart_join <- copy(M3_chartout_selected)
chart_join[, date := as.integer(date)]
data[chart_join, tagged := i.tagged, on = c("date", "indicator_common_id", CONTROL_CHART_LEVEL)]
rm(chart_join)
data[is.na(tagged), tagged := 0]
data[admin_area_1_lookup, admin_area_1 := i.admin_area_1, on = "facility_id"]
data_disruption <- data
rm(data)
invisible(gc())
mem_usage("After disruption data prep")

# Regression helper: fits the disruption model on `dt` (rows of one indicator, and one unit for the
# subnational levels), returns the expected volume summed by `sum_cols` x date, or NULL if the model fails.
reg_formula <- as.formula(paste(SELECTEDCOUNT, "~ date + factor(month) + tagged"))
fit_expected <- function(dt, sum_cols, cluster_col = NULL) {
  n_clusters <- if (!is.null(cluster_col)) uniqueN(dt[[cluster_col]], na.rm = TRUE) else 0L
  model <- tryCatch(
    if (n_clusters > 1) {
      feols(reg_formula, data = dt, cluster = as.formula(paste0("~", cluster_col)))
    } else {
      feols(reg_formula, data = dt)
    },
    error = function(e) { NULL }
  )
  if (is.null(model) || anyNA(coef(model))) return(NULL)
  # If 'tagged' was dropped, effect is 0, otherwise get it from the model
  disruption_effect <- if ("tagged" %in% names(coef(model))) coef(model)["tagged"] else 0
  dt[, expect := pmax(0, predict(model, newdata = dt) - (tagged * disruption_effect))]
  dt[, .(count_expect_sum = sum(expect, na.rm = TRUE)), by = c(sum_cols, "date")]
}

print("Running panel regressions...")
indicators <- unique(data_disruption$indicator_common_id)
setkey(data_disruption, indicator_common_id)

# Step 4a: Indicator level --------------------------------------------------------------------------------
print("Running regressions at the indicator level...")
expect_admin1 <- list()
for (indicator in indicators) {
  print(paste("Processing:", indicator))
  indicator_data <- data_disruption[.(indicator)][!is.na(get(SELECTEDCOUNT))]
  if (nrow(indicator_data) == 0) { next }
  res <- fit_expected(indicator_data, "admin_area_1", "admin_area_3")
  if (!is.null(res)) { res[, indicator_common_id := indicator]; expect_admin1[[indicator]] <- res }
}
expect_admin1 <- rbindlist(expect_admin1)
print("Indicator-level regression complete")

# Step 4b: Indicator x Province -------------------------------------------------------------------------------
print("Running regressions at the province level...")
expect_admin2 <- list()
for (indicator in indicators) {
  indicator_data <- data_disruption[.(indicator)][!is.na(get(SELECTEDCOUNT))]
  if (nrow(indicator_data) == 0) { next }
  for (province_data in split(indicator_data, by = "admin_area_2")) {
    if (is.na(province_data$admin_area_2[1])) { next }
    print(paste("Processing:", indicator, "in region:", province_data$admin_area_2[1]))
    res <- fit_expected(province_data, "admin_area_2", "admin_area_3")
    if (!is.null(res)) {
      res[, indicator_common_id := indicator]
      expect_admin2[[length(expect_admin2) + 1L]] <- res
    }
  }
}
expect_admin2 <- rbindlist(expect_admin2)
print("Province-level regression complete")

# Step 4c: Indicator x District -------------------------------------------------------------------------------
expect_admin3 <- data.table()
if (RUN_DISTRICT_MODEL) {
  print("Running regressions at the district level...")
  mem_usage("Before district loop starts")
  expect_admin3 <- list()
  for (indicator in indicators) {
    indicator_data <- data_disruption[.(indicator)][!is.na(get(SELECTEDCOUNT)) & !is.na(tagged) & !is.na(date)]
    if (nrow(indicator_data) == 0) { next }
    # Note: a district name shared by two provinces is one model, as before; totals are kept per province x district
    for (district_data in split(indicator_data, by = "admin_area_3")) {
      if (is.na(district_data$admin_area_3[1]) || nrow(district_data) < 10) { next }
      print(paste("Processing:", indicator, "in region:", district_data$admin_area_3[1]))
      res <- fit_expected(district_data, c("admin_area_2", "admin_area_3"), "admin_area_4")
      if (!is.null(res)) {
        res[, indicator_common_id := indicator]
        expect_admin3[[length(expect_admin3) + 1L]] <- res
      }
    }
  }
  expect_admin3 <- rbindlist(expect_admin3)
  mem_usage("After district loop")
  print("District/State-level regression complete")
}

# Step 4d: Indicator x admin area 4 ---------------------------------------------------------------------------
expect_admin4 <- data.table()
if (RUN_ADMIN_AREA_4_ANALYSIS) {
  print("Running regressions at the admin area 4 level...")
  expect_admin4 <- list()
  for (indicator in indicators) {
    indicator_data <- data_disruption[.(indicator)][!is.na(get(SELECTEDCOUNT)) & !is.na(tagged) & !is.na(date)]
    if (nrow(indicator_data) == 0) { next }
    for (admin_unit_data in split(indicator_data, by = "admin_area_4")) {
      if (is.na(admin_unit_data$admin_area_4[1]) || nrow(admin_unit_data) < 8) { next }
      print(paste("Processing:", indicator, "in admin_area_4:", admin_unit_data$admin_area_4[1]))
      res <- fit_expected(admin_unit_data, c("admin_area_2", "admin_area_3", "admin_area_4"), NULL)
      if (!is.null(res)) {
        res[, indicator_common_id := indicator]
        expect_admin4[[length(expect_admin4) + 1L]] <- res
      }
    }
  }
  expect_admin4 <- rbindlist(expect_admin4)
  print("Admin_area_4-level regression complete")
}
mem_usage("After all regressions complete")

#-------------------------------------------------------------------------------
# STEP 3: PREPARE RESULTS FOR VISUALIZATION
#-------------------------------------------------------------------------------
# Observed totals per unit x period x indicator over every row (as before), then the model totals attached;
# units without a successful model get an expected sum of 0, exactly as the previous per-row join produced.
summarise_level <- function(unit_cols, expect_dt) {
  obs <- data_disruption[, .(count_sum = sum(get(VISUALIZATIONCOUNT), na.rm = TRUE)),
                         by = c(unit_cols, "period_id", "indicator_common_id")]
  if (nrow(expect_dt) > 0) {
    expect_dt <- copy(expect_dt)
    expect_dt[period_dates, period_id := i.period_id, on = "date"]
    expect_dt <- expect_dt[, .(count_expect_sum = sum(count_expect_sum)), by = c(unit_cols, "period_id", "indicator_common_id")]
    obs[expect_dt, count_expect_sum := i.count_expect_sum, on = c(unit_cols, "period_id", "indicator_common_id")]
    obs[is.na(count_expect_sum), count_expect_sum := 0]
  } else {
    obs[, count_expect_sum := 0]
  }
  obs[, count_expected_if_above_diff_threshold := ifelse(
    abs(100 * (count_expect_sum - count_sum) / count_expect_sum) > DIFFPERCENT,
    count_expect_sum, count_sum)]
  setorderv(obs, c(unit_cols, "period_id", "indicator_common_id"))
  as.data.frame(obs)
}

key_messages <- function(summary_df, unit_cols) {
  summary_df %>%
    mutate(
      shortfall_absolute = pmax(0, count_expect_sum - count_sum, na.rm = TRUE),
      shortfall_percent = ifelse(count_expect_sum > 0,
                                 (count_expect_sum - count_sum) / count_expect_sum * 100, 0),
      surplus_absolute = pmax(0, count_sum - count_expect_sum, na.rm = TRUE),
      surplus_percent = ifelse(count_expect_sum > 0,
                               (count_sum - count_expect_sum) / count_expect_sum * 100, 0)
    ) %>%
    select(all_of(unit_cols), indicator_common_id, period_id,
           count_sum, count_expect_sum,
           shortfall_absolute, shortfall_percent,
           surplus_absolute, surplus_percent) %>%
    arrange(indicator_common_id, period_id)
}

disruption_export <- function(summary_df, unit_cols) {
  summary_df %>% select(all_of(unit_cols), indicator_common_id, period_id,
                        count_sum, count_expect_sum, count_expected_if_above_diff_threshold)
}

empty_key_messages <- function(unit_cols) {
  as.data.frame(c(setNames(replicate(length(unit_cols), character(0), simplify = FALSE), unit_cols),
                  list(indicator_common_id = character(0), period_id = integer(0),
                       count_sum = numeric(0), count_expect_sum = numeric(0),
                       shortfall_absolute = numeric(0), shortfall_percent = numeric(0),
                       surplus_absolute = numeric(0), surplus_percent = numeric(0))))
}
empty_disruption <- function(unit_cols) {
  as.data.frame(c(setNames(replicate(length(unit_cols), character(0), simplify = FALSE), unit_cols),
                  list(indicator_common_id = character(0), period_id = integer(0),
                       count_sum = numeric(0), count_expect_sum = numeric(0),
                       count_expected_if_above_diff_threshold = numeric(0))))
}

print("Creating summary disruptions at national level...")
summary_disruption_admin1 <- summarise_level("admin_area_1", expect_admin1)
print("Saving disruptions data for external analysis...")
write.csv(key_messages(summary_disruption_admin1, "admin_area_1"),
          "M3_all_indicators_shortfalls_admin_area_1.csv", row.names = FALSE)

# Admin_area_2 level: produced when at least one province model succeeded (as before)
if (nrow(expect_admin2) > 0) {
  print("Creating summary disruptions at admin_area_2 level...")
  summary_disruption_admin2 <- summarise_level("admin_area_2", expect_admin2)
  write.csv(key_messages(summary_disruption_admin2, "admin_area_2"),
            "M3_all_indicators_shortfalls_admin_area_2.csv", row.names = FALSE)
} else {
  print("No admin_area_2 data available - creating empty key messages file for compatibility...")
  write.csv(empty_key_messages("admin_area_2"), "M3_all_indicators_shortfalls_admin_area_2.csv", row.names = FALSE)
}

# Admin_area_3 level
if (RUN_DISTRICT_MODEL && nrow(expect_admin3) > 0) {
  print("Creating summary disruptions at admin_area_3 level...")
  summary_disruption_admin3 <- summarise_level(c("admin_area_2", "admin_area_3"), expect_admin3)
  write.csv(key_messages(summary_disruption_admin3, c("admin_area_2", "admin_area_3")),
            "M3_all_indicators_shortfalls_admin_area_3.csv", row.names = FALSE)
} else {
  print("No admin_area_3 data available or RUN_DISTRICT_MODEL=FALSE - creating empty key messages file for compatibility...")
  write.csv(empty_key_messages(c("admin_area_2", "admin_area_3")), "M3_all_indicators_shortfalls_admin_area_3.csv", row.names = FALSE)
}

# Admin_area_4 level
if (RUN_ADMIN_AREA_4_ANALYSIS && nrow(expect_admin4) > 0) {
  print("Creating summary disruptions at admin_area_4 level...")
  summary_disruption_admin4 <- summarise_level(c("admin_area_2", "admin_area_3", "admin_area_4"), expect_admin4)
  write.csv(key_messages(summary_disruption_admin4, c("admin_area_2", "admin_area_3", "admin_area_4")),
            "M3_all_indicators_shortfalls_admin_area_4.csv", row.names = FALSE)
} else {
  print("No admin_area_4 data available or RUN_ADMIN_AREA_4_ANALYSIS=FALSE - creating empty key messages file for compatibility...")
  write.csv(empty_key_messages(c("admin_area_2", "admin_area_3", "admin_area_4")), "M3_all_indicators_shortfalls_admin_area_4.csv", row.names = FALSE)
}

print("=== KEY MESSAGES DATASETS EXPORT COMPLETE ===")

rm(data_disruption)
invisible(gc())

# Save Result Objects ----------------------------------------------------------
print("Saving results...")

# Adjusted data pass-through for viz: a straight file copy (no read/write of the largest file)
print("Writing service utilization output...")
file.copy("M2_adjusted_data.csv", "M3_service_utilization.csv", overwrite = TRUE)

# Export control chart results (period_id only)
print("Saving control chart results...")
M3_chartout_export <- M3_chartout_selected %>%
  mutate(period_id = as.integer(format(date, "%Y%m"))) %>%
  select(
    !!sym(CONTROL_CHART_LEVEL),  # dynamic column
    indicator_common_id,
    period_id,
    tagged
  )
write.csv(M3_chartout_export, "M3_chartout.csv", row.names = FALSE)

# Export summary disruptions
print("Saving national level disruption analysis...")
write.csv(disruption_export(summary_disruption_admin1, "admin_area_1"),
          "M3_disruptions_analysis_admin_area_1.csv", row.names = FALSE)

if (exists("summary_disruption_admin2")) {
  print("Saving regional level disruption analysis...")
  write.csv(disruption_export(summary_disruption_admin2, "admin_area_2"),
            "M3_disruptions_analysis_admin_area_2.csv", row.names = FALSE)
} else {
  print("Creating empty admin_area_2 file for compatibility...")
  write.csv(empty_disruption("admin_area_2"), "M3_disruptions_analysis_admin_area_2.csv", row.names = FALSE)
}

if (RUN_DISTRICT_MODEL && exists("summary_disruption_admin3")) {
  print("Saving state level disruption analysis...")
  write.csv(disruption_export(summary_disruption_admin3, c("admin_area_2", "admin_area_3")),
            "M3_disruptions_analysis_admin_area_3.csv", row.names = FALSE)
} else {
  print("Creating empty admin_area_3 file for compatibility...")
  write.csv(empty_disruption(c("admin_area_2", "admin_area_3")), "M3_disruptions_analysis_admin_area_3.csv", row.names = FALSE)
}

if (RUN_ADMIN_AREA_4_ANALYSIS && exists("summary_disruption_admin4")) {
  print("Saving district level disruption analysis...")
  write.csv(disruption_export(summary_disruption_admin4, c("admin_area_2", "admin_area_3", "admin_area_4")),
            "M3_disruptions_analysis_admin_area_4.csv", row.names = FALSE)
} else {
  print("Creating empty admin_area_4 file for compatibility...")
  write.csv(empty_disruption(c("admin_area_2", "admin_area_3", "admin_area_4")), "M3_disruptions_analysis_admin_area_4.csv", row.names = FALSE)
}

print("=== ANALYSIS COMPLETE ===")
mem_usage("Final")
print("=== MODULE 3 COMPLETE ===")
