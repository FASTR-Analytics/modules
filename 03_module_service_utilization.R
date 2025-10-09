COUNTRY_ISO3 <- "SOM"

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

RUN_DISTRICT_MODEL <- FALSE             # Set to TRUE to run regressions at the lowest geographic level (admin_area_3).
                                       # Set to FALSE for faster runtime.

RUN_ADMIN_AREA_4_ANALYSIS <- FALSE     # Set to TRUE to run finest-level analysis (admin_area_4)
                                       # Warning: This can be very slow for large datasets


PROJECT_DATA_HMIS <- "hmis_sierraleone.csv"
#-------------------------------------------------------------------------------------------------------------
# CB - R code FASTR PROJECT
# Last edit: 2025 Oct 9
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
# Load required libraries
library(lubridate)
library(zoo)
library(MASS)    # For the rlm() function >> robust regression
library(fixest)  # For panel regressions (alternative to 'xtreg' in Stata)
library(stringr)
library(dplyr)
library(tidyr)

# Ensure dplyr::select is used (MASS::select masks it)
select <- dplyr::select

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
raw_data <- read.csv(PROJECT_DATA_HMIS)
outlier_data <- read.csv("M1_output_outliers.csv")
data <- read.csv("M2_adjusted_data.csv")


admin_area_1_lookup <- raw_data %>%
  distinct(facility_id, admin_area_1)
rm(raw_data)
gc()

print("Preparing data for the control chart analysis...")

data <- data %>%
  left_join(
    outlier_data[, c("facility_id", "indicator_common_id", "period_id", "outlier_flag")],
    by = c("facility_id", "indicator_common_id", "period_id")
  ) %>%
  mutate(
    outlier_flag = coalesce(outlier_flag, 0L),
    year = as.integer(substr(period_id, 1, 4)),
    month = as.integer(substr(period_id, 5, 6)),
    date = as.Date(sprintf("%04d-%02d-01", year, month)),
    is_missing = is.na(count_final_none),
    count_model = .data[[SELECTEDCOUNT]]
  ) %>%
  filter(outlier_flag != 1)
rm(outlier_data)
gc()

print("Assigning panel IDs...")
data <- data %>%
  group_by(indicator_common_id, !!sym(CONTROL_CHART_LEVEL)) %>%
  mutate(panelvar = cur_group_id()) %>%
  ungroup()

print(paste("Aggregating data to", CONTROL_CHART_LEVEL, "level..."))
province_data <- data %>%
  group_by(indicator_common_id, !!sym(CONTROL_CHART_LEVEL), date) %>%
  summarise(count_model = sum(count_model, na.rm = TRUE), .groups = "drop") %>%
  group_by(indicator_common_id, !!sym(CONTROL_CHART_LEVEL)) %>%
  mutate(panelvar = cur_group_id()) %>%
  ungroup() %>%
  rename(count_original = count_model)

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
  panel_data <- panel_data %>%
    mutate(count_predict = if (!is.null(mod)) {
      predict(mod, newdata = panel_data)
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

# Run for all panels -----------------------------------------------------------
panel_list <- unique(province_data$panelvar)

results_list <- list()
for(panel in panel_list) {
  print(paste0("Processing panel ID: ", panel))
  panel_data <- province_data %>% filter(panelvar == panel)
  panel_results <- robust_control_chart(panel_data, "count")
  results_list[[as.character(panel)]] <- panel_results
}

M3_chartout <- bind_rows(results_list)

rm(results_list)
gc()
print("Control chart analysis complete")

#-------------------------------------------------------------------------------------------------------------
# STEP 2: DISRUPTION REGRESSION ANALYSIS (MEMORY-OPTIMIZED)
#-------------------------------------------------------------------------------------------------------------
# Step 1: Load and Prepare Data (Your original code, unchanged)
print("Loading and preparing data for disruption analysis...")
M3_chartout_selected <- M3_chartout %>%
  select(date, indicator_common_id, !!sym(CONTROL_CHART_LEVEL), tagged)

rm(M3_chartout)
data_disruption <- data %>%
  left_join(M3_chartout_selected, by = c("date", "indicator_common_id", CONTROL_CHART_LEVEL)) %>%
  mutate(tagged = replace_na(tagged, 0))

# Step 2: Run Panel Regressions (Your original code, unchanged)
print("Running panel regressions...")
indicators <- unique(data_disruption$indicator_common_id)
districts <- unique(data_disruption$admin_area_3)
provinces <- unique(data_disruption$admin_area_2)


# Step 4a: Run Regression for Each Indicator -----------------------------------
print("Running regressions at the indicator level...")
# Initialize results list
indicator_results_list <- list()

for (indicator in indicators) {
  print(paste("Processing:", indicator))
  indicator_data <- data_disruption %>%
    filter(indicator_common_id == indicator) %>%
    drop_na(!!sym(SELECTEDCOUNT))

  if (nrow(indicator_data) == 0) { next }

  # Determine whether clustering is valid
  n_clusters <- indicator_data %>% pull(admin_area_3) %>% n_distinct(na.rm = TRUE)

  # Fit model: clustered if >1 cluster, unclustered otherwise
  model <- tryCatch(
    if (n_clusters > 1) {
      feols(as.formula(paste(SELECTEDCOUNT, "~ date + factor(month) + tagged")),
            data = indicator_data, cluster = ~admin_area_3)
    } else {
      feols(as.formula(paste(SELECTEDCOUNT, "~ date + factor(month) + tagged")),
            data = indicator_data)
    },
    error = function(e) { NULL }
  )
  
  if (!is.null(model) && !anyNA(coef(model))) {
    # --- ROBUSTNESS CHECK ---
    # Check if 'tagged' was dropped. If so, effect is 0, otherwise get it from the model.
    disruption_effect <- if ("tagged" %in% names(coef(model))) coef(model)["tagged"] else 0
    
    # Calculate predictions, removing the disruption effect
    indicator_data <- indicator_data %>%
      mutate(expect_admin_area_1 = predict(model, newdata = .) - (tagged * disruption_effect))
    
    # Calculate coefficients
    indicator_data <- indicator_data %>%
      group_by(date) %>%
      mutate(
        diff = expect_admin_area_1 - !!sym(SELECTEDCOUNT),
        diffmean = mean(diff, na.rm = TRUE),
        predictmean = mean(expect_admin_area_1, na.rm = TRUE),
        b_admin_area_1 = ifelse(tagged == 1 & predictmean != 0, -1 * (diffmean / predictmean), NA_real_),
        b_trend_admin_area_1 = coef(model)["date"]
      ) %>%
      ungroup()
    
    # Calculate p-value ONLY if the coefficient exists
    if ("tagged" %in% names(coef(model))) {
      tagged_se <- sqrt(diag(vcov(model)))["tagged"]
      if(!is.na(tagged_se) && tagged_se > 0) {
        p_val_1 <- 2 * pt(abs(disruption_effect / tagged_se), df.residual(model), lower.tail = FALSE)
        indicator_data$p_admin_area_1 <- p_val_1
      } else {
        indicator_data$p_admin_area_1 <- NA_real_
      }
    } else {
      # If 'tagged' was dropped, p-value is not applicable
      indicator_data$p_admin_area_1 <- NA_real_
    }
    
    # Memory-optimized save to list
    indicator_results_list[[indicator]] <- indicator_data %>%
      select(facility_id, date, indicator_common_id,
             expect_admin_area_1, b_admin_area_1,
             b_trend_admin_area_1, p_admin_area_1)
  }
}

# Combine and merge results
indicator_results_long <- bind_rows(indicator_results_list)
data_disruption <- data_disruption %>%
  left_join(indicator_results_long, by = c("facility_id", "date", "indicator_common_id"))
rm(indicator_results_list, indicator_results_long); gc()
print("Indicator-level regression complete.")


# Step 4b: Run Regression for Each Indicator × Province ------------------------
print("Running regressions at the province level...")
province_results_list <- list()

for (indicator in indicators) {
  for (province in provinces) {
    print(paste("Processing:", indicator, "in region:", province))
    province_data <- data_disruption %>%
      filter(indicator_common_id == indicator, admin_area_2 == province) %>%
      drop_na(!!sym(SELECTEDCOUNT))

    if (nrow(province_data) < 1) { next }

    # Determine whether clustering is valid
    n_clusters <- province_data %>% pull(admin_area_3) %>% n_distinct(na.rm = TRUE)

    # Build regression formula
    reg_formula <- as.formula(paste(SELECTEDCOUNT, "~ date + factor(month) + tagged"))

    # Fit model: clustered if >1 cluster, unclustered otherwise
    model_province <- tryCatch(
      if (n_clusters > 1) {
        feols(reg_formula, data = province_data, cluster = ~admin_area_3)
      } else {
        feols(reg_formula, data = province_data)
      },
      error = function(e) { NULL }
    )
    
    if (!is.null(model_province) && !anyNA(coef(model_province))) {
      disruption_effect <- if ("tagged" %in% names(coef(model_province))) coef(model_province)["tagged"] else 0
      
      province_data <- province_data %>%
        mutate(expect_admin_area_2 = predict(model_province, newdata = .) - (tagged * disruption_effect))
      
      province_data <- province_data %>%
        group_by(date) %>%
        mutate(
          b_admin_area_2 = ifelse(tagged == 1 & mean(expect_admin_area_2, na.rm=T) != 0,
                                  -1 * (mean(expect_admin_area_2 - !!sym(SELECTEDCOUNT), na.rm=T) / mean(expect_admin_area_2, na.rm=T)),
                                  NA_real_),
          b_trend_admin_area_2 = coef(model_province)["date"]
        ) %>%
        ungroup()
      
      if ("tagged" %in% names(coef(model_province))) {
        tagged_se <- sqrt(diag(vcov(model_province)))["tagged"]
        if(!is.na(tagged_se) && tagged_se > 0) {
          p_val <- 2 * pt(abs(disruption_effect / tagged_se), df.residual(model_province), lower.tail = FALSE)
          province_data$p_admin_area_2 <- p_val
        } else { province_data$p_admin_area_2 <- NA_real_ }
      } else { province_data$p_admin_area_2 <- NA_real_ }
      
      province_results_list[[paste(indicator, province, sep = "_")]] <- province_data %>%
        select(facility_id, date, indicator_common_id,
               expect_admin_area_2, b_admin_area_2,
               p_admin_area_2, b_trend_admin_area_2)
    }
  }
}

if (length(province_results_list) > 0) {
  province_results_long <- bind_rows(province_results_list)
  data_disruption <- data_disruption %>%
    left_join(province_results_long, by = c("facility_id", "date", "indicator_common_id"))
  rm(province_results_long)
}
rm(province_results_list); gc()
print("Province-level regression complete.")



# Step 4c: Run Regression for Each Indicator × District ------------------------
if (RUN_DISTRICT_MODEL) {
  print("Running regressions at the district level...")
  district_results_list <- list()
  
  for (indicator in indicators) {
    for (district in districts) {
      print(paste("Processing:", indicator, "in region:", district))
      district_data <- data_disruption %>%
        filter(indicator_common_id == indicator, admin_area_3 == district) %>%
        drop_na(!!sym(SELECTEDCOUNT), tagged, date)
      
      if (nrow(district_data) < 10) { next }
      
      model_district <- tryCatch(
        feols(as.formula(paste(SELECTEDCOUNT, "~ date + tagged")), data = district_data, cluster = ~admin_area_4),
        error = function(e) { NULL }
      )
      
      if (!is.null(model_district) && !anyNA(coef(model_district))) {
        disruption_effect <- if ("tagged" %in% names(coef(model_district))) coef(model_district)["tagged"] else 0
        
        district_data <- district_data %>%
          mutate(expect_admin_area_3 = predict(model_district, newdata = .) - (tagged * disruption_effect))
        
        district_data <- district_data %>%
          group_by(date) %>%
          mutate(
            b_admin_area_3 = ifelse(tagged == 1 & mean(expect_admin_area_3, na.rm = TRUE) != 0,
                                    -1 * (mean(expect_admin_area_3 - !!sym(SELECTEDCOUNT), na.rm = TRUE) / mean(expect_admin_area_3, na.rm = TRUE)),
                                    NA_real_)
          ) %>%
          ungroup()
        
        if ("tagged" %in% names(coef(model_district))) {
          tagged_se <- sqrt(diag(vcov(model_district)))["tagged"]
          if (!is.na(tagged_se) && tagged_se > 0) {
            p_val <- 2 * pt(abs(disruption_effect / tagged_se), df.residual(model_district), lower.tail = FALSE)
            district_data$p_admin_area_3 <- p_val
          } else { district_data$p_admin_area_3 <- NA_real_ }
        } else { district_data$p_admin_area_3 <- NA_real_ }
        
        district_results_list[[paste(indicator, district, sep = "_")]] <- district_data %>%
          select(facility_id, date, indicator_common_id, admin_area_3,
                 expect_admin_area_3, b_admin_area_3, p_admin_area_3)
      }
    }
  }
  
  if (length(district_results_list) > 0) {
    district_results_long <- bind_rows(district_results_list)
    data_disruption <- data_disruption %>%
      left_join(district_results_long, by = c("facility_id", "date", "indicator_common_id", "admin_area_3"))
    rm(district_results_long)
  }
  rm(district_results_list); gc()
  print("District/State-level regression complete.")
}

# Step 4d: Run Regression for Each Indicator × admin area 4 --------------------
if (RUN_ADMIN_AREA_4_ANALYSIS) {
  print("Running regressions at the admin area 4 level...")
  admin_area_4_units <- unique(data_disruption$admin_area_4)
  admin_area_4_results_list <- list()
  
  for (indicator in indicators) {
    for (admin_unit in admin_area_4_units) {
      print(paste("Processing:", indicator, "in admin_area_4:", admin_unit))
      admin_unit_data <- data_disruption %>%
        filter(indicator_common_id == indicator, admin_area_4 == admin_unit) %>%
        drop_na(!!sym(SELECTEDCOUNT), tagged, date)
      
      if (nrow(admin_unit_data) < 8) { next }
      
      model_admin_area_4 <- tryCatch(
        feols(as.formula(paste(SELECTEDCOUNT, "~ date + tagged")), data = admin_unit_data),
        error = function(e) { NULL }
      )
      
      if (!is.null(model_admin_area_4) && !anyNA(coef(model_admin_area_4))) {
        disruption_effect <- if ("tagged" %in% names(coef(model_admin_area_4))) coef(model_admin_area_4)["tagged"] else 0
        
        admin_unit_data <- admin_unit_data %>%
          mutate(expect_admin_area_4 = predict(model_admin_area_4, newdata = .) - (tagged * disruption_effect))
        
        admin_unit_data <- admin_unit_data %>%
          group_by(date) %>%
          mutate(
            b_admin_area_4 = ifelse(tagged == 1 & mean(expect_admin_area_4, na.rm=T) != 0,
                                    -1 * (mean(expect_admin_area_4 - !!sym(SELECTEDCOUNT), na.rm=T) / mean(expect_admin_area_4, na.rm=T)),
                                    NA_real_),
            b_trend_admin_area_4 = coef(model_admin_area_4)["date"]
          ) %>%
          ungroup()
        
        if ("tagged" %in% names(coef(model_admin_area_4))) {
          tagged_se <- sqrt(diag(vcov(model_admin_area_4)))["tagged"]
          if (!is.na(tagged_se) && tagged_se > 0) {
            p_val <- 2 * pt(abs(disruption_effect / tagged_se), df.residual(model_admin_area_4), lower.tail = FALSE)
            admin_unit_data$p_admin_area_4 <- p_val
          } else { admin_unit_data$p_admin_area_4 <- NA_real_ }
        } else { admin_unit_data$p_admin_area_4 <- NA_real_ }
        
        admin_area_4_results_list[[paste(indicator, admin_unit, sep = "_")]] <- admin_unit_data %>%
          select(facility_id, date, indicator_common_id, admin_area_4,
                 expect_admin_area_4, b_admin_area_4, p_admin_area_4, b_trend_admin_area_4)
      }
    }
  }
  
  if (length(admin_area_4_results_list) > 0) {
    admin_area_4_results_long <- bind_rows(admin_area_4_results_list)
    data_disruption <- data_disruption %>%
      left_join(admin_area_4_results_long, by = c("facility_id", "date", "indicator_common_id", "admin_area_4"))
    rm(admin_area_4_results_long)
  }
  rm(admin_area_4_results_list); gc()
  print("Admin_area_4-level regression complete.")
}

#-------------------------------------------------------------------------------
# STEP 3: PREPARE RESULTS FOR VISUALIZATION
#-------------------------------------------------------------------------------

data_disruption <- data_disruption %>%
  left_join(admin_area_1_lookup, by = "facility_id") %>%
  mutate(period_id = as.integer(format(as.Date(date), "%Y%m")))


print("Creating summary disruptions at national level...")
summary_disruption_admin1 <- data_disruption %>%
  group_by(admin_area_1, period_id, indicator_common_id) %>%
  summarise(
    count_original     = mean(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE),
    count_expect       = mean(expect_admin_area_1, na.rm = TRUE),
    b                  = mean(b_admin_area_1, na.rm = TRUE),
    p                  = mean(p_admin_area_1, na.rm = TRUE),
    num_obs            = n(),
    count_expect_sum   = sum(expect_admin_area_1, na.rm = TRUE),
    count_sum          = sum(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE),
    diff_percent       = 100 * (mean(expect_admin_area_1, na.rm = TRUE) - mean(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE)) /
      mean(expect_admin_area_1, na.rm = TRUE),
    diff_percent_sum   = 100 * (count_expect_sum - count_sum) / count_expect_sum,
    count_expected_if_above_diff_threshold = ifelse(
      abs(100 * (count_expect_sum - count_sum) / count_expect_sum) > DIFFPERCENT,
      count_expect_sum, count_sum),
    count_expect_diff = ifelse(
      abs(100 * (mean(expect_admin_area_1, na.rm = TRUE) - mean(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE)) /
            mean(expect_admin_area_1, na.rm = TRUE)) > DIFFPERCENT,
      mean(expect_admin_area_1, na.rm = TRUE), mean(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE)
    ),
    .groups = "drop"
  )

gc()


# Export data for external key message calculation - NATIONAL LEVEL
print("Saving disruptions data for external analysis...")

key_messages_dataset_admin1 <- summary_disruption_admin1 %>%
  mutate(
    year = as.integer(period_id %/% 100),
    month = period_id %% 100,
    shortfall_absolute = pmax(0, count_expect_sum - count_sum, na.rm = TRUE),
    shortfall_percent = ifelse(count_expect_sum > 0, 
                               (count_expect_sum - count_sum) / count_expect_sum * 100, 0),
    surplus_absolute = pmax(0, count_sum - count_expect_sum, na.rm = TRUE),
    surplus_percent = ifelse(count_expect_sum > 0,
                             (count_sum - count_expect_sum) / count_expect_sum * 100, 0)
  ) %>%
  select(admin_area_1, indicator_common_id, period_id,
         count_sum, count_expect_sum, 
         shortfall_absolute, shortfall_percent,
         surplus_absolute, surplus_percent) %>%
  arrange(indicator_common_id, period_id)

write.csv(key_messages_dataset_admin1, "M3_all_indicators_shortfalls_admin_area_1.csv", row.names = FALSE)


# CREATE SUMMARY OBJECTS FIRST (before exporting)
# Admin_area_2 level
if ("expect_admin_area_2" %in% names(data_disruption)) {
  print("Creating summary disruptions at admin_area_2 level...")

  summary_disruption_admin2 <- data_disruption %>%
    group_by(admin_area_2, period_id, indicator_common_id) %>%
    summarise(
      count_original     = mean(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE),
      count_expect       = mean(expect_admin_area_2, na.rm = TRUE),
      b                  = mean(b_admin_area_2, na.rm = TRUE),
      p                  = mean(p_admin_area_2, na.rm = TRUE),
      num_obs            = n(),
      count_expect_sum   = sum(expect_admin_area_2, na.rm = TRUE),
      count_sum          = sum(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE),
      diff_percent       = 100 * (mean(expect_admin_area_2, na.rm = TRUE) - mean(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE)) /
        mean(expect_admin_area_2, na.rm = TRUE),
      diff_percent_sum   = 100 * (count_expect_sum - count_sum) / count_expect_sum,
      count_expected_if_above_diff_threshold = ifelse(
        abs(100 * (count_expect_sum - count_sum) / count_expect_sum) > DIFFPERCENT,
        count_expect_sum, count_sum),
      count_expect_diff = ifelse(
        abs(100 * (mean(expect_admin_area_2, na.rm = TRUE) - mean(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE)) /
              mean(expect_admin_area_2, na.rm = TRUE)) > DIFFPERCENT,
        mean(expect_admin_area_2, na.rm = TRUE), mean(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE)
      ),
      .groups = "drop"
    )
} else {
  print("Warning: expect_admin_area_2 columns not found - will create empty admin_area_2 file")
}

# Admin_area_3 level
if (RUN_DISTRICT_MODEL & "expect_admin_area_3" %in% names(data_disruption)) {
    print("Creating summary disruptions at admin_area_3 level...")

    summary_disruption_admin3 <- data_disruption %>%
      group_by(admin_area_3, period_id, indicator_common_id) %>%
      summarise(
        count_original     = mean(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE),
        count_expect       = mean(expect_admin_area_3, na.rm = TRUE),
        b                  = mean(b_admin_area_3, na.rm = TRUE),
        p                  = mean(p_admin_area_3, na.rm = TRUE),
        num_obs            = n(),
        count_expect_sum   = sum(expect_admin_area_3, na.rm = TRUE),
        count_sum          = sum(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE),
        diff_percent       = 100 * (mean(expect_admin_area_3, na.rm = TRUE) - mean(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE)) /
          mean(expect_admin_area_3, na.rm = TRUE),
        diff_percent_sum   = 100 * (count_expect_sum - count_sum) / count_expect_sum,
        count_expected_if_above_diff_threshold = ifelse(
          abs(100 * (count_expect_sum - count_sum) / count_expect_sum) > DIFFPERCENT,
          count_expect_sum, count_sum),
        count_expect_diff = ifelse(
          abs(100 * (mean(expect_admin_area_3, na.rm = TRUE) - mean(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE)) /
                mean(expect_admin_area_3, na.rm = TRUE)) > DIFFPERCENT,
          mean(expect_admin_area_3, na.rm = TRUE), mean(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE)
        ),
        .groups = "drop"
      )
  } else {
    print("Warning: expect_admin_area_3 columns not found - will create empty admin_area_3 file")
  }

# Admin_area_4 level
if ("expect_admin_area_4" %in% names(data_disruption)) {
    print("Creating summary disruptions at admin_area_4 level...")
    summary_disruption_admin4 <- data_disruption %>%
      group_by(admin_area_4, period_id, indicator_common_id) %>%
      summarise(
        count_original     = mean(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE),
        count_expect       = mean(expect_admin_area_4, na.rm = TRUE),
        b                  = mean(b_admin_area_4, na.rm = TRUE),
        p                  = mean(p_admin_area_4, na.rm = TRUE),
        num_obs            = n(),
        count_expect_sum   = sum(expect_admin_area_4, na.rm = TRUE),
        count_sum          = sum(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE),
        diff_percent       = 100 * (mean(expect_admin_area_4, na.rm = TRUE) - mean(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE)) /
          mean(expect_admin_area_4, na.rm = TRUE),
        diff_percent_sum   = 100 * (count_expect_sum - count_sum) / count_expect_sum,
        count_expected_if_above_diff_threshold = ifelse(
          abs(100 * (count_expect_sum - count_sum) / count_expect_sum) > DIFFPERCENT,
          count_expect_sum, count_sum),
        count_expect_diff = ifelse(
          abs(100 * (mean(expect_admin_area_4, na.rm = TRUE) - mean(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE)) /
                mean(expect_admin_area_4, na.rm = TRUE)) > DIFFPERCENT,
          mean(expect_admin_area_4, na.rm = TRUE), mean(!!sym(VISUALIZATIONCOUNT), na.rm = TRUE)
        ),
        .groups = "drop"
      )
  } else {
    print("Warning: expect_admin_area_4 columns not found - will create empty admin_area_4 file")
  }

# NOW EXPORT THE KEY MESSAGES DATASETS
# Export data for external key message calculation - ADMIN_AREA_2 LEVEL
if (exists("summary_disruption_admin2") && nrow(summary_disruption_admin2) > 0) {
  print("Saving admin area 2 level key messages dataset...")
  
  key_messages_dataset_admin2 <- summary_disruption_admin2 %>%
    mutate(
      year = as.integer(period_id %/% 100),
      month = period_id %% 100,
      shortfall_absolute = pmax(0, count_expect_sum - count_sum, na.rm = TRUE),
      shortfall_percent = ifelse(count_expect_sum > 0, 
                                 (count_expect_sum - count_sum) / count_expect_sum * 100, 0),
      surplus_absolute = pmax(0, count_sum - count_expect_sum, na.rm = TRUE),
      surplus_percent = ifelse(count_expect_sum > 0,
                               (count_sum - count_expect_sum) / count_expect_sum * 100, 0)
    ) %>%
    select(admin_area_2, indicator_common_id, period_id,
           count_sum, count_expect_sum, 
           shortfall_absolute, shortfall_percent,
           surplus_absolute, surplus_percent) %>%
    arrange(indicator_common_id, period_id)
  
  write.csv(key_messages_dataset_admin2, "M3_all_indicators_shortfalls_admin_area_2.csv", row.names = FALSE)
  
} else {
  print("No admin_area_2 data available - creating empty key messages file for compatibility...")
  dummy_key_messages_admin2 <- data.frame(
    admin_area_2 = character(0),
    indicator_common_id = character(0),
    period_id = integer(0),
    count_sum = numeric(0),
    count_expect_sum = numeric(0),
    shortfall_absolute = numeric(0),
    shortfall_percent = numeric(0),
    surplus_absolute = numeric(0),
    surplus_percent = numeric(0)
  )
  write.csv(dummy_key_messages_admin2, "M3_all_indicators_shortfalls_admin_area_2.csv", row.names = FALSE)
  
}

# Export data for external key message calculation - ADMIN_AREA_3 LEVEL
if (RUN_DISTRICT_MODEL && exists("summary_disruption_admin3") && nrow(summary_disruption_admin3) > 0) {
  print("Saving admin area 3 level key messages dataset...")
  
  key_messages_dataset_admin3 <- summary_disruption_admin3 %>%
    mutate(
      year = as.integer(period_id %/% 100),
      month = period_id %% 100,
      shortfall_absolute = pmax(0, count_expect_sum - count_sum, na.rm = TRUE),
      shortfall_percent = ifelse(count_expect_sum > 0, 
                                 (count_expect_sum - count_sum) / count_expect_sum * 100, 0),
      surplus_absolute = pmax(0, count_sum - count_expect_sum, na.rm = TRUE),
      surplus_percent = ifelse(count_expect_sum > 0,
                               (count_sum - count_expect_sum) / count_expect_sum * 100, 0)
    ) %>%
    select(admin_area_3, indicator_common_id, period_id,
           count_sum, count_expect_sum, 
           shortfall_absolute, shortfall_percent,
           surplus_absolute, surplus_percent) %>%
    arrange(indicator_common_id, period_id)
  
  write.csv(key_messages_dataset_admin3, "M3_all_indicators_shortfalls_admin_area_3.csv", row.names = FALSE)
  
} else {
  print("No admin_area_3 data available or RUN_DISTRICT_MODEL=FALSE - creating empty key messages file for compatibility...")
  dummy_key_messages_admin3 <- data.frame(
    admin_area_3 = character(0),
    indicator_common_id = character(0),
    period_id = integer(0),
    count_sum = numeric(0),
    count_expect_sum = numeric(0),
    shortfall_absolute = numeric(0),
    shortfall_percent = numeric(0),
    surplus_absolute = numeric(0),
    surplus_percent = numeric(0)
  )
  write.csv(dummy_key_messages_admin3, "M3_all_indicators_shortfalls_admin_area_3.csv", row.names = FALSE)
  
}

# Export data for external key message calculation - ADMIN_AREA_4 LEVEL 
if (RUN_ADMIN_AREA_4_ANALYSIS && exists("summary_disruption_admin4") && nrow(summary_disruption_admin4) > 0) {
  print("Saving admin area 4 level key messages dataset...")
  
  key_messages_dataset_admin4 <- summary_disruption_admin4 %>%
    mutate(
      year = as.integer(period_id %/% 100),
      month = period_id %% 100,
      shortfall_absolute = pmax(0, count_expect_sum - count_sum, na.rm = TRUE),
      shortfall_percent = ifelse(count_expect_sum > 0, 
                                 (count_expect_sum - count_sum) / count_expect_sum * 100, 0),
      surplus_absolute = pmax(0, count_sum - count_expect_sum, na.rm = TRUE),
      surplus_percent = ifelse(count_expect_sum > 0,
                               (count_sum - count_expect_sum) / count_expect_sum * 100, 0)
    ) %>%
    select(admin_area_4, indicator_common_id, period_id,
           count_sum, count_expect_sum, 
           shortfall_absolute, shortfall_percent,
           surplus_absolute, surplus_percent) %>%
    arrange(indicator_common_id, period_id)
  
  write.csv(key_messages_dataset_admin4, "M3_all_indicators_shortfalls_admin_area_4.csv", row.names = FALSE)
  
} else {
  print("No admin_area_4 data available or RUN_ADMIN_AREA_4_ANALYSIS=FALSE - creating empty key messages file for compatibility...")
  dummy_key_messages_admin4 <- data.frame(
    admin_area_4 = character(0),
    indicator_common_id = character(0),
    period_id = integer(0),
    count_sum = numeric(0),
    count_expect_sum = numeric(0),
    shortfall_absolute = numeric(0),
    shortfall_percent = numeric(0),
    surplus_absolute = numeric(0),
    surplus_percent = numeric(0)
  )
  write.csv(dummy_key_messages_admin4, "M3_all_indicators_shortfalls_admin_area_4.csv", row.names = FALSE)
  
}

### MUST CHANGE M3_all_indicators_shortfalls.csv to >>> M3_all_indicators_shortfalls_admin_area_1.csv <<<

print("=== KEY MESSAGES DATASETS EXPORT COMPLETE ===")
print("Generated key messages files:")

print("- M3_all_indicators_shortfalls_admin_area_1.csv (National level)")
print("- M3_all_indicators_shortfalls_admin_area_2.csv (Regional level)")
print("- M3_all_indicators_shortfalls_admin_area_3.csv (State level)")
print("- M3_all_indicators_shortfalls_admin_area_4.csv (District level)")

rm(data_disruption)

# Save Result Objects ----------------------------------------------------------
print("Saving results...")

# Load and save adjusted data (pass-through for viz)
print("Reloading adjusted data and writing service utilization output...")
data_adjusted <- read.csv(
  "M2_adjusted_data.csv",
  colClasses = c(
    facility_id = "character",
    indicator_common_id = "character",
    period_id = "integer",
    admin_area_2 = "character",
    admin_area_3 = "character",
    count_final_none = "numeric",
    count_final_completeness = "numeric"
  )
)
write.csv(data_adjusted, "M3_service_utilization.csv", row.names = FALSE)
gc(); rm(data_adjusted)

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

# Export summary disruptions at national level (admin_area_1)
print("Saving national level disruption analysis...")
summary_disruption_admin1_export <- summary_disruption_admin1 %>%
  select(
    admin_area_1,
    indicator_common_id,
    period_id,
    count_sum,
    count_expect_sum,
    count_expected_if_above_diff_threshold
  )
write.csv(summary_disruption_admin1_export, "M3_disruptions_analysis_admin_area_1.csv", row.names = FALSE)

# Export summary disruptions at admin_area_2 level 
if (exists("summary_disruption_admin2")) {
  print("Saving regional level disruption analysis...")
  summary_disruption_admin2_export <- summary_disruption_admin2 %>%
    select(
      admin_area_2,
      indicator_common_id,
      period_id,
      count_sum,
      count_expect_sum,
      count_expected_if_above_diff_threshold
    )
  write.csv(summary_disruption_admin2_export, "M3_disruptions_analysis_admin_area_2.csv", row.names = FALSE)
} else {
  print("Creating empty admin_area_2 file for compatibility...")
  dummy_disruption_admin2 <- data.frame(
    admin_area_2 = character(0),
    indicator_common_id = character(0),
    period_id = integer(0),
    count_sum = numeric(0),
    count_expect_sum = numeric(0),
    count_expected_if_above_diff_threshold = numeric(0)
  )
  write.csv(dummy_disruption_admin2, "M3_disruptions_analysis_admin_area_2.csv", row.names = FALSE)
}

# Export summary disruptions at admin_area_3 level (conditional)
if (RUN_DISTRICT_MODEL & exists("summary_disruption_admin3")) {
  print("Saving state level disruption analysis...")
  summary_disruption_admin3_export <- summary_disruption_admin3 %>%
    select(
      admin_area_3,
      indicator_common_id,
      period_id,
      count_sum,
      count_expect_sum,
      count_expected_if_above_diff_threshold
    )
  write.csv(summary_disruption_admin3_export, "M3_disruptions_analysis_admin_area_3.csv", row.names = FALSE)
} else {
  print("Creating empty admin_area_3 file for compatibility...")
  dummy_disruption_admin3 <- data.frame(
    admin_area_3 = character(0),
    indicator_common_id = character(0),
    period_id = integer(0),
    count_sum = numeric(0),
    count_expect_sum = numeric(0),
    count_expected_if_above_diff_threshold = numeric(0)
  )
  write.csv(dummy_disruption_admin3, "M3_disruptions_analysis_admin_area_3.csv", row.names = FALSE)
}

# Export summary disruptions at admin_area_4 level (conditional)
if (RUN_ADMIN_AREA_4_ANALYSIS & exists("summary_disruption_admin4")) {
  print("Saving district level disruption analysis...")
  summary_disruption_admin4_export <- summary_disruption_admin4 %>%
    select(
      admin_area_4,
      indicator_common_id,
      period_id,
      count_sum,
      count_expect_sum,
      count_expected_if_above_diff_threshold
    )
  write.csv(summary_disruption_admin4_export, "M3_disruptions_analysis_admin_area_4.csv", row.names = FALSE)
} else {
  print("Creating empty admin_area_4 file for compatibility...")
  dummy_disruption_admin4 <- data.frame(
    admin_area_4 = character(0),
    indicator_common_id = character(0),
    period_id = integer(0),
    count_sum = numeric(0),
    count_expect_sum = numeric(0),
    count_expected_if_above_diff_threshold = numeric(0)
  )
  write.csv(dummy_disruption_admin4, "M3_disruptions_analysis_admin_area_4.csv", row.names = FALSE)
}

print("=== ANALYSIS COMPLETE ===")
