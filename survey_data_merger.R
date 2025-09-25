# ========================================
# SURVEY DATA MERGER TOOL
# Author: Data Integration Specialist
# Purpose: Intelligently merge new survey data with existing datasets
# ========================================

library(dplyr)
library(data.table)
library(readr)
library(lubridate)
library(stringr)

# ========================================
# CONFIGURATION
# ========================================

# File paths
SURVEY_DATA_FILE <- "/Users/claireboulange/Desktop/modules/survey_data_unified.csv"
POPULATION_DATA_FILE <- "/Users/claireboulange/Desktop/modules/population_estimates_only.csv"
BACKUP_DIR <- "/Users/claireboulange/Desktop/modules/backups"
CLEANED_DATA_DIR <- "/Users/claireboulange/Desktop/modules/cleaned_data_downloads"

# Create directories if they don't exist
if(!dir.exists(BACKUP_DIR)) {
  dir.create(BACKUP_DIR, recursive = TRUE)
  message("📁 Created backup directory: ", BACKUP_DIR)
}

if(!dir.exists(CLEANED_DATA_DIR)) {
  dir.create(CLEANED_DATA_DIR, recursive = TRUE)
  message("📁 Created cleaned data directory: ", CLEANED_DATA_DIR)
}

# ========================================
# UTILITY FUNCTIONS
# ========================================

create_backup <- function(file_path) {
  if(file.exists(file_path)) {
    backup_name <- paste0(basename(file_path), "_backup_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    backup_path <- file.path(BACKUP_DIR, backup_name)
    file.copy(file_path, backup_path)
    message("💾 Backup created: ", backup_path)
    return(backup_path)
  }
  return(NULL)
}

validate_data_structure <- function(df, expected_cols) {
  missing_cols <- setdiff(expected_cols, names(df))
  if(length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  extra_cols <- setdiff(names(df), expected_cols)
  if(length(extra_cols) > 0) {
    message("⚠️  Extra columns found (will be ignored): ", paste(extra_cols, collapse = ", "))
  }
  
  return(TRUE)
}

standardize_data_types <- function(df) {
  df %>%
    mutate(
      admin_area_1 = as.character(admin_area_1),
      admin_area_2 = as.character(admin_area_2),
      year = as.integer(year),
      indicator_id = as.character(indicator_id),
      indicator_common_id = as.character(indicator_common_id),
      indicator_type = as.character(indicator_type),
      survey_value = as.numeric(survey_value),
      source = as.character(source),
      source_detail = as.character(source_detail),
      survey_type = as.character(survey_type)
    )
}

# ========================================
# DATA LOADING FUNCTIONS
# ========================================

load_existing_data <- function(file_path, data_type = "survey") {
  if(!file.exists(file_path)) {
    message("📄 File doesn't exist yet: ", file_path)
    return(data.frame())
  }
  
  message("📥 Loading existing ", data_type, " data from: ", basename(file_path))
  
  tryCatch({
    df <- fread(file_path)
    message("✅ Loaded ", nrow(df), " existing records")
    
    # Validate structure
    expected_cols <- c("admin_area_1", "admin_area_2", "year", "indicator_id", 
                       "indicator_common_id", "indicator_type", "survey_value", 
                       "source", "source_detail", "survey_type")
    validate_data_structure(df, expected_cols)
    
    # Standardize data types
    df <- standardize_data_types(df)
    
    # Report summary
    message("📊 Data summary:")
    print(table(df$source))
    message("Countries: ", length(unique(df$admin_area_1)))
    message("Years: ", paste(range(df$year, na.rm = TRUE), collapse = " - "))
    
    return(df)
    
  }, error = function(e) {
    stop("❌ Error loading ", file_path, ": ", e$message)
  })
}

find_new_data_files <- function(download_dir = CLEANED_DATA_DIR) {
  # Look for recently downloaded cleaned data files
  pattern <- "cleaned_survey_data_.*\\.csv$"
  files <- list.files(download_dir, pattern = pattern, full.names = TRUE)
  
  if(length(files) == 0) {
    message("❌ No cleaned data files found in ", download_dir)
    message("💡 Looking for files matching pattern: ", pattern)
    message("📁 Expected location: ", CLEANED_DATA_DIR)
    message("💾 To use this folder, download cleaned data from Shiny app to: ", basename(CLEANED_DATA_DIR))
    return(character(0))
  }
  
  # Sort by modification time (most recent first)
  file_info <- file.info(files)
  files <- files[order(file_info$mtime, decreasing = TRUE)]
  
  message("📁 Found ", length(files), " potential data files in ", basename(download_dir), ":")
  for(i in seq_along(files)) {
    mod_time <- format(file_info[files[i], "mtime"], "%Y-%m-%d %H:%M")
    message("  ", i, ". ", basename(files[i]), " (", mod_time, ")")
  }
  
  return(files)
}

load_new_data <- function(file_path = NULL) {
  if(is.null(file_path)) {
    # Auto-find the most recent file
    files <- find_new_data_files()
    if(length(files) == 0) {
      stop("❌ No new data files found")
    }
    file_path <- files[1]  # Most recent
    message("🎯 Auto-selected most recent file: ", basename(file_path))
  }
  
  if(!file.exists(file_path)) {
    stop("❌ File not found: ", file_path)
  }
  
  message("📥 Loading new data from: ", basename(file_path))
  
  tryCatch({
    df <- fread(file_path)
    message("✅ Loaded ", nrow(df), " new records")
    
    # Validate structure
    expected_cols <- c("admin_area_1", "admin_area_2", "year", "indicator_id", 
                       "indicator_common_id", "indicator_type", "survey_value", 
                       "source", "source_detail", "survey_type")
    validate_data_structure(df, expected_cols)
    
    # Standardize data types
    df <- standardize_data_types(df)
    
    # Report summary
    message("📊 New data summary:")
    print(table(df$source))
    print(table(df$indicator_type))
    message("Countries: ", paste(unique(df$admin_area_1), collapse = ", "))
    message("Years: ", paste(unique(df$year), collapse = ", "))
    
    return(df)
    
  }, error = function(e) {
    stop("❌ Error loading new data: ", e$message)
  })
}

# ========================================
# DATA MERGING FUNCTIONS
# ========================================

create_record_key <- function(df) {
  df %>%
    mutate(
      record_key = paste(admin_area_1, admin_area_2, year, indicator_common_id, source, sep = "_")
    )
}

identify_duplicates <- function(existing_df, new_df) {
  existing_keys <- create_record_key(existing_df)$record_key
  new_with_keys <- create_record_key(new_df)
  
  duplicates <- new_with_keys$record_key %in% existing_keys
  
  message("🔍 Duplicate analysis:")
  message("  - New records: ", nrow(new_df))
  message("  - Duplicates found: ", sum(duplicates))
  message("  - New unique records: ", sum(!duplicates))
  
  if(sum(duplicates) > 0) {
    dup_summary <- new_with_keys[duplicates, ] %>%
      group_by(source, admin_area_1) %>%
      summarise(count = n(), .groups = "drop")
    
    message("📋 Duplicates by source and country:")
    print(dup_summary)
  }
  
  return(list(
    duplicates = duplicates,
    unique_new = new_with_keys[!duplicates, ],
    duplicate_records = new_with_keys[duplicates, ]
  ))
}

merge_datasets <- function(existing_df, new_df, allow_updates = FALSE) {
  message("\n🔄 MERGING DATASETS")
  message("=" , paste(rep("=", 50), collapse = ""))
  
  # Identify duplicates
  dup_analysis <- identify_duplicates(existing_df, new_df)
  
  if(nrow(dup_analysis$unique_new) == 0) {
    message("ℹ️  No new unique records to add")
    return(existing_df)
  }
  
  # Clean up the unique new records (remove the record_key column)
  unique_new_clean <- dup_analysis$unique_new %>%
    select(-record_key)
  
  # Merge the data
  merged_df <- bind_rows(existing_df, unique_new_clean)
  
  message("✅ Merge completed!")
  message("  - Original records: ", nrow(existing_df))
  message("  - New records added: ", nrow(unique_new_clean))
  message("  - Total records: ", nrow(merged_df))
  
  # Handle duplicates if requested
  if(allow_updates && nrow(dup_analysis$duplicate_records) > 0) {
    message("⚠️  Note: ", nrow(dup_analysis$duplicate_records), " duplicate records were ignored")
    message("💡 To update existing records, implement update logic here")
  }
  
  return(merged_df)
}

separate_by_indicator_type <- function(df) {
  population_data <- df %>%
    filter(indicator_type == "population_estimate")
  
  survey_data <- df %>%
    filter(indicator_type != "population_estimate")
  
  message("📊 Data separation:")
  message("  - Population estimates: ", nrow(population_data), " records")
  message("  - Survey data: ", nrow(survey_data), " records")
  
  return(list(
    population = population_data,
    survey = survey_data
  ))
}

# ========================================
# MAIN MERGER FUNCTION
# ========================================

merge_survey_data <- function(new_data_file = NULL, create_backups = TRUE, dry_run = FALSE) {
  
  message("\n🚀 SURVEY DATA MERGER TOOL")
  message("=" , paste(rep("=", 60), collapse = ""))
  message("Started at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  
  # Step 1: Create backups
  if(create_backups) {
    message("\n1️⃣ CREATING BACKUPS")
    create_backup(SURVEY_DATA_FILE)
    create_backup(POPULATION_DATA_FILE)
  }
  
  # Step 2: Load existing data
  message("\n2️⃣ LOADING EXISTING DATA")
  existing_survey <- load_existing_data(SURVEY_DATA_FILE, "survey")
  existing_population <- load_existing_data(POPULATION_DATA_FILE, "population")
  
  # Combine existing data for analysis
  existing_combined <- bind_rows(existing_survey, existing_population)
  message("📊 Total existing records: ", nrow(existing_combined))
  
  # Step 3: Load new data
  message("\n3️⃣ LOADING NEW DATA")
  new_data <- load_new_data(new_data_file)
  
  # Step 4: Separate new data by type
  message("\n4️⃣ SEPARATING DATA BY TYPE")
  separated_new <- separate_by_indicator_type(new_data)
  
  # Step 5: Merge survey data
  message("\n5️⃣ MERGING SURVEY DATA")
  merged_survey <- merge_datasets(existing_survey, separated_new$survey)
  
  # Step 6: Merge population data
  message("\n6️⃣ MERGING POPULATION DATA")
  merged_population <- merge_datasets(existing_population, separated_new$population)
  
  # Step 7: Save results (unless dry run)
  if(!dry_run) {
    message("\n7️⃣ SAVING UPDATED FILES")
    
    if(nrow(merged_survey) > nrow(existing_survey)) {
      fwrite(merged_survey, SURVEY_DATA_FILE)
      message("💾 Updated survey data file: ", nrow(merged_survey), " total records")
    } else {
      message("ℹ️  No changes to survey data")
    }
    
    if(nrow(merged_population) > nrow(existing_population)) {
      fwrite(merged_population, POPULATION_DATA_FILE)
      message("💾 Updated population data file: ", nrow(merged_population), " total records")
    } else {
      message("ℹ️  No changes to population data")
    }
  } else {
    message("\n7️⃣ DRY RUN - NO FILES SAVED")
    message("Would save:")
    message("  - Survey data: ", nrow(merged_survey), " records")
    message("  - Population data: ", nrow(merged_population), " records")
  }
  
  # Step 8: Final summary
  message("\n✅ MERGE COMPLETE!")
  message("=" , paste(rep("=", 60), collapse = ""))
  message("Final totals:")
  message("  📋 Survey data: ", nrow(merged_survey), " records")
  message("  📈 Population data: ", nrow(merged_population), " records")
  message("  🌍 Countries: ", length(unique(c(merged_survey$admin_area_1, merged_population$admin_area_1))))
  message("  📅 Years: ", paste(range(c(merged_survey$year, merged_population$year), na.rm = TRUE), collapse = " - "))
  
  # Return results for inspection
  return(list(
    survey_data = merged_survey,
    population_data = merged_population,
    summary = list(
      survey_records = nrow(merged_survey),
      population_records = nrow(merged_population),
      new_survey_added = nrow(merged_survey) - nrow(existing_survey),
      new_population_added = nrow(merged_population) - nrow(existing_population)
    )
  ))
}

# ========================================
# INTERACTIVE FUNCTIONS
# ========================================

interactive_merge <- function() {
  cat("\n🎯 INTERACTIVE SURVEY DATA MERGER\n")
  cat("=====================================\n")
  
  # Find available files
  files <- find_new_data_files()
  
  if(length(files) == 0) {
    cat("❌ No new data files found in Downloads folder\n")
    cat("💡 Please ensure you have downloaded cleaned data from the fetcher\n")
    return(invisible())
  }
  
  # Let user choose file
  if(length(files) > 1) {
    cat("\nSelect a file to merge:\n")
    for(i in seq_along(files)) {
      mod_time <- format(file.info(files[i])$mtime, "%Y-%m-%d %H:%M")
      cat("  ", i, ". ", basename(files[i]), " (", mod_time, ")\n")
    }
    
    choice <- as.integer(readline("Enter file number: "))
    if(is.na(choice) || choice < 1 || choice > length(files)) {
      cat("❌ Invalid choice\n")
      return(invisible())
    }
    chosen_file <- files[choice]
  } else {
    chosen_file <- files[1]
    cat("📁 Using: ", basename(chosen_file), "\n")
  }
  
  # Ask about dry run
  dry_run <- tolower(readline("Do you want to do a dry run first? (y/n): ")) == "y"
  
  # Run the merge
  result <- merge_survey_data(new_data_file = chosen_file, dry_run = dry_run)
  
  if(dry_run) {
    proceed <- tolower(readline("\nProceed with actual merge? (y/n): ")) == "y"
    if(proceed) {
      result <- merge_survey_data(new_data_file = chosen_file, dry_run = FALSE)
    }
  }
  
  return(result)
}

# ========================================
# HELPER FUNCTIONS FOR FILE MANAGEMENT
# ========================================

open_cleaned_data_folder <- function() {
  system(paste("open", shQuote(CLEANED_DATA_DIR)))
  message("📁 Opened cleaned data folder: ", CLEANED_DATA_DIR)
}

check_downloads_folder <- function() {
  downloads_dir <- "~/Downloads"
  pattern <- "cleaned_survey_data_.*\\.csv$"
  files <- list.files(downloads_dir, pattern = pattern, full.names = TRUE)
  
  if(length(files) == 0) {
    message("📁 No cleaned data files found in Downloads folder")
    return(invisible())
  }
  
  message("📁 Found ", length(files), " cleaned data file(s) in Downloads:")
  for(i in seq_along(files)) {
    mod_time <- format(file.info(files[i])$mtime, "%Y-%m-%d %H:%M")
    message("  ", i, ". ", basename(files[i]), " (", mod_time, ")")
  }
  
  move_all <- tolower(readline("Move all files to cleaned data folder? (y/n): ")) == "y"
  
  if(move_all) {
    for(file in files) {
      move_file_to_cleaned_folder(file)
    }
    message("✅ All files moved to cleaned data folder")
  }
  
  return(files)
}

move_file_to_cleaned_folder <- function(file_path) {
  if(!file.exists(file_path)) {
    stop("❌ File not found: ", file_path)
  }
  
  file_name <- basename(file_path)
  dest_path <- file.path(CLEANED_DATA_DIR, file_name)
  
  success <- file.copy(file_path, dest_path, overwrite = TRUE)
  if(success) {
    message("✅ Moved file to: ", dest_path)
    file.remove(file_path)
    message("🗑️  Removed original file from: ", file_path)
  } else {
    stop("❌ Failed to move file")
  }
}

# Quick status check
check_data_status <- function() {
  message("📊 CURRENT DATA STATUS")
  message("=" , paste(rep("=", 40), collapse = ""))
  
  if(file.exists(SURVEY_DATA_FILE)) {
    survey_df <- fread(SURVEY_DATA_FILE)
    message("📋 Survey data: ", nrow(survey_df), " records")
    message("   Sources: ", paste(names(table(survey_df$source)), collapse = ", "))
    message("   Countries: ", length(unique(survey_df$admin_area_1)))
    message("   Years: ", paste(range(survey_df$year, na.rm = TRUE), collapse = " - "))
  } else {
    message("📋 Survey data: File not found")
  }
  
  if(file.exists(POPULATION_DATA_FILE)) {
    pop_df <- fread(POPULATION_DATA_FILE)
    message("📈 Population data: ", nrow(pop_df), " records")
    message("   Sources: ", paste(names(table(pop_df$source)), collapse = ", "))
    message("   Countries: ", length(unique(pop_df$admin_area_1)))
    message("   Years: ", paste(range(pop_df$year, na.rm = TRUE), collapse = " - "))
  } else {
    message("📈 Population data: File not found")
  }
  
  # Check for new files
  new_files <- find_new_data_files()
  if(length(new_files) > 0) {
    message("🆕 Available new files: ", length(new_files))
  } else {
    message("🆕 No new files found")
  }
}

# ========================================
# USAGE EXAMPLES
# ========================================

cat("
🚀 SURVEY DATA MERGER TOOL - USAGE GUIDE
=========================================

## Quick Start:
interactive_merge()          # Interactive mode with prompts

## Manual merge:
result <- merge_survey_data()  # Auto-find most recent file
result <- merge_survey_data('path/to/file.csv')  # Specify file

## Utilities:
check_data_status()          # Check current data status
find_new_data_files()       # Find available files
check_downloads_folder()    # Check Downloads and offer to move files
open_cleaned_data_folder()  # Open the cleaned data folder in Finder
move_file_to_cleaned_folder('~/Downloads/file.csv')  # Move file to correct location

## Options:
merge_survey_data(dry_run = TRUE)     # Preview changes only
merge_survey_data(create_backups = FALSE)  # Skip backups

## Examples:
# 1. Quick interactive merge
interactive_merge()

# 2. Check what you have
check_data_status()

# 3. Merge latest file
result <- merge_survey_data()

📁 Files managed:
- Survey data: ", basename(SURVEY_DATA_FILE), "
- Population data: ", basename(POPULATION_DATA_FILE), "
- Backups: ", basename(BACKUP_DIR), "
- Cleaned downloads: ", basename(CLEANED_DATA_DIR), "

💡 IMPORTANT: Save your cleaned data downloads to:
   ", CLEANED_DATA_DIR, "

")