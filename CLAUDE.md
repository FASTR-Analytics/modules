# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains R modules for the FASTR (Fast and Sustained Technical Response) project. The modules process HMIS (Health Management Information System) data through a pipeline of data quality assessment, adjustments, service utilization analysis, and coverage estimates.

## Module Files

The main analysis modules are numbered sequentially:
- **01_module_data_quality_assessment.R** - Evaluates HMIS data reliability through outlier detection, completeness checks, and consistency analysis
- **02_module_data_quality_adjustments.R** - Applies corrections based on DQA results
- **03_module_service_utilization.R** - Analyzes service utilization patterns and disruptions
- **04_module_coverage_estimates_part1.R** - Coverage estimation (part 1)
- **04_module_coverage_estimates_part2.R** - Coverage estimation (part 2)

## Code Editing Guidelines

### IMPORTANT: Always Update the "Last edit" Date

**When editing any module file (*.R), you MUST update the date in the header comment:**

```r
# Last edit: YYYY MMM DD
```

**Format**: Use 4-digit year, 3-letter month abbreviation, 2-digit day (e.g., "2025 Nov 04")

**Example**:
```r
# CB - R code FASTR PROJECT
# Last edit: 2025 Nov 04
# Module: DATA QUALITY ASSESSMENT
```

This helps track when changes were made and ensures proper version documentation.

### Other Guidelines

- Module files use parameters defined at the top of each script (COUNTRY_ISO3, thresholds, etc.)
- Output files follow naming convention: `M{module_number}_output_{description}.csv`
- Each module is designed to run independently but follows the sequential pipeline
- Geographic columns follow pattern: `admin_area_1`, `admin_area_2`, etc.

## Data Pipeline Flow

1. **Module 01**: Data Quality Assessment → Outputs outlier flags, completeness scores, consistency checks
2. **Module 02**: Data Adjustments → Applies corrections based on DQA results
3. **Module 03**: Service Utilization → Analyzes patterns and disruptions
4. **Module 04**: Coverage Estimates → Calculates coverage using adjusted data and denominators

## TODO: Coverage Module Improvements

### Code Changes Needed
1. **Refactor source filter in `process_survey_data()`**: Currently only accepts `source %in% c("dhs", "mics")`. Should accept configurable list or broader set of sources (slhds, unicef, etc.)

2. **Refactor crudebr handling**: Code hardcodes `crudebr_unwpp` at lines 638-641. Should handle crudebr from any source, not just UNWPP.

3. **Refactor poptot handling**: Similarly hardcodes `poptot_unwpp`. Should be source-agnostic.

### Database Sources to Standardize
Current sources in `survey_data_unified.csv`:
- `DHS Sub-national` (14,470 rows) - works (maps to "dhs")
- `DHS National` (2,190 rows) - works (maps to "dhs")
- `MICS` (2,262 rows) - works
- `UNWPP` (174 rows) - works
- `DHS` (27 rows) - works
- `UNICEF Other` (18 rows) - **needs code support**

### Temporary Workarounds Applied
- **Somaliland**: Survey data uses `source="DHS"` instead of "SLHDS"
- **Somaliland crudebr**: Uses `source="UNWPP"` to create `crudebr_unwpp` column
