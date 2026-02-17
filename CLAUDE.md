# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains R modules for the FASTR (Fast and Sustained Technical Response) project. The modules process HMIS (Health Management Information System) data through a pipeline of data quality assessment, adjustments, service utilization analysis, and coverage estimates.

## Module Files

The main analysis modules are numbered sequentially:
- **m001_module_data_quality_assessment.R** - Evaluates HMIS data reliability through outlier detection, completeness checks, and consistency analysis
- **m002_module_data_quality_adjustments.R** - Applies corrections based on DQA results
- **m003_module_service_utilization.R** - Analyzes service utilization patterns and disruptions
- **m004_module_coverage_estimates_part1.R** - Coverage estimation (part 1)
- **m004_module_coverage_estimates_part2.R** - Coverage estimation (part 2)
- **m007_module_nhss_scorecard.R** - National Health Sector Scorecard (NHSS) quarterly indicators

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
5. **Module 07**: NHSS Scorecard → Quarterly scorecard indicators from M2 adjusted data + asset files

## Coverage Module Design Notes

### Survey Source Handling (COMPLETED 2025 Jan 14)
1. **Source filter expanded**: `process_survey_data()` now accepts ALL survey sources (not just dhs/mics). Priority order: DHS > other sources > MICS
2. **Console messages added**: Reports which sources are available

### UNWPP Denominator Design
The `dwpp_*` denominators (pregnancy, livebirth, dpt, measles) are specifically for UNWPP-based population estimates. By design:
- **ONLY use UNWPP source columns** (`crudebr_unwpp`, `poptot_unwpp`, `totu1pop_unwpp`)
- **No fallback** to DHS/MICS population data - these are conceptually different denominators
- **Robust handling**: Code checks if columns exist before using them - won't crash if missing
- **Informational messages**: Reports what UNWPP data is available and what's missing

### Database Sources Status
Current sources in `survey_data_unified.csv`:
- `DHS Sub-national` (14,470 rows) - works (maps to "dhs")
- `DHS National` (2,190 rows) - works (maps to "dhs")
- `MICS` (2,262 rows) - works
- `UNWPP` (174 rows) - works
- `DHS` (27 rows) - works
- Any other source - **now supported** (priority: DHS > other > MICS)

### Somaliland Data
- **Survey data**: Uses `source="DHS"` in survey_data_unified.csv
- **Population data**: `crudebr`, `poptot`, `totu1pop` in `population_estimates_only.csv` with `source="UNWPP"`

## Module 3: Disruption Detection Method Comparison Framework

### Comparison Framework Files (Added 2025 Jan 30)

A comprehensive methodological comparison framework for testing detection and expected modeling methods:

- **03_method_comparison_engine.R** - Core functions for detection × expected modeling
- **03_run_full_comparison.R** - Run all 12 combinations for key indicators
- **03_shiny_full_comparison.R** - Interactive Shiny app for team walkthrough
- **03_comparison_report.Rmd** - Static HTML report for sharing

### Study Design Matrix

Tests 3 detection methods × 4 expected modeling methods = 12 combinations:

```
                        DETECTION METHOD
                   ┌─────────┬─────────┬─────────┐
                   │   rlm   │   STL   │  bsts   │
                   │(current)│  (IQR)  │ (prob)  │
┌──────────────────┼─────────┼─────────┼─────────┤
│ feols (current)  │   A1*   │   A2    │   A3    │
├──────────────────┼─────────┼─────────┼─────────┤
│ STL decomposition│   B1    │   B2    │   B3    │
├──────────────────┼─────────┼─────────┼─────────┤
│ bsts forecast    │   C1    │   C2    │   C3    │
├──────────────────┼─────────┼─────────┼─────────┤
│ CausalImpact     │   D1    │   D2    │   D3    │
└──────────────────┴─────────┴─────────┴─────────┘

* A1 = Current Module 3 approach (rlm detection + feols expected)
```

### Key Indicators for Comparison

6 indicators covering high-volume and rare events:
- `penta1` (high volume ~62k/month)
- `anc1` (high volume ~68k/month)
- `opd` (very high volume ~2.2M/month)
- `Maternal_deaths` (rare event ~57/month)
- `still_births` (low volume ~660/month)
- `measles1` (high volume ~61k/month)

### Usage

```r
# 1. Run all 12 combinations
source("03_run_full_comparison.R")

# 2. Generate HTML report
rmarkdown::render("03_comparison_report.Rmd")

# 3. Launch interactive Shiny app
source("03_shiny_full_comparison.R")
launch_full_comparison_app()
```

### Academic References

- **rlm**: Huber (1981); Venables & Ripley (2002)
- **STL**: Cleveland et al. (1990) Journal of Official Statistics
- **bsts**: Scott & Varian (2014) Google
- **CausalImpact**: Brodersen et al. (2015) Annals of Applied Statistics
- **feols**: Berge (2018) fixest package

## Module 7: National Health Sector Scorecard (NHSS)

### Overview

The NHSS (formerly "RMNCAH+N SWAp-Aligned Scorecard") produces quarterly performance indicators at state and national level. Renamed per TWG recommendation to reflect its broader scope beyond reproductive/maternal/child health. Used during quarterly performance dialogues.

**File**: `m007_module_nhss_scorecard.R`
**Parameters**: `SCORECARD_YEAR`, `SCORECARD_QUARTER`, `COUNTRY_ISO3`, `SELECTED_COUNT_VARIABLE`
**Outputs**: `M7_scorecard_combined.csv`, `M7_scorecard_summary_Q{Q}_{YEAR}.csv`

### Data Requirements Stock Take

The scorecard needs data from **two sources**: the standard HMIS pipeline (M2 adjusted data) and four separate asset CSV files. Everything must match on `admin_area_2` (state name) and cover the target quarter.

#### Source 1: `M2_adjusted_data_admin_area.csv` (from Module 02 pipeline)

Required columns: `period_id` (YYYYMM), `admin_area_3` (state), `indicator_common_id`, `{SELECTED_COUNT_VARIABLE}`

The following 27 `indicator_common_id` values are needed. If an indicator is missing the scorecard column goes to NA — it does not crash.

| `indicator_common_id` | Used in scorecard indicator(s) | Category |
|----------------------|-------------------------------|----------|
| `anc1` | B (ANC 1-visit coverage), C (ANC4/ANC1 ratio), N (IPTp3 coverage) | ANC |
| `anc4` | C (ANC4/ANC1 ratio) | ANC |
| `delivery` | D (Skilled birth attendance), E (Uterotonics), F (Fistula rate) | Delivery |
| `sba` | D (Skilled birth attendance) | Delivery |
| `uterotonics` | E (Uterotonics coverage) | Delivery |
| `obstetric_fistula` | F (Fistula per 10000 deliveries) | Delivery |
| `birth_asphyxia` | G (Newborn resuscitation) | Newborn |
| `neonatal_resuscitation` | G (Newborn resuscitation) | Newborn |
| `live_births` | H (Postnatal visits), I (LBW KMC) | Newborn |
| `pnc` | H (Postnatal visits within 3 days) | Postnatal |
| `kmc` | I (LBW KMC coverage) | Newborn |
| `child_registration` | J (Birth registration) | Child |
| `modern_fp` | K (Modern contraceptive use) | FP |
| `pneumonia` | L (Pneumonia treatment) | Child health |
| `pneumonia_treatment` | L (Pneumonia treatment) | Child health |
| `diarrhoea` | M (Diarrhea ORS+zinc treatment) | Child health |
| `ors_zinc` | M (Diarrhea ORS+zinc treatment) | Child health |
| `iptp3` | N (IPTp3 coverage) | Malaria |
| `llin` | O (Under-5 LLIN coverage) | Malaria |
| `fully_immunized` | O (LLIN, denominator), Q (Fully immunized coverage) | Immunization |
| `bcg` | P (BCG coverage) | Immunization |
| `ebf` | S (Exclusive breastfeeding) | Nutrition |
| `nutrition_screening` | T (Growth monitoring) | Nutrition |
| `gbv_cases` | U (GBV care coverage) | GBV |
| `gbv_care` | U (GBV care coverage) | GBV |
| `opd` | V (Facility utilisation per 100 person-years) | Utilisation |
| `measles1` | W (Under-1 fully immunised, MCV1 proxy) | Immunization |

#### Source 2: Asset CSV files (from `scorecard_module_prepare_assets.R`)

Each asset file must have columns: `admin_area_2`, `period_id` (YYYYMM), `count`

| Asset file | Scorecard indicator | Notes |
|-----------|-------------------|-------|
| `vaccine_stockout_pct.csv` | R: Vaccine stockout % | Averaged across quarter months per state |
| `nhmis_2019_reporting_rate.csv` | X: NHMIS reporting rate | Averaged across quarter months per state |
| `nhmis_data_timeliness.csv` | Y: NHMIS data timeliness | Averaged across quarter months per state |
| `total_population.csv` | Denominator for B, D-F, I-K, P-Q, S-T, V-W | Fallback: current year → previous year Q4 average |

**Note**: These are NOT from the standard DHIS2 facility-month pipeline. Ideally they should migrate there. Currently prepared manually via `scorecard_module_prepare_assets.R`.

#### Population Assumptions (hardcoded quarterly rates)

These convert total population into age/sex-specific denominators:

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `PREGNANT_WOMEN_PCT` | 0.05 × 0.25 = 0.0125 | Quarterly pregnant women as fraction of pop |
| `BIRTHS_PCT` | 0.04 × 0.25 = 0.01 | Quarterly births as fraction of pop |
| `WOMEN_15_49_PCT` | 0.22 × 0.25 = 0.055 | Quarterly women 15-49 as fraction of pop |
| `CHILDREN_U5_PCT` | 0.12 × 0.25 = 0.03 | Quarterly children under 5 as fraction of pop |
| `INFANTS_0_6M_PCT` | 0.02 | Infants 0-6 months as fraction of pop |

### Current Indicator List (23 indicators)

| # | Variable | Formula | Source |
|---|----------|---------|--------|
| B | `anc_coverage_1_visit` | anc1 / (pop × PREGNANT_WOMEN_PCT) × 100 | HMIS |
| C | `anc4_anc1_ratio` | anc4 / anc1 × 100 | HMIS |
| D | `skilled_birth_attendance` | sba / delivery × 100 | HMIS |
| E | `uterotonics_coverage` | uterotonics / delivery × 100 | HMIS |
| F | `fistula_per_1000_deliveries` | obstetric_fistula / (delivery / 10000) | HMIS |
| G | `newborn_resuscitation` | neonatal_resuscitation / birth_asphyxia × 100 | HMIS |
| H | `postnatal_visits_3d` | pnc / live_births × 100 | HMIS |
| I | `lbw_kmc_coverage` | kmc / (live_births × 0.15) × 100 | HMIS |
| J | `birth_registration` | child_registration / (pop × BIRTHS_PCT) × 100 | HMIS |
| K | `modern_contraceptive_use` | modern_fp / (pop × WOMEN_15_49_PCT) × 100 | HMIS |
| L | `pneumonia_antibiotic_treatment` | pneumonia_treatment / pneumonia × 100 | HMIS |
| M | `diarrhea_ors_zinc_treatment` | ors_zinc / diarrhoea × 100 | HMIS |
| N | `iptp3_coverage` | iptp3 / anc1 × 100 | HMIS |
| O | `under5_llin_coverage` | llin / fully_immunized × 100 | HMIS |
| P | `bcg_coverage` | bcg / (pop × BIRTHS_PCT) × 100 | HMIS |
| Q | `fully_immunized_coverage` | fully_immunized / (pop × BIRTHS_PCT) × 100 | HMIS |
| R | `vaccine_stockout_percentage` | from asset file | Asset |
| S | `exclusive_breastfeeding_rate` | ebf / (pop × INFANTS_0_6M_PCT) × 100 | HMIS |
| T | `growth_monitoring_coverage` | nutrition_screening / (pop × 0.20 × 0.25) × 100 | HMIS |
| U | `gbv_care_coverage` | gbv_care / gbv_cases × 100 | HMIS |
| V | `facility_utilisation_opd` | opd / (pop × 0.25) × 100 | HMIS |
| W | `under1_fully_immunised_mcv1` | measles1 / (pop × BIRTHS_PCT) × 100 | HMIS |
| X | `nhmis_reporting_rate_final` | from asset file | Asset |
| Y | `nhmis_data_timeliness_final` | from asset file | Asset |

### TWG Change Log

- **2026 Feb 17**: TWG-approved changes — renamed from "RMNCAH+N SWAp-Aligned Scorecard" to "National Health Sector Scorecard (NHSS)". Removed: Malaria ACT treatment rate, Penta3 coverage. Added: Facility utilisation (OPD per 100 person-years), Under-1 fully immunised (MCV1 proxy). Fixed: `delivery` column name (was `deliveries`). Confirmed institutional denominators for SBA (D) and ANC4/ANC1 (C). Parameterized country ISO3 and made summary table quarter-dynamic.
