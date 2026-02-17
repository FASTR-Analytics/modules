# National Health Sector Scorecard (NHSS) — Changes Following Q4 NM&E TWG

**Date**: February 2026
**Based on**: Report of Q4 NM&E TWG meeting recommendations
**Updated by**: Claire Boulange

---

## 1. What changed?

The TWG reviewed the scorecard (previously called the "RMNCAH+N SWAp-Aligned Scorecard") and approved a set of changes. The code has now been updated to reflect those decisions.

### Rename

The scorecard is now called the **National Health Sector Scorecard (NHSS)** to better reflect its scope — it covers more than just reproductive, maternal, newborn, and child health.

### Indicators removed (2)

| Removed indicator | Why |
|---|---|
| **Malaria ACT treatment rate** (ACT treatments / confirmed malaria cases) | TWG decision — not needed on this scorecard |
| **Penta3 coverage** (Penta3 doses / estimated births) | TWG decision — BCG and fully immunised already capture immunisation performance |

### Indicators added (2)

| New indicator | What it measures | How it's calculated |
|---|---|---|
| **Facility utilisation (OPD)** | How much the population uses outpatient services | Total OPD visits in the quarter, divided by the population, expressed per 100 people per year |
| **Under-1 fully immunised (MCV1 proxy)** | Proxy for full immunisation of infants | Measles 1st dose (MCV1) in the quarter, divided by the estimated number of births |

### Bug fix

The code was looking for a column called `deliveries` but the actual data column is called `delivery` (singular). This meant that three indicators — Skilled Birth Attendance, Uterotonics Coverage, and Fistula Rate — were silently returning blank values. Now fixed.

---

## 2. How does the scorecard turn monthly facility data into quarterly state-level numbers?

The HMIS data comes in at the **facility level, every month**. The scorecard needs **one number per state, per quarter**. Here's how:

1. **Filter to the target quarter** — e.g. for Q1 2025, keep only January, February, and March 2025
2. **Sum across all facilities in each state** — all facility counts for a given indicator are added up within each state
3. **Calculate the indicator** — divide one count by another (e.g. ANC4 visits / ANC1 visits) to get a rate or percentage

For indicators that use population as the denominator (like "ANC coverage"), the population is adjusted to represent **one quarter's worth** rather than a full year. For example, if roughly 4% of the population gives birth per year, then in one quarter that's about 1% (= 4% divided by 4 quarters). This is why you see population multiplied by small fractions like 0.01 — it's the annual rate divided by 4.

The national row is calculated by summing all state-level raw counts first, then computing the indicator — so it's a proper weighted national figure, not an average of state percentages.

---

## 3. What are the "asset files" and why are they separate?

Most scorecard indicators are calculated directly from the HMIS data that flows through the standard pipeline (facility-level, monthly, extracted from DHIS2). But **four pieces of information** do not come from that pipeline:

| Asset file | What it contains | Why it's separate |
|---|---|---|
| **Vaccine stockout data** | Percentage of facilities with vaccine stockouts, by state | This is a facility-level stockout indicator, not a service count — it's structured differently from routine HMIS |
| **NHMIS reporting rate** | Percentage of expected reports that were actually submitted, by state | This is a system performance metric about DHIS2 itself, not a health service count |
| **NHMIS data timeliness** | Percentage of reports submitted on time, by state | Same as above — measures reporting behaviour, not service delivery |
| **Total population** | State population estimates | Population projections come from national statistics / UNWPP, not from facility reporting |

These four files are prepared separately (using `scorecard_module_prepare_assets.R`) and loaded alongside the HMIS data when the scorecard runs. Each file has the same simple structure: state name, month (YYYYMM format), and a value.

**In the future**, it would be cleaner to integrate these into the standard DHIS2 data pipeline so everything comes from one source. For now they are maintained as separate CSV files.

---

## 4. Full indicator list after changes (23 indicators)

| # | Indicator | How it's calculated | Data source |
|---|-----------|-------------------|-------------|
| B | ANC coverage (1st visit) | ANC1 visits / estimated pregnant women | HMIS + population |
| C | ANC4 / ANC1 ratio | ANC4 visits / ANC1 visits | HMIS only |
| D | Skilled birth attendance | Skilled births / total deliveries | HMIS only |
| E | Uterotonics coverage | Uterotonics given / total deliveries | HMIS only |
| F | Fistula rate | Obstetric fistula cases per 10,000 deliveries | HMIS only |
| G | Newborn resuscitation | Resuscitations / birth asphyxia cases | HMIS only |
| H | Postnatal visits (within 3 days) | PNC visits / live births | HMIS only |
| I | LBW KMC coverage | KMC cases / (live births x 15% LBW estimate) | HMIS only |
| J | Birth registration | Registrations / estimated births | HMIS + population |
| K | Modern contraceptive use | Modern FP users / estimated women 15-49 | HMIS + population |
| L | Pneumonia treatment | Treatments / pneumonia cases | HMIS only |
| M | Diarrhea ORS+zinc treatment | ORS+zinc / diarrhea cases | HMIS only |
| N | IPTp3 coverage | IPTp3 doses / ANC1 visits | HMIS only |
| O | Under-5 LLIN coverage | LLINs distributed / fully immunised children | HMIS only |
| P | BCG coverage | BCG doses / estimated births | HMIS + population |
| Q | Fully immunised coverage | Fully immunised / estimated births | HMIS + population |
| R | Vaccine stockout % | From asset file | Asset file |
| S | Exclusive breastfeeding | EBF cases / estimated infants 0-6 months | HMIS + population |
| T | Growth monitoring | Nutrition screenings / estimated children under 5 | HMIS + population |
| U | GBV care coverage | GBV care received / GBV cases reported | HMIS only |
| V | **Facility utilisation (OPD)** *(NEW)* | OPD visits per 100 person-years | HMIS + population |
| W | **Under-1 fully immunised — MCV1 proxy** *(NEW)* | MCV1 doses / estimated births | HMIS + population |
| X | NHMIS reporting rate | From asset file | Asset file |
| Y | NHMIS data timeliness | From asset file | Asset file |

**Removed** (previously on the scorecard):
- ~~Malaria ACT treatment rate~~
- ~~Penta3 coverage~~

---

## 5. Parameters to set before running

These are all set at the top of the script. Change them before each run.

### Run parameters

| Parameter | Current value | What to change when |
|---|---|---|
| `SCORECARD_YEAR` | 2025 | Change to the year you're reporting on |
| `SCORECARD_QUARTER` | 1 | Change to the quarter (1-4) |
| `COUNTRY_ISO3` | "NGA" | Change if running for a different country |
| `SELECTED_COUNT_VARIABLE` | "count_final_none" | Which adjusted count to use from the pipeline |

### Population assumptions (demographic fractions)

Several indicators need to estimate a target population from the total state population — for example, "how many pregnant women are there in this state this quarter?" Since we don't have that number directly, we estimate it as a fraction of the total population.

Each fraction is the **annual rate divided by 4** to get a quarterly figure (except infants 0-6 months which is already a point-in-time proportion).

| Parameter | Annual assumption | Quarterly value | What it represents | Which indicators use it |
|---|---|---|---|---|
| `PREGNANT_WOMEN_PCT` | 5% of pop pregnant per year | 0.0125 | Estimated pregnant women this quarter | B: ANC coverage |
| `BIRTHS_PCT` | 4% of pop gives birth per year | 0.01 | Estimated births this quarter | J: Birth registration, P: BCG, Q: Fully immunised, W: MCV1 proxy |
| `WOMEN_15_49_PCT` | 22% of pop are women 15-49 | 0.055 | Estimated women of reproductive age this quarter | K: Modern contraceptive use |
| `INFANTS_0_6M_PCT` | 2% of pop at any time | 0.02 | Estimated infants 0-6 months (not annualised — point-in-time) | S: Exclusive breastfeeding |
| `CHILDREN_U5_PCT` | 12% of pop are under 5 | 0.03 | Estimated children under 5 this quarter | *Legacy parameter — kept for future use, not currently used by any indicator* |

**Note on `CHILDREN_U5_PCT`**: This parameter exists from legacy code. The LLIN indicator (O) logically could use it, but currently follows the original Excel formula which uses `fully_immunized` as its denominator instead. The parameter is kept in case a future revision switches LLIN or another indicator to use a proper population-based under-5 denominator.
