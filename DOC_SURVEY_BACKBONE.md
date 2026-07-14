# DOC: Survey & Population Backbone — Semantics, Traps, and Maintenance

Reference for `survey_data_unified.csv` and `population_estimates_only.csv`.

Read this before adding, renaming, or deleting anything in either file. It exists because
a full audit on **2026 Jul 14** found several indicators whose *name did not match what the
column actually held*. None of them crashed anything. That is exactly why they survived.

---

## 1. The two files

| File | Holds | Rows (2026 Jul 14) |
|---|---|---|
| `survey_data_unified.csv` | survey indicators | 26,312 |
| `population_estimates_only.csv` | population denominators | 2,689 |

### Both files use the SAME writer — always

```r
write.csv(df, path, row.names = FALSE, na = "")   # fully quoted
```

**Never use `readr::write_csv` on either database.** Two reasons:

1. Minimal quoting rewrites every line of a 26k-row file, so the real change disappears inside a
   spurious full-file diff.
2. `readr::write_csv` writes the literal string **`NA`** for missing values unless you pass `na=""`.
   Swapping the writer doesn't just churn quotes — it can turn missing data into the text "NA".

> **History (2026 Jul 14).** The two files used to use *different* writers — survey via `write.csv`,
> population via `readr::write_csv` — and the rule was documented but not enforced. It was already
> being violated: `moz_fix_routing.R` wrote **both** files with `readr::write_csv`, and so did
> `append_to_databases()` in `R/integration_functions.R` — *the main path the Shiny app uses*. The
> app itself would have reformatted the entire survey file. Both databases were unified on
> `write.csv` and all 10 integration scripts updated. **The rule is now "one writer, both files."
> Do not reintroduce a per-file style.**

**Schema (13 cols, both files):** `admin_area_1, admin_area_2, year, indicator_id,
indicator_common_id, indicator_type, survey_value, source, source_detail, survey_type,
country_name, iso2_code, iso3_code`

### Routing

These eight `POP_INDICATORS` go in the population file, everything else in the survey file:

```
poptot, popu5, totu1pop, totu5pop, livebirth, womenrepage, popgrowth, crudebr
```

This applies **even to DHS-sourced rows**. It is enforced by `append_to_databases()` in the
fetcher; a hand-rolled append must replicate it.

### The fetcher APPENDS, and can REPLACE — but it never regenerates

`append_to_databases()` does two things:

1. **Appends** genuinely new rows — `bind_rows(survey_db, new_survey_data)`.
2. **Replaces** rows the user marks for update — `anti_join` the old row out, then bind the new one.

It never rewrites rows you did not target. So the consequence people get wrong still holds:
**a bad row already in the CSV can only be fixed in the CSV.** Fixing the fetcher stops *new* bad
rows; it does nothing to the ~26k already there. Most fixes need doing in **both** places. Fixing
only one is the commonest failure mode here.

### The replace key MUST include `source`

The DB keeps **one value per `(admin_area_1, admin_area_2, year, indicator_common_id, source)`**.
Different sources are *designed* to coexist: DHS + WUENIC `penta3`, DHS + UNWPP `crudebr`,
MICS + UNWPP `poptot`. **353 keys / 706 rows** in the current DBs hold two or more sources.

> **Bug fixed 2026 Jul 14.** `detect_duplicates()` built its `composite_key` from only four columns
> — **no `source`** — and `append_to_databases()` `anti_join`ed on the same four. Three consequences:
> a new WUENIC row was mis-flagged as a duplicate of the existing DHS row; the `left_join` fanned one
> new record out across every source's row; and choosing "update" **deleted every source at that key**
> and wrote back only the new one. Verified against the live DB: updating Afghanistan's 2015 DHS
> `bcg` value silently destroyed the WUENIC `bcg` row for the same year. Both keys now include
> `source`.

---

## 2. Indicators whose name does not mean what you think

Every one of these was a real bug found in the audit. Read this section before trusting a column.

### `womenrepage` is a population estimate, NOT a sample count

Until 2026 Jul 14 this single name held **two quantities four orders of magnitude apart**:

| | rows | range | what it really was |
|---|---|---|---|
| UNWPP (`52`) | 331 | 35,796 – 57,393,553 | female population aged 15–49 |
| DHS (`FE_FRTY_W_NPG`) | 1,045 | 27 – 42,221 | **number of women interviewed** |

A denominator built from that column without filtering on `source` would have been catastrophically
wrong. Now split:

- **`womenrepage`** — UNWPP female population 15–49. `population_estimate`. **Population file.** A real denominator.
- **`women_interviewed`** — DHS unweighted sample count. `survey_count`. **Survey file.** Never a denominator.

The DHS favourite `FE_FRTY_W_NPG` was commented "Women of reproductive age", which is what made
it look interchangeable. It is not.

### `fp` was never family planning coverage — the name is now RETIRED

`fp` held `FP_SRCM_W_TOT`, the **total** of the "source of modern method" distribution, which is
100% by construction. All 138 rows sat between 0.993 and 1.000. It was a denominator check
someone mistook for a prevalence — the fetcher comment literally read `# mCPR`.

**Deleted.** Do not reintroduce the name. Correct codes:

| Concept | DHS code | common_id |
|---|---|---|
| Modern contraceptive use, **currently married women** | `FP_CUSM_W_MOD` | (choose one, see below) |
| Modern contraceptive use, **all women** | `FP_CUSA_W_MOD` | `contraceptive_modern` |
| Any method | `FP_CUSA_W_ANY` | `contraceptive_any` |
| Unmet need | `FP_NADA_W_UNT` | `unmet_need` |
| Demand satisfied | `MNCH_DEMAND_FP` | `fp_demand_satisfied` |

⚠️ **Open decision before any FP re-fetch.** `FP_CUSA_*` is *all women*; `FP_CUSM_*` is *currently
married women*. They are not comparable and DHS reports both. EDHS 2024-25 leads with the married-women
figure (mCPR 34.5%). **Pick one denominator and hold it across all countries.** Mixing them silently
produces a league table of nonsense. There is currently **no DHS survey-based FP coverage in the
backbone at all** — `mcpr` (331 rows) is UNWPP *modelled* data, not survey data.

### `anc_none` is the INVERSE of ANC coverage (was `anc1_old`)

`rh_ancn_w_n01` is the "**zero** ANC visits" category. Median 0.036, against `anc1` median 0.894.
It was called `anc1_old`, which invites reading it as an older ANC1 definition. **Renamed to
`anc_none`.** It is `1 - coverage`. Do not plot it on the same axis as `anc1`.

### `stillbirth` holds COUNTS, not rates

`indicator_type = number`. Harari and Gambella record `0`; Nigeria records 567. A stillbirth
*rate* is ~15–25 per 1,000 and can never be zero. These are **numbers of events**.

This is why the Ethiopia EDHS 2024-25 stillbirth figure (9 per 1,000, a genuine **rate**) was
**not** imported. Writing it into this column would have turned Amhara's 62 stillbirths into "11",
and the result would have looked entirely plausible. If you want a stillbirth *rate*, create a new
`indicator_common_id`. Do not reuse this one.

### `u5mr` is a rate even where it was typed `percent`

24 Madagascar rows were typed `percent` while holding values from 31 to **107.6**. A percent of
107.6 is impossible; they were rates. Fixed to `rate`. Any code dividing percents by 100 would have
shrunk Madagascar's under-5 mortality by two orders of magnitude.

---

## 3. Canonical vocabulary — one concept, one name

Three concepts were each stored under two names, split along the national/subnational boundary.
Verified as true duplicates (where they overlapped, values were **identical**: 20/20, 6/6, 3/3),
then merged.

| Concept | Retired name | **Canonical name** | DHS code |
|---|---|---|---|
| Fertility | ~~`tfr`~~ | **`total_fertility_rate`** | `FE_FRTR_W_TFR` |
| Birth rate | ~~`crude_birth_rate`~~ | **`crudebr`** | `FE_FRTR_W_CBR` |
| Stillbirths | ~~`still`~~ | **`stillbirth`** | `CM_PNMR_C_NSB` |
| FP source total | ~~`fp`~~ | **`fp_source_total`** (not fetched) | `FP_SRCM_W_TOT` |
| No ANC | ~~`anc1_old`~~ | **`anc_none`** | `rh_ancn_w_n01` |
| Women interviewed | ~~`womenrepage`(DHS)~~ | **`women_interviewed`** | `FE_FRTY_W_NPG` |

**Never reintroduce a retired name.** The fetcher's mapping table has been updated so it cannot.

---

## 4. Known-broken, NOT fixed here (needs an owner)

### Rotavirus is dead on both sides

- m004 asks for `rota1` / `rota2` (script.R lines 174, 189, 204, 298).
- The backbone stores the DHS names `rotavirus1` / `rotavirus2` (596 rows, 15 countries).
- m004's `recode()` block at line 288 translates `polio1→opv1`, `vitamina→vitaminA`… **but has no
  rota entry.** So the filter matches zero rows.
- Separately, **no `hmis_*.csv` in the repo carries a rota dose count at all.**

So even adding the two recode lines yields nothing — there is no HMIS numerator to pair against.
Rota coverage has returned empty for every country since the indicator was added, without warning.
Fix needs **two lines in m004 (module code) AND a DHIS2 pull for rota doses.** Not fixable from the
data side. **Do not "fix" this by renaming `rotavirus1`→`rota1` in the CSV** — that diverges from the
fetcher and from the polio/vitamin A pattern, and still produces nothing.

### Admin joins that silently produce zero coverage

m004 joins survey→HMIS on `admin_area_2` with **no fuzzy matching and no stripping**. One character
off = zero coverage for that region, silently.

| Country | Issue | Verdict |
|---|---|---|
| **MRT** | survey has one `Nouakchott`; HMIS has `Nouakchott Sud` + `Nord` + `Ouest` | **Do not auto-split.** Splitting one survey value across three districts invents data the survey never collected. |
| **TCD** | `hmis_TCD.csv` has `WADIFIRA`; the backbone **and** the survey both say `WADI FIRA` | The **HMIS extract** is wrong, not the survey. Fix belongs to whoever owns the Chad DHIS2 pull. |
| **AFG** | HMIS keeps provinces at `admin_area_4`; survey puts them at `admin_area_2` | No subnational join possible. |
| **MWI** | HMIS `admin_area_2` = 5 zones; districts live at `admin_area_3`; survey puts districts at `admin_area_2` | No subnational join possible. |
| **ETH** | legacy `SNNPR` rows | Region no longer exists (split 2020+). Left as historical record; joins to nothing. |
| **CIV, TCD, LBR** | `Rural`, `Other urban`, `Zone 1…8`, `Foret Rurale`, `Rest of Country` | DHS **analytical strata, not places.** Never meant to join. Not a bug. |

`"N'Djaména 1997"` in the Chad rows is a useful reminder: **DHS year-stamps region names when
boundaries move.** Two entries with the same name and different years are not necessarily the same
geographic entity.

### `06_survey_data_fetcher/assets/ethiopia_backbone.csv` is STALE

13 entries, and it never got the SNNPR split. Missing `Central Ethiopian region` and
`South West Ethiopia Region`. Has trailing whitespace (`'Amhara Region '`), and calls Benishangul a
`Regional Health Bureau`. **Any future ETH fetch harmonising to this backbone will produce wrong names.**

Until it is refreshed: **match the real HMIS (`01_hmis_data/hmis_ethiopia_q2.csv`), not the backbone.**

---

## 5. Ethiopia specifics

### The `year` column is in the ETHIOPIAN CALENDAR

Ethiopia is the **only** country in the file with this offset. Every ET survey is stored at
**Gregorian year − 8**, and the HMIS periods run `200912`–`201705` to match.

| Survey | `source_detail` | stored `year` |
|---|---|---|
| EDHS 2000 | `ET2000DHS` | 1992 |
| EDHS 2011 | `ET2011DHS` | 2003 |
| EDHS 2019 | `ET2019DHS` | 2011 |
| **EDHS 2024-25** | **`ET2025DHS`** | **2017** |

Writing `2025` would break the HMIS join silently. (Other countries show −1/−2 offsets; that is just
DHS fieldwork start year, not a calendar difference.)

### Region structure

EDHS 2024-25 uses the **new 14-region structure** (SNNPR split into Central Ethiopia, South Ethiopia,
South West Ethiopia, Sidama). The Ethiopian HMIS already matches it. Names in the survey file are the
**HMIS spellings**, which differ from the EDHS report's:

| EDHS report | survey file / HMIS |
|---|---|
| Tigrai | `Tigray Region` |
| Benishangul-Gumuz | `Benishangul Gumuz Region` |
| Central Ethiopia | `Central Ethiopian region` *(lowercase "region" — sic)* |
| Addis Ababa | `Addis Ababa City Administration` |

Tigray **was** surveyed and has full estimates. (The dashboard spells it "Tigrai", which is easy to
miss when grepping.)

---

## 6. Pre-write validation checklist

Run **all** of these and show the numbers before writing. Every one of them caught something real.

- [ ] schema is the 13 cols; `country_name` / `iso2_code` / `admin_area_1` consistent per `iso3_code`
- [ ] 0 NAs in `admin_area_1, year, indicator_common_id, survey_value, source`
- [ ] 0 duplicate keys on `(admin_area_1, admin_area_2, year, indicator_common_id, source)`
      — **note `admin_area_1` is in the key**; a country-name case difference will *hide* duplicates
      from you (this is exactly how 236 duplicate Madagascar rows went unnoticed)
- [ ] `percent` values within [0, 1]; `rate` / `number` values sane
- [ ] **one `indicator_type` per `indicator_common_id`** — a mixed type means two concepts share a name
- [ ] subnational `admin_area_2` ⊆ the HMIS/backbone (0 unmatched)
- [ ] `POP_INDICATORS` only in the population file, and vice versa
- [ ] correct per-file quoting (see §1) — diff should be the size of your change, not 25k lines
- [ ] **reconcile**: replay your transforms from the pristine backup and assert the result is
      byte-identical to what you wrote. Anything unexplained is a bug you have not noticed yet.

Always take a timestamped backup first: `*.csv.bak-YYYYMMDD-HHMMSS`.

---

## 7. Traps that have already bitten

1. **`parse()` passing means nothing for a positional table.** `get_indicator_mapping()` in
   `cleaning_functions.R` is a `data.frame()` of two *positionally paired* vectors (`original_id`
   ↔ `common_id`). Replacing one element with three silently mis-maps every indicator after it.
   R syntax-checks fine. **Always `source()` the file and assert `nrow()` plus spot-check the
   entries after your edit point.**

2. **Moving rows between the two files carries un-harmonised metadata with them.** Rows moved out of
   the population file reintroduced a `Central African Republic` country name and a bad
   `IRS Nzérékoré` apostrophe that had *already been fixed* in the survey file. **Re-run the
   harmonisation after any cross-file move**, then re-run the audit.

3. **"It's the same concept" is not the same as "it's a duplicate."** Before merging two names,
   check whether they ever cover the same `(country, area, year, source)` and whether they *agree*
   there. If they disagree, they are different indicators wearing similar names — do not merge.

4. **A plausible comment on the wrong code is the hardest bug to see.** `FP_SRCM_W_TOT # mCPR` and
   `FE_FRTY_W_NPG # Women of reproductive age` both survived for years because the comment described
   what someone *wanted* the code to be. Verify the code against the DHS definition, not the comment.

---

## 8. Change log

### 2026 Jul 14 — backbone audit + Ethiopia EDHS 2024-25

**Ethiopia EDHS 2024-25 loaded** — 330 rows, 22 indicators × 15 areas (14 regions + NATIONAL),
`year=2017` (Ethiopian calendar), `source_detail=ET2025DHS`. Values verified against the official
**Final Report FR399 (July 2026)**: 11 of 12 national figures match exactly (facility delivery differs
by 0.1). `stillbirth` deliberately **excluded** — EDHS reports a rate, the column holds counts (§2).

`imr` (40) and `nmr` (26) **were** imported despite not matching a published national table exactly;
FR399 publishes 39 and 25 for the *5 years* before the survey and no 10-year national total. The
values sit correctly on Ethiopia's own trend and the worst case is a mislabelled reference period,
not a wrong magnitude. Flagged for Nicole Danfakah / Dan Kabtyimer (ARO Data & Digital) to confirm.

**Data fixes** (survey file 25,920 → 26,312 rows; population file 3,166 → 2,689):

| Change | Rows |
|---|---|
| deleted `fp` (`FP_SRCM_W_TOT`, always ~100%) | −138 |
| deleted rows with empty `survey_value` (all CAF "RS 6" 1994) | −10 |
| renamed `anc1_old` → `anc_none` | 1,182 |
| split `womenrepage` → `women_interviewed` (DHS sample counts), moved to survey file | 1,045 |
| merged `tfr`→`total_fertility_rate`, `crude_birth_rate`→`crudebr`, `still`→`stillbirth` | 540 |
| deduped after merge + Madagascar name harmonisation | −267 |
| fixed `u5mr` `indicator_type` percent→rate (MDG) | 24 |
| fixed `indicator_type` `other`→real type (`fully_immunized`, `pnc1`) | 169 |
| fixed GIN `IRS Nzérékoré` → `IRS N'zérékoré` (backbone is source of truth) | 71 |
| fixed `iso2_code` holding an ISO-3 (MRT, CIV) | 216 |
| harmonised `country_name` / `admin_area_1` (CAF, CIV, MDG) | 1,834 |

**Fetcher fixes** (`06_survey_data_fetcher`, local commits — not pushed):

- `data_functions.R` — removed `FP_SRCM_W_TOT` from DHS favourites (it was labelled `# mCPR`);
  corrected the `FE_FRTY_W_NPG` comment.
- `cleaning_functions.R` — `fe_frty_w_npg` → `women_interviewed`; `fp_srcm_w_tot` → `fp_source_total`;
  `MNCH_DEMAND_FP` → `fp_demand_satisfied`; TFR unified to `total_fertility_rate`; Guinea province
  mapping target corrected to `IRS N'zérékoré` (**the fetcher was writing the bad name on purpose**);
  removed the label fallback that bucketed any FP/contraceptive label into `fp`.
- `indicator_mappings.R` — UNWPP `"19"` and the TFR pattern now emit `total_fertility_rate`.

**Not fixed, needs an owner:** rota (m004 + DHIS2), MRT Nouakchott, AFG/MWI admin levels,
stale `ethiopia_backbone.csv`, `hmis_TCD.csv` `WADIFIRA` typo, FP denominator decision + re-fetch.
