# Chain-Based Denominator Selection — Team Notes

## What we do now

We pick the "best" denominator **separately for each population group**. The pregnancy denominator might come from ANC1, the DPT denominator from deliveries, the livebirth denominator from penta1. Each group picks whatever source fits the survey data best, independently.

This means our population estimates can be **internally inconsistent**. We might estimate 100,000 pregnancies (from ANC1) but 120,000 live births (from deliveries). That's demographically impossible — you can't have more live births than pregnancies.

## What the chain-based approach would do

Pick **one HMIS indicator** as the anchor for the whole country. Use it to derive ALL population estimates through the demographic cascade:

```
anchor → pregnancies → deliveries → births → live births → DPT-eligible → measles-eligible
```

For example, if "delivery" is the selected chain:
- Pregnancies = delivery count / delivery survey rate, adjusted backward through cascade
- Live births = delivery count / delivery survey rate
- DPT-eligible = live births × (1 - neonatal mortality)
- Measles-eligible = DPT-eligible × (1 - post-neonatal mortality)

Everything flows from one number. The population estimates are guaranteed to be consistent with each other.

### Which anchors are eligible?

| Anchor | What it measures | Can derive all groups? |
|--------|-----------------|----------------------|
| ANC1 | First antenatal visits → pregnancies | Yes (forward through full cascade) |
| Delivery | Institutional deliveries → live births | Yes (backward to pregnancies, forward to DPT/measles) |
| SBA | Skilled birth attendance → live births | Yes (same as delivery, uses delivery survey rate) |
| Live births | Reported live births | Yes (same as delivery) |

Not eligible:
- **Penta1**: Only gives DPT + measles groups (can't go backward to pregnancies/livebirths)
- **BCG**: National-level only (no subnational data)
- **UNWPP**: National-level only; used as a scoring benchmark instead

### How we'd pick the best chain

Score each chain on two things:

**Score A — Are the numbers in the right ballpark?**
Compare chain denominators against UNWPP population estimates. A chain whose denominators are 50% off from UNWPP scores worse than one that's 10% off.

**Score B — Does the coverage match DHS/MICS?**
For each indicator, compare coverage (using the chain's denominator) against the survey value. **Exclude the anchor indicator itself** (see "circularity" below).

Normalize both scores to 0–1, combine 50/50. Lowest total score wins.

---

## How chain selection would work in practice

### Step 1: Calculate all denominators (same as now)

The code already calculates denominators from every anchor — `danc1_pregnancy`, `ddelivery_livebirth`, `dpenta1_dpt`, etc. This doesn't change. We just stop picking winners per group and instead evaluate whole chains.

### Step 2: Score each chain

For each eligible chain (anc1, delivery, sba, livebirths), calculate two scores:

**Score A — Magnitude check against external benchmarks**

For each population group (pregnancy, livebirth, DPT, measles), compare the chain's denominator value against UNWPP:

```
relative_error = |chain_denominator - UNWPP_estimate| / UNWPP_estimate
```

Average across all available groups. A chain whose denominators are consistently close to UNWPP scores low (good). One that's way off scores high (bad).

This catches obvious problems — e.g., a chain that implies 2 million pregnancies in a country of 3 million people.

**Score B — Coverage fit against DHS/MICS**

For each indicator that has survey data, calculate the squared error between coverage (using the chain's denominator) and the survey value:

```
squared_error = (HMIS_count / chain_denominator - survey_value)²
```

**Critical rule:** Exclude the anchor indicator from this score. If the chain is "delivery", don't score delivery coverage — it's circular and always perfect (see Issue #1 below).

Average across all non-anchor indicators.

**Combined score**

Normalize each score to 0–1 across chains (so the best chain gets 0, worst gets 1), then average:

```
total = 0.5 × normalized_magnitude + 0.5 × normalized_survey_fit
```

Lowest total wins. Tiebreaker: prefer anc1 > delivery > sba > livebirths (by typical HMIS reporting quality).

### Step 3: Apply the selected chain everywhere

Once a chain is selected at national level:
- National uses that chain's denominators for all indicators
- Admin2 uses the same chain
- Admin3 uses the same chain
- No fallback logic needed — the eligible chains all work at subnational level

### What the code would look like

The scoring logic replaces the current `ranked_denom_per_group` block. Roughly:

```r
# For each chain, gather its denominators across all population groups
# Score A: compare against UNWPP
# Score B: compare coverage against survey (excluding anchor indicator)
# Normalize, combine, pick winner
# Then: filter all coverage data to only use the winning chain's denominators
```

m006 replaces its 7 `DENOM_*` parameters with a single `SELECTED_CHAIN` parameter (or reads the chain from the M5 output).

---

## How to present chain selection visually

### Chain comparison dot plot

A single plot that shows, for each population group in the cascade, what each chain estimates vs what external benchmarks say.

```
                       Population estimate (thousands)
                  0    200    400    600    800    1000
                  |------|------|------|------|------|
Pregnancies       ●-----------●--------○---------◇
Live births            ●------●--------○------◇
DPT-eligible           ●-----●--------○-----◇
Measles-eligible       ●-----●-------○-----◇

● = chain estimates (one colour per chain)
○ = UNWPP benchmark
◇ = DHIS2/census benchmark (if available)
```

**What to look for:**
- Chains that cluster near the benchmarks = good
- Chains that drift further from benchmarks as you go down the cascade = error amplification
- Chains where all dots are on the same side of the benchmark (all too high, or all too low) = systematic bias

The selected chain gets highlighted (bold/larger dot). The legend shows each chain's combined score.

### Coverage comparison panel

For each indicator, show the coverage estimate from each chain alongside the DHS/MICS survey value:

```
         Coverage (%)
         0%    25%    50%    75%   100%
ANC1     |------●=====●======▲-----|     ← anchor for anc1 chain (circular)
Penta1   |----------●--●---●-▲----|
Penta3   |--------●--●---●---▲----|
BCG      |----------●---●--●-▲----|
Delivery |------●====●=======▲----|     ← anchor for delivery chain (circular)

● = coverage from each chain (colour-coded)
▲ = DHS/MICS survey value
```

This shows: does the coverage from my selected chain actually match DHS? Circular values (anchor indicator for its own chain) will sit exactly on the survey triangle — you can see that visually.

### Year-by-year denominator trend

For the selected chain, plot the denominator values over time alongside UNWPP and any census data:

```
Estimated live births over time

800k ─                              ╱ UNWPP
     │                    ○ ○ ○ ○ ○
700k ─           ● ● ● ●
     │     ● ● ●                    chain estimate
600k ─ ● ●
     │         ◇                     census (if available)
500k ─
     └──────────────────────────────
     2010  2012  2014  2016  2018  2020  2022
```

This shows whether the chain's population estimates are tracking the expected growth trajectory or diverging.

---

## External data we could use to validate chain selection

### What we have now

| Benchmark | Coverage | Level | Used for |
|-----------|----------|-------|----------|
| **UNWPP** (UN World Population Prospects) | All countries | National only | Score A — magnitude check |
| **DHS/MICS surveys** | Most countries | National + some subnational | Score B — coverage fit |

### What we could add

#### 1. DHIS2 population targets

Most countries load population targets into their DHIS2 for their own coverage calculations. These are census-based projections from the national statistics office — not UN modeled numbers.

**Why useful:**
- Available at **district level** (admin3/admin4) — gives us subnational validation that UNWPP can't
- Represents what the MoH actually plans against — operationally relevant
- Often disaggregated by age group (under-1, under-5, women 15-49)

**What we need:** Extract from each country's DHIS2: total population, expected pregnancies, expected live births, population under 1, population under 5. At lowest available geographic level, for all available years.

**Current status:**

| Country | DHIS2 pop data | Status |
|---------|---------------|--------|
| Nigeria | Total population by facility, monthly 2020-2025 | Have total pop; need age breakdown |
| Zambia | National total population, annual 2010-2023 | Have total pop; need subnational + age groups |
| Guinea | Not extracted | Need to pull |
| Sierra Leone | Not extracted | Need to pull |

#### 2. National census data

When a recent census exists, it's the most reliable population count available. Census data can validate whether chain denominators are reasonable in absolute terms.

**Examples:**
- Zambia 2022 census: 19.6M total (vs UNWPP projection of 20.2M — 2.7% gap)
- Nigeria 2006 census + projections

Census data is typically total population only, but combined with known birth rates it can cross-check livebirth denominators.

#### 3. WHO/UNICEF Estimates of National Immunization Coverage (WUENIC)

WUENIC provides annual national coverage estimates for key vaccines (BCG, DPT3, MCV1, etc.) derived from both survey and administrative data. Could serve as a third benchmark for Score B — does our chain's coverage match the international consensus estimate?

**Advantage:** Available for all countries, annually updated
**Limitation:** National level only, and itself derived from similar inputs (so not fully independent)

#### 4. Civil registration and vital statistics (CRVS)

Countries with functioning birth/death registration systems have direct counts of births and deaths. Where available, these are a strong check on livebirth denominators.

**Current reality:** CRVS completeness varies widely. Many FASTR countries have <50% registration completeness, limiting usefulness. But where completeness is high (e.g., South Africa), it's very valuable.

### How additional benchmarks would integrate into scoring

Score A currently only uses UNWPP. With additional benchmarks:

```
For each chain:
  unwpp_error       = mean relative error vs UNWPP (national)
  dhis2_nat_error   = mean relative error vs DHIS2 targets (national)
  dhis2_sub_error   = mean relative error vs DHIS2 targets (across all admin2/3 units)
  census_error      = relative error vs census-derived estimate (if available)

  Score A = average of whichever components are available
```

The subnational DHIS2 component naturally carries more weight because it averages across many geographic units — this is good because subnational fit is a harder test.

Score B (survey coverage fit) stays unchanged — still DHS/MICS based.

---

## Issues and trade-offs

### 1. The anchor indicator loses its independent coverage estimate

This is the biggest thing to accept.

If "ANC1" is the selected chain, then ANC1 coverage is calculated as:

```
ANC1 coverage = ANC1 count / ANC1 denominator
              = ANC1 count / (ANC1 count / ANC1 survey rate)
              = ANC1 survey rate
```

It just returns the survey value. The HMIS data contributes nothing to the ANC1 coverage estimate — it cancels out.

**What this means in practice:**
- The anchor indicator's coverage line in the output will just be the DHS/MICS value carried forward. No HMIS-based trend for that indicator.
- For all OTHER indicators, coverage is a genuine combination of HMIS counts and survey-calibrated denominators.

**Is this acceptable?** The anchor is the indicator we trust most for population sizing — that's why it was selected. We're effectively saying "we believe the ANC1 survey rate, and we use it to calibrate everything else." The trade-off is that we can't also independently estimate ANC1 coverage from the same data.

### 2. One bad HMIS indicator can pull everything off

If the selected chain's anchor has a data quality problem in one year (say deliveries drop due to a reporting gap), ALL denominators for that year shift. Under the current per-group approach, only the delivery-derived denominators would be affected.

**Mitigation:** This is partly why we run M1 (outlier detection) and M2 (adjustments) before M4/M5. The adjusted counts should already have reporting gaps corrected. But it's still a concentration of risk.

### 3. The cascade amplifies errors

Small errors in the anchor compound through each demographic adjustment step:

```
pregnancies → (× pregnancy loss rate) → deliveries → (× twin rate) → births → (× stillbirth rate) → live births → ...
```

Each step multiplies by a demographic rate that is itself an estimate. A chain that enters at "pregnancies" (ANC1) has more steps to reach "DPT-eligible" than one that enters at "live births" (delivery). More steps = more accumulated error.

**Counter-argument:** The current per-group approach has this same cascade math — it just picks different entry points per group, hiding the accumulated error behind inconsistent population estimates.

### 4. Fewer denominators available at subnational level

BCG and UNWPP denominators are national-only. Under the current approach, if the best denominator is national-only, we fall back to the second-best for subnational. Under chain-based, the chain is fixed — if the selected chain doesn't work well at subnational level, we're stuck with it.

**Mitigation:** Chain selection could incorporate subnational fit into the scoring (e.g., if DHIS2 population targets are available at district level, score how well each chain matches them across all districts).

### 5. What do we show for the anchor indicator?

Options:
- **Option A**: Show the survey carry value as the coverage estimate. Be transparent that it's not HMIS-derived. Label it differently in the output (e.g., "survey_reference" instead of "coverage_cov").
- **Option B**: Don't output coverage for the anchor indicator at all. Only show it for indicators where the denominator is independent.
- **Option C**: Show it but flag it. Add a column like `is_anchor = TRUE/FALSE` so downstream users know.

### 6. SBA as a chain depends on delivery survey data

SBA doesn't have its own survey coverage rate in most DHS datasets. When SBA is the selected chain, it uses the delivery survey rate as a proxy. This means the SBA chain is really just the delivery chain with a different numerator — not truly independent.

---

## What stays the same

- The demographic cascade formulas (stillbirth rate, NMR, twin rate, etc.) don't change
- UNWPP denominators still exist as benchmarks, just not as selectable denominators
- Survey data processing (DHS/MICS priority, carry-forward) stays the same
- Subnational analysis still uses the national-level denominator selection
- All three output levels (national, admin2, admin3) still produced

## What changes

| Current | Chain-based |
|---------|------------|
| Best denominator picked per population group independently | One chain selected for all groups |
| m006 has 7 `DENOM_*` parameters | Single `SELECTED_CHAIN` parameter |
| Coverage estimates exist for all indicators | Anchor indicator coverage = survey value (circular) |
| Population estimates may be inconsistent across groups | All population estimates flow from one source |
| National-only fallback logic for subnational | Same chain everywhere, no fallback needed |

## Questions for the team

1. Is demographic coherence important enough to accept losing the independent coverage estimate for one indicator?
2. Should we show the anchor indicator's (circular) coverage in the output, or suppress/flag it?
3. Are we comfortable concentrating risk on a single HMIS indicator for all denominators?
4. Do we want to keep the per-group approach as a fallback or comparison mode?
