# Special Metrics Synchronization

This document explains the relationship between metric definitions in this repository and the hardcoded const arrays in `wb-fastr`.

## The Const File

**Location:** `wb-fastr/client/src/generate_visualization/get_style_from_po/_0_conditional_consts.ts`

This file contains arrays of metric IDs that require special rendering behavior in the app:

```typescript
export const SPECIAL_COVERAGE_CHART_METRICS = [...]
export const SPECIAL_PERCENT_CHANGE_CHART_METRICS = [...]
export const SPECIAL_DISRUPTIONS_CHART_METRICS = [...]
export const METRICS_WITH_NEGATIVE_PCT_VALUES = [...]
```

## Why These Arrays Exist

The app needs to know which metrics require special chart rendering before it processes vizPreset configs. These arrays allow the rendering pipeline to apply metric-specific logic (axis scaling, color schemes, data label formatting, etc.) based on the metric ID.

## What Each Array Represents

### SPECIAL_COVERAGE_CHART_METRICS

Metrics with `specialCoverageChart: true` in their vizPreset `config.s`.

**Current metrics:** m4-01-01, m6-01-01, m6-02-01, m6-03-01

**Rendering behavior:** Coverage charts show multiple estimation types (survey-based, projected, HMIS-derived) with specific styling for comparing coverage sources over time.

**How to identify:** Search for `specialCoverageChart: true` in metric files.

### SPECIAL_PERCENT_CHANGE_CHART_METRICS

Metrics with `specialBarChart: true` in their vizPreset `config.s`.

**Current metrics:** m3-01-01

**Rendering behavior:** Bar charts showing year-over-year percent change with specific axis and label formatting.

**How to identify:** Search for `specialBarChart: true` in metric files.

**Note:** The const name says "PERCENT_CHANGE" but the module flag is `specialBarChart`. This naming difference is intentional.

### SPECIAL_DISRUPTIONS_CHART_METRICS

Metrics with `specialDisruptionsChart: true` in their vizPreset `config.s`.

**Current metrics:** m3-02-01, m3-03-01, m3-04-01, m3-05-01

**Rendering behavior:** Disruption charts comparing actual vs expected service volumes with specific visual treatment for showing deviations from baseline.

**How to identify:** Search for `specialDisruptionsChart: true` in metric files.

### METRICS_WITH_NEGATIVE_PCT_VALUES

Metrics whose `postAggregationExpression` can produce negative values.

**Current metrics:** m3-02-02, m3-03-02, m3-04-02, m3-05-02

**Rendering behavior:** Axis scaling and color formatting must handle negative percentages (e.g., -15% meaning service volume 15% below expected).

**How to identify:** Check `postAggregationExpression.expression` for subtraction without `ABS()` wrapper:
- **Negative-capable:** `pct_diff = (count_sum - count_expect_sum)/count_expect_sum`
- **Always positive:** `pct_change = ABS(count_final_none-count_final_outliers)/count_final_none`

## Keeping In Sync

When adding or modifying metrics:

1. **Adding a coverage metric:** If it has `specialCoverageChart: true`, add its ID to `SPECIAL_COVERAGE_CHART_METRICS`

2. **Adding a percent change metric:** If it has `specialBarChart: true`, add its ID to `SPECIAL_PERCENT_CHANGE_CHART_METRICS`

3. **Adding a disruptions metric:** If it has `specialDisruptionsChart: true`, add its ID to `SPECIAL_DISRUPTIONS_CHART_METRICS`

4. **Adding a metric with negative values:** If its `postAggregationExpression.expression` contains subtraction without `ABS()`, add its ID to `METRICS_WITH_NEGATIVE_PCT_VALUES`

5. **Removing a metric:** Remove its ID from any const arrays it appears in

## Verification

Run this from the `wb-fastr-modules` root to verify alignment:

```bash
echo "=== specialCoverageChart ===" && grep -rl "specialCoverageChart: true" m*/_metrics/ | xargs basename -a | sed 's/.ts//' | sort
echo "=== specialBarChart ===" && grep -rl "specialBarChart: true" m*/_metrics/ | xargs basename -a | sed 's/.ts//' | sort
echo "=== specialDisruptionsChart ===" && grep -rl "specialDisruptionsChart: true" m*/_metrics/ | xargs basename -a | sed 's/.ts//' | sort
echo "=== negative pct (no ABS) ===" && grep -rl "postAggregationExpression:" m*/_metrics/*.ts | while read f; do grep -A15 "postAggregationExpression:" "$f" | grep "expression:" | grep -v "ABS" | grep -v "null" && basename "$f" .ts; done
```

Compare output against the arrays in `_0_conditional_consts.ts`.

## Consequences of Misalignment

- **Missing from const array:** Metric renders with default styling instead of special treatment. Charts may have wrong axis bounds, missing visual cues, or incorrect data label formatting.

- **Extra in const array:** App applies special rendering to a metric that doesn't need it. May cause visual glitches or incorrect interpretation.

- **Negative values not declared:** Axis may clip negative values or color scale may not handle below-zero percentages correctly.
