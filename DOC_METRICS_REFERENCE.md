# Metrics Reference

How metrics work in wb-fastr, and how to author them correctly.

---

## Mental Model

A **metric** defines how to aggregate columns from a results object into meaningful
values, plus how to visualize those values. When a user views a metric:

1. **SQL query** aggregates `valueProps` using `valueFunc`, grouped by disaggregation dimensions
2. **Post-aggregation** (optional) applies a formula to the aggregated values
3. **Visualization** renders the data according to a vizPreset's config

```
Results CSV        Metric Definition           SQL Query                    Visualization
┌────────────┐     ┌──────────────────┐        ┌─────────────────────┐      ┌─────────────┐
│ facility_id│     │ valueProps:      │        │ SELECT              │      │ Table/Chart │
│ period_id  │     │   ["cases"]      │        │   admin_area_2,     │      │ showing     │
│ admin_area │ ──► │ valueFunc: "SUM" │ ──►    │   SUM(cases)        │ ──►  │ aggregated  │
│ cases      │     │ disaggregateBy:  │        │ FROM results        │      │ cases by    │
│ tests      │     │   [admin_area_2] │        │ GROUP BY admin_area │      │ region      │
└────────────┘     └──────────────────┘        └─────────────────────┘      └─────────────┘
```

---

## Core Fields

### valueProps

Array of column names from the results object to aggregate.

```typescript
valueProps: ["cases", "tests"]
```

These become the "values" in the visualization. If you have multiple valueProps,
each becomes a separate series/column in charts and tables.

### valueFunc

How to aggregate the valueProps. One of:

| Function   | SQL              | Use when                                    |
|------------|------------------|---------------------------------------------|
| `"SUM"`    | `SUM(col)`       | Counting totals (cases, visits, doses)      |
| `"AVG"`    | `AVG(col)`       | Proportions, rates, averages                |
| `"COUNT"`  | `COUNT(col)`     | Counting rows                               |
| `"MIN"`    | `MIN(col)`       | Finding minimums                            |
| `"MAX"`    | `MAX(col)`       | Finding maximums                            |
| `"identity"` | `col` (no agg) | Pre-aggregated values, single-row results  |

**Important:** `valueFunc` applies to ALL valueProps equally. If you need different
aggregation functions for different columns, use `postAggregationExpression`.

### resultsObjectId

Which results object (CSV file) this metric reads from.

```typescript
resultsObjectId: "M1_output_outliers.csv"
```

Must match an `id` in `_results_objects.ts`.

### formatAs

How to format the aggregated values:

- `"percent"` — multiply by 100, show % symbol (use for proportions 0-1)
- `"number"` — show as-is with thousands separators

---

## Disaggregation

### requiredDisaggregationOptions

Dimensions that MUST be included in any query for this metric. The metric
doesn't make sense without them.

```typescript
requiredDisaggregationOptions: ["facility_type"]
```

**Example:** A metric showing "cases by facility type" must always group by
`facility_type`. Users can add more disaggregations, but this one is mandatory.

### disaggregateBy (in vizPresets)

Which dimensions to show in THIS specific visualization, and how to display them.

```typescript
disaggregateBy: [
  { disOpt: "indicator_common_id", disDisplayOpt: "col" },
  { disOpt: "admin_area_2", disDisplayOpt: "row" },
]
```

### disDisplayOpt values

Where the dimension appears visually:

| Value        | Table                      | Timeseries Chart         | Bar Chart            |
|--------------|----------------------------|--------------------------|----------------------|
| `"row"`      | Row headers                | Y-axis grouping          | Y-axis grouping      |
| `"rowGroup"` | Grouped row headers        | —                        | —                    |
| `"col"`      | Column headers             | X-axis panes             | X-axis panes         |
| `"colGroup"` | Grouped column headers     | —                        | —                    |
| `"series"`   | —                          | Multiple lines           | Grouped bars         |
| `"cell"`     | Cell grouping              | Grid of small charts     | Grid of small charts |
| `"indicator"`| —                          | —                        | Bar categories       |
| `"replicant"`| Separate tables            | Small multiples          | Small multiples      |
| `"mapArea"`  | —                          | —                        | Map regions (maps)   |

**Common patterns:**

```typescript
// Table: indicators across columns, regions down rows
disaggregateBy: [
  { disOpt: "indicator_common_id", disDisplayOpt: "col" },
  { disOpt: "admin_area_2", disDisplayOpt: "row" },
]

// Timeseries: multiple lines per indicator
disaggregateBy: [
  { disOpt: "indicator_common_id", disDisplayOpt: "series" },
]

// Small multiples: separate chart per region
disaggregateBy: [
  { disOpt: "admin_area_2", disDisplayOpt: "replicant" },
  { disOpt: "indicator_common_id", disDisplayOpt: "series" },
]
```

---

## Post-Aggregation

### postAggregationExpression

Computes a derived value AFTER the main aggregation. Use when you need to
combine multiple columns with different aggregation functions.

```typescript
postAggregationExpression: {
  ingredientValues: [
    { prop: "cases", func: "SUM" },
    { prop: "tests", func: "SUM" },
  ],
  expression: "cases / NULLIF(tests, 0)",
}
```

**How it works:**
1. SQL aggregates each ingredientValue with its own func
2. Results wrapped in subquery
3. Outer query applies the expression

**Result SQL:**
```sql
SELECT admin_area_2, (cases / NULLIF(tests, 0)) AS result
FROM (
  SELECT admin_area_2, SUM(cases) AS cases, SUM(tests) AS tests
  FROM results GROUP BY admin_area_2
) AS subq
```

**When to use:**
- Ratios: `numerator / NULLIF(denominator, 0)`
- Differences: `col_a - col_b`
- Weighted averages
- Any calculation that needs columns aggregated differently

**Note:** When using `postAggregationExpression`, set `valueProps` to the output
column name(s) you want to display, and `valueFunc` typically to `"identity"`.

---

## vizPresets

Each vizPreset defines a default visualization for the metric.

### Structure

```typescript
vizPresets: [
  {
    id: "outlier-table",                    // Unique within metric
    label: { en: "...", fr: "..." },        // Shown in preset picker
    description: { en: "...", fr: "..." },  // Tooltip
    importantNotes: null,                   // Or { en: "...", fr: "..." }
    needsReplicant: false,                  // True if requires replicant selection
    allowedFilters: ["indicator_common_id", "admin_area_2"],
    createDefaultVisualizationOnInstall: "uuid-here",  // Or null
    config: {
      d: { /* data config */ },
      s: { /* style config */ },
      t: { /* text config */ },
    },
  },
]
```

### config.d (Data Config)

Controls what data to fetch and how to arrange it:

```typescript
d: {
  type: "table",                           // "table" | "timeseries" | "chart" | "map"
  timeseriesGrouping: "period_id",         // "period_id" | "quarter_id" | "year"
  valuesDisDisplayOpt: "col",              // Where values appear
  valuesFilter: ["cases"],                 // Optional: subset of valueProps
  disaggregateBy: [
    { disOpt: "admin_area_2", disDisplayOpt: "row" },
  ],
  filterBy: [                              // Fixed filters
    { disOpt: "facility_type", values: ["hospital"] },
  ],
  periodFilter: { filterType: "last_n_months", nMonths: 12 },
}
```

### config.s (Style Config)

Visual styling — all fields optional, merged with defaults:

```typescript
s: {
  content: "lines",                        // "lines" | "bars" | "points" | "areas"
  decimalPlaces: 1,                        // 0 | 1 | 2 | 3
  colorScale: "pastel-discrete",
  hideLegend: false,
  showDataLabels: false,
  barsStacked: false,
  forceYMax1: true,                        // Force Y-axis to 100% for percentages
  // Conditional formatting (scale-based)
  cfMode: "scale",                         // "none" | "scale" | "thresholds"
  cfScalePaletteKind: "custom",            // "preset" | "custom"
  cfScaleCustomFrom: "#fee0d2",
  cfScaleCustomTo: "#de2d26",
  // Conditional formatting (threshold-based)
  // cfMode: "thresholds",
  // cfThresholdCutoffs: [0.9, 0.8, 0.7],
  // cfThresholdBuckets: [{ color: "#22c55e" }, { color: "#eab308" }, ...],
}
```

### config.t (Text Config)

Captions and footnotes:

```typescript
t: {
  caption: { en: "Title", fr: "Titre" },
  captionRelFontSize: null,                // Or number for custom size
  subCaption: { en: "Subtitle DATE_RANGE", fr: "..." },
  subCaptionRelFontSize: null,
  footnote: { en: "Methodology note...", fr: "..." },
  footnoteRelFontSize: null,
}
```

**Note:** `DATE_RANGE` is a magic string replaced with the actual date range at render time.

---

## Other Fields

### valueLabelReplacements

Maps column names to display-friendly labels:

```typescript
valueLabelReplacements: {
  outlier_flag: "Outlier",
  completeness_flag: "Complete",
}
```

Used in chart legends and table headers.

### variantLabel

When multiple metrics share the same `label`, each must have a `variantLabel`
to distinguish them in lists:

```typescript
// Two metrics with same label
{ id: "m1-01-01", label: { en: "Coverage", ... }, variantLabel: { en: "DTP3", ... } }
{ id: "m1-01-02", label: { en: "Coverage", ... }, variantLabel: { en: "MCV1", ... } }
```

### hide

If `true`, metric is computed but not shown in the UI. Useful for intermediate
calculations or data that powers other metrics.

```typescript
hide: true
```

### aiDescription

Structured metadata for humans and AI systems. All fields are translatable strings:

| Field                    | Purpose                                      |
|--------------------------|----------------------------------------------|
| `summary`                | One-sentence description                     |
| `methodology`            | How it's calculated                          |
| `interpretation`         | What values mean, what's good/bad            |
| `typicalRange`           | Expected value range                         |
| `caveats`                | Important limitations (nullable)             |
| `disaggregationGuidance` | How to best disaggregate this metric         |

---

## Valid Disaggregation Options

From `disaggregation_options.ts`:

| Option               | Description                    |
|----------------------|--------------------------------|
| `period_id`          | Month (YYYYMM format)          |
| `quarter_id`         | Quarter (YYYYQ format)         |
| `year`               | Year                           |
| `indicator_common_id`| Standardized indicator code    |
| `admin_area_1`       | First admin level (province)   |
| `admin_area_2`       | Second admin level (district)  |
| `admin_area_3`       | Third admin level (subdistrict)|
| `facility_id`        | Individual facility            |
| `facility_type`      | Type of facility               |
| `facility_ownership` | Public/private/NGO etc         |
| `urbanicity`         | Urban/rural classification     |

---

## Common Patterns

### Simple proportion metric

```typescript
{
  id: "m1-01-01",
  resultsObjectId: "M1_output_outliers.csv",
  valueProps: ["outlier_flag"],           // Binary 0/1 column
  valueFunc: "AVG",                        // AVG of binary = proportion
  formatAs: "percent",
  label: { en: "Proportion of outliers", fr: "..." },
  // ...
}
```

### Count metric

```typescript
{
  id: "m1-01-00",
  resultsObjectId: "M1_output_outliers.csv",
  valueProps: ["facility_id"],
  valueFunc: "COUNT",
  formatAs: "number",
  label: { en: "Number of records", fr: "..." },
  hide: true,                              // Often hidden, used for context
  // ...
}
```

### Ratio metric with postAggregation

```typescript
{
  id: "m2-01-01",
  resultsObjectId: "M2_output_testing.csv",
  valueProps: ["positivity_rate"],         // Output column name
  valueFunc: "identity",                   // Don't re-aggregate
  formatAs: "percent",
  postAggregationExpression: {
    ingredientValues: [
      { prop: "positive_cases", func: "SUM" },
      { prop: "total_tests", func: "SUM" },
    ],
    expression: "positive_cases / NULLIF(total_tests, 0)",
  },
  // ...
}
```

---

## Validation Checklist

Before finalizing a metric:

- [ ] `resultsObjectId` matches an ID in `_results_objects.ts`
- [ ] `valueProps` columns exist in the results object's `createTableStatementPossibleColumns`
- [ ] `valuesFilter` (if used) only includes values from `valueProps`
- [ ] `postAggregationExpression.ingredientValues[].prop` columns exist in results object
- [ ] Metrics with duplicate `label` all have `variantLabel`
- [ ] `allowedFilters` uses valid disaggregation options
- [ ] `disaggregateBy[].disOpt` uses valid disaggregation options
