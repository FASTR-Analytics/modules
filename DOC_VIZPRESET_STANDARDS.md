# VizPreset Standards

Guidelines for writing consistent, clean viz presets across modules.

## IDs

### Viz Preset ID (`id`)

**Scope:** Only needs to be unique within a single metric's `vizPresets` array. Lookups are metric-scoped (`metric.vizPresets.find(p => p.id === presetId)`).

**Pattern:** `{metric-concept}-{viz-type}[-{parent-context}]`

- Use kebab-case
- Include viz type: `table`, `map`, `chart`, `timeseries`
- Include parent context when there's a replicant filter (e.g., `single-admin-area-2-multiple-admin-area-3`)
- Base case (no parent filter) can omit the parent context

**Examples:**
```
disruption-differences-table                                          # base table, no replicant
disruption-differences-table-single-admin-area-2-multiple-admin-area-3  # table with AA2 replicant
disruption-differences-map-admin-area-2                               # map showing AA2
disruption-chart-single-indicator                                     # chart filtered to single indicator
```

### createDefaultVisualizationOnInstall

**Scope:** Must be globally unique across all modules. This UUID becomes the presentation object's primary key in the database.

- Always use a valid UUID v4
- Never use placeholders like `"MAKE THIS FOR ME"`
- Set to `null` if the preset should not create a default visualization on install

## Labels

**Purpose:** Short, scannable text for card UI (192px min width).

**Pattern:** `{Metric concept} [map] by Admin Area {N}`

- Keep concise - labels appear in small preview cards
- Include admin area level for differentiation
- Add "map" for map visualizations
- No parentheticals or verbose qualifiers

**Examples:**
```
"Service disruption by Admin Area 2"
"Service disruption map by Admin Area 2"
"Service disruption by Admin Area 3"
```

## Descriptions

**Purpose:** Verbose context shown below the label in card UI.

**Pattern:** `{Viz type} showing {what} [for {parent context}], with {dimensions}`

- Be specific about what's displayed
- Include parent/child context (single/multiple admin areas)
- Use "actual vs expected" for difference comparisons
- Do NOT mention implementation details like "conditional formatting"

**Examples:**
```
"Table showing percentage difference between actual vs expected service volume, with multiple Admin Areas 2 and multiple indicators"
"Table showing percentage difference between actual vs expected service volume for a single Admin Area 2, with multiple Admin Areas 3 and multiple indicators"
"Map showing percentage difference between actual vs expected service volume by Admin Area 2"
```

## Config Structure

### d (data config)

```typescript
d: {
  type: "table" | "map" | "timeseries",
  timeseriesGrouping: "period_id",
  valuesDisDisplayOpt: "col" | "cell",
  disaggregateBy: [...],
  selectedReplicantValue: "",  // Required when using replicant
  filterBy: [],
  periodFilter: {
    filterType: "last_n_calendar_quarters",
    nQuarters: 1,
  },
}
```

**Key rules:**
- Always include `selectedReplicantValue` when any disaggregation uses `disDisplayOpt: "replicant"`
- Use `""` for tables (user selects), `"anc1"` or specific value for maps
- Snapshot visualizations use `nQuarters: 1`

### s (style config)

```typescript
s: {
  ...CF_DIVERGING_5_10_20,  // Use presets from cf_presets.ts
}
```

**Available presets:**
- `CF_DIVERGING_5_10_20` - 7-bucket diverging (±5%, ±10%, ±20%) - matches app preset "fmt-thresholds-5-10-20"
- `CF_NEG10_POS10` - 3-bucket diverging (±10%)
- `CF_90_80`, `CF_80_70` - Coverage thresholds (higher is better)
- `CF_01_03`, `CF_05_10`, `CF_10_20` - Rate thresholds (lower is better)

### t (text config)

**subCaption must always include date range:**

- English: Use `DATE_RANGE` placeholder
- French: Use `PLAGE_DE_DATES` placeholder

Minimal subCaption (all presets should have at least this):

```typescript
subCaption: {
  en: "DATE_RANGE",
  fr: "PLAGE_DE_DATES",
}
```

For coverage metrics with disclaimer:
```typescript
import { SUBCAPTION_DISCLAIMER } from "../../_shared/text_presets.ts";
// ...
subCaption: SUBCAPTION_DISCLAIMER,
```

For maps with descriptive subCaption:

```typescript
subCaption: {
  en: "Percentage difference between actual and expected service volume by Admin Area 2, DATE_RANGE",
  fr: "Différence en pourcentage entre le volume de services réel et attendu par Zone administrative 2, PLAGE_DE_DATES",
}
```

## Text Conventions

### Terminology
- Use "Admin Area 2/3/4" consistently, not "region", "district", "sub-district"
- Use "actual vs expected" for difference/comparison metrics
- Use "service volume" not "service delivery volume"

### Footnotes
- Describe data interpretation, not UI implementation
- Do NOT mention colors (e.g., "darker colors indicate...")
- Focus on what values mean

**Good:**
```
"Negative values indicate service delivery below expected levels. Expected values are based on pre-disruption trends."
```

**Bad:**
```
"Darker colors indicate larger disruptions."
```

### French translations
- "Admin Area N" → "Zone administrative N"
- "Service disruption" → "Perturbation des services"
- "percentage difference" → "différence en pourcentage"
- "actual vs expected" → "réel et attendu"

## Symmetry Checklist

When reviewing a set of related metrics (e.g., m3-03-02, m3-04-02, m3-05-02):

- [ ] IDs follow consistent pattern across admin levels
- [ ] Labels follow same structure, differing only by admin level
- [ ] Descriptions follow same structure
- [ ] All use same CF preset
- [ ] All use same periodFilter settings
- [ ] All tables with replicants have `selectedReplicantValue`
- [ ] All maps have caption/subCaption/footnote (or all don't)
- [ ] Map subCaptions use consistent "Admin Area N" terminology
- [ ] Footnotes don't reference UI implementation details

## Shared Text Presets

Common text blocks are defined in `_shared/text_presets.ts` to avoid duplication:

| Constant                               | Usage                                                    |
| -------------------------------------- | -------------------------------------------------------- |
| `SUBCAPTION_DISCLAIMER`                | Coverage metrics subCaption with DATE_RANGE + disclaimer |
| `FOOTNOTE_COVERAGE_BASE`               | Coverage footnote with MICS/DHS attribution              |
| `FOOTNOTE_COVERAGE_WITH_CURRENT_YEAR`  | Extended footnote with current year scaling note         |

Import and use:

```typescript
import {
  FOOTNOTE_COVERAGE_WITH_CURRENT_YEAR,
  SUBCAPTION_DISCLAIMER,
} from "../../_shared/text_presets.ts";

// In vizPreset config:
t: {
  caption: { en: "...", fr: "..." },
  subCaption: SUBCAPTION_DISCLAIMER,
  footnote: FOOTNOTE_COVERAGE_WITH_CURRENT_YEAR,
  // ...
}
```

## Adding New CF Presets

When adding a new conditional formatting preset to `.validation/cf_presets.ts`:

1. Check if an equivalent exists in `wb-fastr/lib/legacy_cf_presets.ts`
2. If yes, use exact same colors to ensure UI shows it as a preset (not "Custom")
3. Use `{ key: "base200" }` for neutral colors when matching app presets
4. Document which app preset it matches in a comment

```typescript
// Diverging 7-bucket: matches app preset "fmt-thresholds-5-10-20"
export const CF_DIVERGING_5_10_20 = thresholds(
  [-0.2, -0.1, -0.05, 0.05, 0.1, 0.2],
  [
    { color: "#e73535" },
    { color: CF_LIGHTER_RED },
    { color: "#f8c4c4" },
    { color: { key: "base200" } },
    { color: "#b4e3c8" },
    { color: CF_LIGHTER_GREEN },
    { color: "#3ea46a" },
  ],
);
```
