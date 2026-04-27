# WB-FASTR Modules

This repository contains the source definitions for analysis modules used by
**wb-fastr** (World Bank FASTR - Frequent Assessments and System Tools for
Resilience). Each module defines an R-based data processing pipeline with
configurable parameters, metrics, and output schemas.

---

## Relationship to wb-fastr

```
wb-fastr-modules (this repo)          wb-fastr (main app)
┌─────────────────────────┐           ┌─────────────────────────────────┐
│ m001/                   │           │ server/module_loader/           │
│   _core.ts              │  build    │   load_module.ts                │
│   _metrics.ts           │  ──────►  │     - fetches definition.json   │
│   _parameters.ts        │           │     - fetches script.R          │
│   _results_objects.ts   │           │     - validates against schema  │
│   script.R              │           │                                 │
│   definition.json ◄─────┼───────────┼── GitHub raw fetch              │
└─────────────────────────┘           │                                 │
                                      │ server/worker_routines/         │
                                      │   run_module/                   │
                                      │     - injects parameters        │
                                      │     - executes R in Docker      │
                                      │     - imports CSV results       │
                                      └─────────────────────────────────┘
```

**Flow:**
1. Module authors edit TypeScript definitions and R scripts in this repo
2. Build script compiles TypeScript into `definition.json`
3. Push to GitHub (FASTR-Analytics/modules)
4. wb-fastr fetches `definition.json` + `script.R` from GitHub on demand
5. User installs module in a project → definition stored in PostgreSQL
6. User runs module → R script executed, CSV output imported to database

---

## Module Structure

Each module lives in a directory named `m###` (e.g., `m001`, `m002`). Required files:

| File                   | Purpose                                           |
|------------------------|---------------------------------------------------|
| `_core.ts`             | Label, prerequisites, data sources, script type   |
| `_metrics.ts`          | Metric definitions with aggregation and viz       |
| `_parameters.ts`       | User-configurable parameters injected into R      |
| `_results_objects.ts`  | Output CSV schemas (column names and SQL types)   |
| `script.R`             | The R script that performs analysis               |
| `definition.json`      | **Generated** - do not edit directly              |

### _core.ts

```typescript
import type { ModuleDefinitionCore } from "../.validation/_module_definition_github.ts";

export const core: ModuleDefinitionCore = {
  label: {
    en: "M1. Data quality assessment",
    fr: "M1. Évaluation de la qualité des données",
  },
  prerequisites: [],                    // Module IDs that must run first
  scriptGenerationType: "template",     // "template" or "hfa"
  dataSources: [
    {
      sourceType: "dataset",
      replacementString: "PROJECT_DATA_HMIS",  // Replaced in R script
      datasetType: "hmis",
    },
  ],
  assetsToImport: [],                   // Additional files copied to sandbox
};
```

### _metrics.ts

Metrics define how results are aggregated and visualized. Each metric references
a results object and specifies which columns to aggregate.

```typescript
import type { MetricDefinitionGithub } from "../.validation/_module_definition_github.ts";

export const metrics: MetricDefinitionGithub[] = [
  {
    id: "m1-01-01",
    resultsObjectId: "M1_output_outliers.csv",
    valueProps: ["outlier_flag"],
    valueFunc: "AVG",                   // COUNT, SUM, AVG, MIN, MAX, identity
    formatAs: "percent",
    label: { en: "Proportion of outliers", fr: "Proportion de valeurs aberrantes" },
    variantLabel: null,
    requiredDisaggregationOptions: [],
    valueLabelReplacements: {},
    postAggregationExpression: null,
    aiDescription: null,
    importantNotes: null,
    hide: false,
    vizPresets: [
      {
        id: "outlier-timeseries",
        label: { en: "Trend over time", fr: "Tendance dans le temps" },
        description: { en: "Line chart by indicator", fr: "Graphique linéaire par indicateur" },
        importantNotes: null,
        needsReplicant: false,
        allowedFilters: ["indicator_common_id"],
        createDefaultVisualizationOnInstall: null,
        config: {
          d: {
            type: "timeseries",
            timeseriesGrouping: "period_id",
            valuesDisDisplayOpt: "series",
            disaggregateBy: [{ disOpt: "indicator_common_id", disDisplayOpt: "series" }],
            filterBy: [],
          },
          s: { content: "lines" },
          t: { caption: null, captionRelFontSize: null, subCaption: null, subCaptionRelFontSize: null, footnote: null, footnoteRelFontSize: null },
        },
      },
    ],
  },
];
```

### _parameters.ts

Parameters are injected into the R script via string replacement before execution.

```typescript
import type { ModuleParameter } from "../.validation/_module_definition_github.ts";

export const parameters: ModuleParameter[] = [
  {
    replacementString: "OUTLIER_THRESHOLD",
    description: {
      en: "Threshold for outlier detection",
      fr: "Seuil de détection des valeurs aberrantes",
    },
    input: {
      inputType: "number",
      defaultValue: "3",
    },
  },
];
```

In the R script, use the replacement string as a placeholder:
```r
threshold <- OUTLIER_THRESHOLD
```

### _results_objects.ts

Defines the schema for each CSV file the R script produces. Column names must
match exactly what the R script writes.

```typescript
import type { ResultsObjectDefinitionGithub } from "../.validation/_module_definition_github.ts";

export const resultsObjects: ResultsObjectDefinitionGithub[] = [
  {
    id: "M1_output_outliers.csv",
    createTableStatementPossibleColumns: {
      facility_id: "TEXT NOT NULL",
      period_id: "INTEGER NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      outlier_flag: "INTEGER NOT NULL",
    },
  },
];
```

### script.R

The R script receives:
- Data source tables via replacement strings (e.g., `PROJECT_DATA_HMIS`)
- Parameter values via replacement strings
- Working directory set to the module sandbox

Must output CSV files matching the results objects:
```r
write.csv(result_df, "M1_output_outliers.csv", row.names = FALSE)
```

---

## Build Process

### Building Definitions

```bash
deno run --allow-read --allow-write build_definitions.ts
```

This:
1. Scans for directories matching `m\d+`
2. Imports all four TypeScript files from each module
3. Validates against Zod schema in `.validation/`
4. Writes combined `definition.json` to each module directory

### Validating Definitions

```bash
deno run --allow-read --allow-net .validation/validate_definitions.ts
```

### Deploy Script

The `./deploy` script runs build + validation, then commits and pushes:

```bash
./deploy
```

It will prompt for a commit message and warn if validation fails.

---

## Adding a New Module

1. Create directory `m###/` with the next available number
2. Create all four TypeScript files following the patterns above
3. Create `script.R` with placeholder replacement strings
4. Run `deno run --allow-read --allow-write build_definitions.ts`
5. Verify `definition.json` was generated
6. Test locally by setting `MODULE_SOURCE=local` in wb-fastr

---

## Module Registry

wb-fastr maintains a registry of available modules in
`lib/types/module_registry.ts`. After adding a module here, it must also be
registered there to appear in the UI.

---

## Local Development

To test modules without pushing to GitHub:

1. In wb-fastr, set environment variable `MODULE_SOURCE=local`
2. Place modules in wb-fastr's `modules/` directory (or symlink)
3. wb-fastr will read from filesystem instead of GitHub

---

## Schema Validation

The `.validation/` directory contains:

| File                              | Purpose                                    |
|-----------------------------------|--------------------------------------------|
| `_module_definition_github.ts`    | Zod schema for module definitions          |
| `validate_definitions.ts`         | Script to validate all definition.json     |
| `cf_presets.ts`                   | Conditional formatting preset definitions  |
| `disaggregation_options.ts`       | Valid disaggregation dimension options     |

The schema here must stay in sync with wb-fastr's
`lib/types/_module_definition_github.ts`. When adding new fields:
1. Update both schemas
2. Rebuild all definitions
3. Test in wb-fastr
