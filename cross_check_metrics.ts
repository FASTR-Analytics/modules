/**
 * Cross-checks metric definitions for internal consistency.
 *
 * Validates relationships between metrics, results objects, and vizPresets
 * that the Zod schema cannot express.
 *
 * Run: deno task crosscheck
 */

import { walk } from "jsr:@std/fs@^1/walk";

const ROOT = new URL(".", import.meta.url).pathname;

interface ResultsObject {
  id: string;
  createTableStatementPossibleColumns: Record<string, string> | false;
}

interface Metric {
  id: string;
  resultsObjectId: string;
  valueProps: string[];
  postAggregationExpression: {
    ingredientValues: { prop: string; func: string }[];
    expression: string;
  } | null;
  vizPresets: {
    id: string;
    config: {
      d: {
        valuesFilter?: string[];
        disaggregateBy?: { disOpt: string }[];
      };
    };
  }[];
}

interface ModuleDefinition {
  resultsObjects: ResultsObject[];
  metrics: Metric[];
}

interface Issue {
  module: string;
  metricId: string;
  presetId?: string;
  message: string;
}

function crossCheckModule(modulePath: string, def: ModuleDefinition): Issue[] {
  const issues: Issue[] = [];
  const moduleName = modulePath.split("/").at(-2) ?? modulePath;

  const resultsObjectColumns = new Map<string, Set<string>>();
  for (const ro of def.resultsObjects) {
    if (ro.createTableStatementPossibleColumns !== false) {
      resultsObjectColumns.set(
        ro.id,
        new Set(Object.keys(ro.createTableStatementPossibleColumns))
      );
    }
  }

  const allMetricIds = new Set(def.metrics.map((m) => m.id));

  for (const metric of def.metrics) {
    const roColumns = resultsObjectColumns.get(metric.resultsObjectId);

    // Check 1: valueProps reference valid columns
    // Skip if postAggregationExpression is used - valueProps are computed output names
    if (roColumns && !metric.postAggregationExpression) {
      for (const prop of metric.valueProps) {
        if (!roColumns.has(prop)) {
          issues.push({
            module: moduleName,
            metricId: metric.id,
            message: `valueProps "${prop}" not found in resultsObject "${metric.resultsObjectId}" columns: [${[...roColumns].join(", ")}]`,
          });
        }
      }
    }

    // Check 2: postAggregationExpression ingredientValues reference valid columns
    if (metric.postAggregationExpression && roColumns) {
      for (const ingredient of metric.postAggregationExpression.ingredientValues) {
        if (!roColumns.has(ingredient.prop)) {
          issues.push({
            module: moduleName,
            metricId: metric.id,
            message: `postAggregationExpression ingredientValue "${ingredient.prop}" not found in resultsObject "${metric.resultsObjectId}" columns`,
          });
        }
      }
    }

    // Check vizPresets
    for (const preset of metric.vizPresets) {
      // Check 4: valuesFilter references valid valueProps
      const valuesFilter = preset.config.d.valuesFilter;
      if (valuesFilter && valuesFilter.length > 0) {
        const valuePropSet = new Set(metric.valueProps);
        for (const filterValue of valuesFilter) {
          if (!valuePropSet.has(filterValue)) {
            issues.push({
              module: moduleName,
              metricId: metric.id,
              presetId: preset.id,
              message: `valuesFilter "${filterValue}" not found in metric valueProps: [${metric.valueProps.join(", ")}]`,
            });
          }
        }
      }

      // Check 5: disaggregateBy disOpt references columns that exist in results object
      // (Only check non-standard disaggregation options that should be columns)
      const standardDisOpts = new Set([
        "period_id",
        "quarter_id",
        "year",
        "indicator_common_id",
        "admin_area_1",
        "admin_area_2",
        "admin_area_3",
        "facility_id",
        "facility_type",
        "facility_ownership",
        "urbanicity",
      ]);

      const disaggregateBy = preset.config.d.disaggregateBy;
      if (disaggregateBy && roColumns) {
        for (const dis of disaggregateBy) {
          if (!standardDisOpts.has(dis.disOpt) && !roColumns.has(dis.disOpt)) {
            issues.push({
              module: moduleName,
              metricId: metric.id,
              presetId: preset.id,
              message: `disaggregateBy disOpt "${dis.disOpt}" is not a standard option and not found in resultsObject columns`,
            });
          }
        }
      }
    }
  }

  return issues;
}

async function main() {
  let totalIssues = 0;
  let modulesChecked = 0;

  console.log("Cross-checking metric definitions...\n");

  for await (const entry of walk(ROOT, {
    match: [/\/definition\.json$/],
    skip: [/\/\./, /\/node_modules\//],
  })) {
    const relativePath = entry.path.replace(ROOT, "");

    let def: ModuleDefinition;
    try {
      def = JSON.parse(await Deno.readTextFile(entry.path));
    } catch (e) {
      console.error(`ERROR  ${relativePath}: Failed to parse JSON: ${e}`);
      continue;
    }

    if (!def.metrics || !def.resultsObjects) {
      continue;
    }

    const issues = crossCheckModule(entry.path, def);
    modulesChecked++;

    if (issues.length === 0) {
      console.log(`OK     ${relativePath}`);
    } else {
      console.log(`ISSUES ${relativePath}`);
      for (const issue of issues) {
        const location = issue.presetId
          ? `${issue.metricId} > ${issue.presetId}`
          : issue.metricId;
        console.log(`       [${location}] ${issue.message}`);
      }
      totalIssues += issues.length;
    }
  }

  console.log("");
  console.log(`Checked ${modulesChecked} modules.`);

  if (totalIssues > 0) {
    console.log(`Found ${totalIssues} issue(s).`);
    Deno.exit(1);
  } else {
    console.log("All cross-checks passed.");
  }
}

main();
