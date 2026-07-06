import {
  ModuleDefinitionJSONSchema,
  type ModuleDefinitionGithub,
  type MetricDefinitionGithub,
} from "./.validation/_module_definition_github.ts";

const root = new URL(".", import.meta.url);
const moduleDirs: string[] = [];
for await (const entry of Deno.readDir(root)) {
  if (entry.isDirectory && /^m\d+$/.test(entry.name)) {
    moduleDirs.push(entry.name);
  }
}
moduleDirs.sort();

async function loadMetrics(dir: string): Promise<{ metrics: MetricDefinitionGithub[]; errors: string[] }> {
  const metricsFolderPath = new URL(`${dir}/_metrics/`, root);
  try {
    await Deno.stat(metricsFolderPath);
  } catch {
    const metricsMod = await import(`./${dir}/_metrics.ts`);
    return { metrics: metricsMod.metrics, errors: [] };
  }
  const metrics: MetricDefinitionGithub[] = [];
  const errors: string[] = [];
  const files: string[] = [];
  for await (const entry of Deno.readDir(metricsFolderPath)) {
    if (entry.isFile && entry.name.endsWith(".ts")) {
      files.push(entry.name);
    }
  }
  files.sort();
  for (const file of files) {
    const expectedId = file.replace(/\.ts$/, "");
    const mod = await import(`./${dir}/_metrics/${file}`);
    if (mod.metric.id !== expectedId) {
      errors.push(`${dir}/_metrics/${file}: metric.id "${mod.metric.id}" does not match filename "${expectedId}"`);
    }
    metrics.push(mod.metric);
  }
  return { metrics, errors };
}

let hadError = false;
const allModules: { dir: string; def: ModuleDefinitionGithub }[] = [];

for (const dir of moduleDirs) {
  const corePath = `./${dir}/_core.ts`;
  const resultsPath = `./${dir}/_results_objects.ts`;
  const parametersPath = `./${dir}/_parameters.ts`;

  try {
    await Deno.stat(new URL(`${dir}/_core.ts`, root));
  } catch {
    console.log(`skip  ${dir} (no _core.ts)`);
    continue;
  }

  const [coreMod, metricsResult, resultsMod, parametersMod] = await Promise.all([
    import(corePath),
    loadMetrics(dir),
    import(resultsPath),
    import(parametersPath),
  ]);

  if (metricsResult.errors.length > 0) {
    hadError = true;
    for (const err of metricsResult.errors) {
      console.error(`FAIL  ${err}`);
    }
    continue;
  }

  const definition = {
    ...coreMod.core,
    resultsObjects: resultsMod.resultsObjects,
    metrics: metricsResult.metrics,
    configRequirements: { parameters: parametersMod.parameters },
  };

  const result = ModuleDefinitionJSONSchema.safeParse(definition);
  if (!result.success) {
    hadError = true;
    console.error(`FAIL  ${dir}`);
    for (const issue of result.error.issues) {
      console.error(`  - ${issue.path.join(".")}: ${issue.message}`);
    }
    continue;
  }

  const outPath = new URL(`${dir}/definition.json`, root);
  await Deno.writeTextFile(
    outPath,
    JSON.stringify(result.data, null, 2) + "\n",
  );
  console.log(`build ${dir}`);
  allModules.push({ dir, def: result.data });
}

if (hadError) Deno.exit(1);

// Build column lookup for cross-checks
const resultsObjectColumns = new Map<string, Set<string>>();
for (const { def } of allModules) {
  for (const ro of def.resultsObjects) {
    if (ro.createTableStatementPossibleColumns !== false) {
      resultsObjectColumns.set(
        ro.id,
        new Set(Object.keys(ro.createTableStatementPossibleColumns))
      );
    }
  }
}

// Standard disaggregation options (not columns)
const standardDisOpts = new Set([
  "period_id",
  "quarter_id",
  "year",
  "indicator_common_id",
  "admin_area_1",
  "admin_area_2",
  "admin_area_3",
  "admin_area_4",
  "facility_id",
  "facility_type",
  "facility_ownership",
  "urbanicity",
]);

// Check for duplicate IDs and run all validations
const resultsObjectIds = new Map<string, string>();
const metricIds = new Map<string, string>();
const defaultVizIds = new Map<string, string>();

for (const { dir, def } of allModules) {
  for (const ro of def.resultsObjects) {
    if (resultsObjectIds.has(ro.id)) {
      console.error(`DUPLICATE resultsObjectId "${ro.id}" in ${dir} (first seen in ${resultsObjectIds.get(ro.id)})`);
      hadError = true;
    } else {
      resultsObjectIds.set(ro.id, dir);
    }
  }

  for (const metric of def.metrics) {
    if (metricIds.has(metric.id)) {
      console.error(`DUPLICATE metricId "${metric.id}" in ${dir} (first seen in ${metricIds.get(metric.id)})`);
      hadError = true;
    } else {
      metricIds.set(metric.id, dir);
    }

    const roColumns = resultsObjectColumns.get(metric.resultsObjectId);

    // Check valueProps reference valid columns (skip if postAggregationExpression)
    if (roColumns && !metric.postAggregationExpression) {
      for (const prop of metric.valueProps) {
        if (!roColumns.has(prop)) {
          console.error(`INVALID valueProps "${prop}" in ${dir}:${metric.id} - not in resultsObject columns`);
          hadError = true;
        }
      }
    }

    // Check postAggregationExpression ingredientValues reference valid columns
    if (metric.postAggregationExpression && roColumns) {
      for (const ingredient of metric.postAggregationExpression.ingredientValues) {
        if (!roColumns.has(ingredient.prop)) {
          console.error(`INVALID ingredientValue "${ingredient.prop}" in ${dir}:${metric.id} - not in resultsObject columns`);
          hadError = true;
        }
      }
    }

    const vizPresetIdsInMetric = new Set<string>();
    for (const preset of metric.vizPresets) {
      const location = `${dir}:${metric.id}:${preset.id}`;
      if (vizPresetIdsInMetric.has(preset.id)) {
        console.error(`DUPLICATE vizPresetId "${preset.id}" within metric ${metric.id} (${dir})`);
        hadError = true;
      } else {
        vizPresetIdsInMetric.add(preset.id);
      }

      // Check requiredDisaggregationOptions are in disaggregateBy
      const disaggregateByOpts = new Set(preset.config.d.disaggregateBy.map(d => d.disOpt));
      const timeOpts = new Set(["year", "period_id", "quarter_id"]);
      const isTimeseries = preset.config.d.type === "timeseries";
      const hasPeriodFilter = !!preset.config.d.periodFilter;
      for (const required of metric.requiredDisaggregationOptions) {
        if ((isTimeseries || hasPeriodFilter) && timeOpts.has(required)) continue;
        if (!disaggregateByOpts.has(required)) {
          console.error(`MISSING requiredDisaggregationOption "${required}" in ${location} disaggregateBy`);
          hadError = true;
        }
      }

      // Check replicant requires selectedReplicantValue (except for admin_area_* replicants)
      const replicantDis = preset.config.d.disaggregateBy.find(d => d.disDisplayOpt === "replicant");
      if (replicantDis && !replicantDis.disOpt.startsWith("admin_area_")) {
        if (!preset.config.d.selectedReplicantValue) {
          console.error(`MISSING selectedReplicantValue in ${location} (has replicant on ${replicantDis.disOpt})`);
          hadError = true;
        }
      }

      // Check valuesFilter references valid valueProps
      const valuesFilter = preset.config.d.valuesFilter;
      if (valuesFilter && valuesFilter.length > 0) {
        const valuePropSet = new Set(metric.valueProps);
        for (const filterValue of valuesFilter) {
          if (!valuePropSet.has(filterValue)) {
            console.error(`INVALID valuesFilter "${filterValue}" in ${location} - not in metric valueProps`);
            hadError = true;
          }
        }
      }

      // Check disaggregateBy disOpt references valid columns (non-standard options)
      if (roColumns) {
        for (const dis of preset.config.d.disaggregateBy) {
          if (!standardDisOpts.has(dis.disOpt) && !roColumns.has(dis.disOpt)) {
            console.error(`INVALID disaggregateBy "${dis.disOpt}" in ${location} - not standard and not in resultsObject columns`);
            hadError = true;
          }
        }
      }

      // Check createDefaultVisualizationOnInstall UUID
      if (preset.createDefaultVisualizationOnInstall) {
        const uuid = preset.createDefaultVisualizationOnInstall;
        if (!/^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/i.test(uuid)) {
          console.error(`INVALID UUID "${uuid}" in ${location} - not a valid UUID format`);
          hadError = true;
        }
        if (defaultVizIds.has(uuid)) {
          console.error(`DUPLICATE createDefaultVisualizationOnInstall "${uuid}" in ${location} (first seen in ${defaultVizIds.get(uuid)})`);
          hadError = true;
        } else {
          defaultVizIds.set(uuid, location);
        }
      }

      // Convention: subCaption must include DATE_RANGE/PLAGE_DE_DATES
      const subCaption = preset.config.t?.subCaption;
      if (subCaption && typeof subCaption === "object") {
        if (subCaption.en && !subCaption.en.includes("DATE_RANGE")) {
          console.error(`CONVENTION subCaption.en missing DATE_RANGE in ${location}`);
          hadError = true;
        }
        if (subCaption.fr && subCaption.fr.includes("DATE_RANGE")) {
          console.error(`CONVENTION subCaption.fr uses DATE_RANGE instead of PLAGE_DE_DATES in ${location}`);
          hadError = true;
        }
        if (subCaption.fr && !subCaption.fr.includes("PLAGE_DE_DATES")) {
          console.error(`CONVENTION subCaption.fr missing PLAGE_DE_DATES in ${location}`);
          hadError = true;
        }
        if (subCaption.pt && subCaption.pt.includes("DATE_RANGE")) {
          console.error(`CONVENTION subCaption.pt uses DATE_RANGE instead of INTERVALO_DE_DATAS in ${location}`);
          hadError = true;
        }
        if (subCaption.pt && !subCaption.pt.includes("INTERVALO_DE_DATAS")) {
          console.error(`CONVENTION subCaption.pt missing INTERVALO_DE_DATAS in ${location}`);
          hadError = true;
        }
      }

      // Convention: labels/descriptions should use "Admin Area N" not "region"/"district"
      const label = preset.label;
      if (label && typeof label === "object") {
        if (label.en && /\b(region|district)\b/i.test(label.en)) {
          console.error(`CONVENTION label.en uses "region" or "district" instead of "Admin Area N" in ${location}`);
          hadError = true;
        }
      }
      const description = preset.description;
      if (description && typeof description === "object") {
        if (description.en && /\b(regions|districts)\b/i.test(description.en)) {
          console.error(`CONVENTION description.en uses "regions" or "districts" instead of "Admin Areas N" in ${location}`);
          hadError = true;
        }
      }
    }
  }
}

if (hadError) {
  console.error("\nBuild failed with errors.");
  Deno.exit(1);
} else {
  console.log(`\nBuilt ${allModules.length} modules successfully.`);
}
