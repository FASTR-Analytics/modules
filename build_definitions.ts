import {
  ModuleDefinitionJSONSchema,
  type ModuleDefinitionGithub,
} from "./.validation/_module_definition_github.ts";

const root = new URL(".", import.meta.url);
const moduleDirs: string[] = [];
for await (const entry of Deno.readDir(root)) {
  if (entry.isDirectory && /^m\d+$/.test(entry.name)) {
    moduleDirs.push(entry.name);
  }
}
moduleDirs.sort();

let hadError = false;
const allModules: { dir: string; def: ModuleDefinitionGithub }[] = [];

for (const dir of moduleDirs) {
  const corePath = `./${dir}/_core.ts`;
  const metricsPath = `./${dir}/_metrics.ts`;
  const resultsPath = `./${dir}/_results_objects.ts`;
  const parametersPath = `./${dir}/_parameters.ts`;

  try {
    await Deno.stat(new URL(`${dir}/_core.ts`, root));
  } catch {
    console.log(`skip  ${dir} (no _core.ts)`);
    continue;
  }

  const [coreMod, metricsMod, resultsMod, parametersMod] = await Promise.all([
    import(corePath),
    import(metricsPath),
    import(resultsPath),
    import(parametersPath),
  ]);

  const definition = {
    ...coreMod.core,
    resultsObjects: resultsMod.resultsObjects,
    metrics: metricsMod.metrics,
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

// Check for duplicate IDs across all modules
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

    const vizPresetIdsInMetric = new Set<string>();
    for (const preset of metric.vizPresets) {
      const location = `${dir}:${metric.id}:${preset.id}`;
      if (vizPresetIdsInMetric.has(preset.id)) {
        console.error(`DUPLICATE vizPresetId "${preset.id}" within metric ${metric.id} (${dir})`);
        hadError = true;
      } else {
        vizPresetIdsInMetric.add(preset.id);
      }

      if (preset.createDefaultVisualizationOnInstall) {
        const uuid = preset.createDefaultVisualizationOnInstall;
        if (defaultVizIds.has(uuid)) {
          console.error(`DUPLICATE createDefaultVisualizationOnInstall "${uuid}" in ${location} (first seen in ${defaultVizIds.get(uuid)})`);
          hadError = true;
        } else {
          defaultVizIds.set(uuid, location);
        }
      }
    }
  }
}

if (hadError) Deno.exit(1);
