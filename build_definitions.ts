import { ModuleDefinitionJSONSchema } from "./.validation/module_definition_validator.ts";

const root = new URL(".", import.meta.url);
const moduleDirs: string[] = [];
for await (const entry of Deno.readDir(root)) {
  if (entry.isDirectory && /^m\d+$/.test(entry.name)) {
    moduleDirs.push(entry.name);
  }
}
moduleDirs.sort();

let hadError = false;

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
}

if (hadError) Deno.exit(1);
