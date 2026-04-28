import { metrics as m007Metrics } from "./m007/_metrics.ts";
import type { MetricDefinitionGithub, VizPreset } from "./.validation/_module_definition_github.ts";

const modules: { dir: string; metrics: MetricDefinitionGithub[] }[] = [
  { dir: "m007", metrics: m007Metrics },
];

function stringifyValue(value: unknown, indent: number): string {
  const spaces = "  ".repeat(indent);
  const innerSpaces = "  ".repeat(indent + 1);

  if (value === null) return "null";
  if (value === undefined) return "undefined";
  if (typeof value === "string") return JSON.stringify(value);
  if (typeof value === "number" || typeof value === "boolean") return String(value);

  if (Array.isArray(value)) {
    if (value.length === 0) return "[]";
    const items = value.map(v => stringifyValue(v, indent + 1));
    if (items.every(i => !i.includes("\n")) && items.join(", ").length < 60) {
      return `[${items.join(", ")}]`;
    }
    return `[\n${innerSpaces}${items.join(`,\n${innerSpaces}`)},\n${spaces}]`;
  }

  if (typeof value === "object") {
    const entries = Object.entries(value);
    if (entries.length === 0) return "{}";
    const lines = entries.map(([k, v]) => {
      const key = /^[a-zA-Z_][a-zA-Z0-9_]*$/.test(k) ? k : JSON.stringify(k);
      return `${key}: ${stringifyValue(v, indent + 1)}`;
    });
    return `{\n${innerSpaces}${lines.join(`,\n${innerSpaces}`)},\n${spaces}}`;
  }

  return String(value);
}

function generateMetricFile(metric: MetricDefinitionGithub): string {
  const lines: string[] = [];

  lines.push(`import type { MetricDefinitionGithub, VizPreset } from "../../.validation/_module_definition_github.ts";`);
  lines.push(``);

  // Generate vizPresets as separate const
  lines.push(`export const vizPresets: VizPreset[] = ${stringifyValue(metric.vizPresets, 0)};`);
  lines.push(``);

  // Generate metric without vizPresets, then add reference
  const metricWithoutVizPresets = { ...metric };
  delete (metricWithoutVizPresets as Record<string, unknown>).vizPresets;

  lines.push(`export const metric: MetricDefinitionGithub = {`);
  for (const [key, value] of Object.entries(metricWithoutVizPresets)) {
    lines.push(`  ${key}: ${stringifyValue(value, 1)},`);
  }
  lines.push(`  vizPresets,`);
  lines.push(`};`);
  lines.push(``);

  return lines.join("\n");
}

for (const { dir, metrics } of modules) {
  const metricsDir = `./${dir}/_metrics`;

  try {
    await Deno.mkdir(metricsDir, { recursive: true });
  } catch {
    // Already exists
  }

  for (const metric of metrics) {
    const filename = `${metricsDir}/${metric.id}.ts`;
    const content = generateMetricFile(metric);
    await Deno.writeTextFile(filename, content);
    console.log(`wrote ${filename}`);
  }
}

console.log("\nDone. You can now delete the old _metrics.ts files.");
