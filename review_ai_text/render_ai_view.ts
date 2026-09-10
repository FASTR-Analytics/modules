// Renders, from each module's built definition.json, the text an AI agent
// receives from wb-fastr's two shared metric tools (get_available_metrics and
// the header of get_metric_data) plus the SPA-only get_available_modules and
// get_module_settings tools. It mirrors the formatters in wb-fastr:
//   lib/ai_tools/format_metrics_list_for_ai.ts
//   lib/ai_tools/format_metric_data_for_ai.ts (formatItemsAsMarkdown, header)
//   client/src/components/project_ai/ai_tools/tools/_internal/format_modules_list_for_ai.ts
//   client/src/components/project_ai/ai_tools/tools/_internal/format_module_settings_for_ai.ts
// and server/runs/disaggregation_availability.ts for the disaggregation list.
// Re-check those files before trusting this output; the rules are copied, not
// imported. Facility options are shown as conditional because they depend on
// the instance's structure config.
//
// Usage: deno run --allow-read --allow-write review_ai_text/render_ai_view.ts <outDir>

import { FROZEN_MODULE_DIRS } from "../_frozen_modules.ts";
import type {
  MetricDefinitionGithub,
  ModuleDefinitionGithub,
  VizPreset,
} from "../.validation/_module_definition_github.ts";

const PHYSICAL_DISAGGREGATION_COLUMNS = [
  "admin_area_2",
  "admin_area_3",
  "admin_area_4",
  "indicator_common_id",
  "denominator",
  "denominator_best_or_survey",
  "source_indicator",
  "target_population",
  "ratio_type",
  "hfa_indicator",
  "hfa_variant_item",
  "hfa_category",
  "hfa_sub_category",
  "hfa_service_category",
  "time_point",
  "iceh_indicator",
  "strat",
  "level",
];
const FACILITY_OPTIONS = ["facility_type", "facility_ownership"];

type Lang = { en: string; fr: string; pt?: string };

const root = new URL("../", import.meta.url);
const outDir = Deno.args[0];
if (!outDir) {
  console.error("Usage: render_ai_view.ts <outDir>");
  Deno.exit(1);
}
await Deno.mkdir(outDir, { recursive: true });

type LoadedModule = { dir: string; def: ModuleDefinitionGithub };
const modules: LoadedModule[] = [];
for await (const entry of Deno.readDir(root)) {
  if (
    entry.isDirectory && /^m\d+$/.test(entry.name) &&
    !FROZEN_MODULE_DIRS.includes(entry.name)
  ) {
    const def = JSON.parse(
      await Deno.readTextFile(new URL(`${entry.name}/definition.json`, root)),
    ) as ModuleDefinitionGithub;
    modules.push({ dir: entry.name, def });
  }
}
modules.sort((a, b) => a.dir.localeCompare(b.dir));

function availableDisaggregations(
  def: ModuleDefinitionGithub,
  metric: MetricDefinitionGithub,
): { options: string[]; note: string } {
  const ro = def.resultsObjects.find((r) => r.id === metric.resultsObjectId);
  if (!ro || ro.createTableStatementPossibleColumns === false) {
    return {
      options: [],
      note: "results object declares no columns; availability cannot be derived here",
    };
  }
  const cols = new Set(Object.keys(ro.createTableStatementPossibleColumns));
  const out: string[] = [];
  for (const d of PHYSICAL_DISAGGREGATION_COLUMNS) if (cols.has(d)) out.push(d);
  let note = "";
  if (cols.has("facility_id")) {
    out.push(...FACILITY_OPTIONS.map((f) => `${f}?`));
    note = "options marked ? need the instance's structure config to enable them";
  }
  if (cols.has("period_id")) out.push("year", "month", "quarter_id", "period_id");
  else if (cols.has("quarter_id")) out.push("quarter_id", "year");
  else if (cols.has("year")) out.push("year");
  return { options: out, note };
}

function replicateByProp(preset: VizPreset): string | undefined {
  for (const dis of preset.config.d.disaggregateBy) {
    if (dis.disDisplayOpt !== "replicant") continue;
    const filter = preset.config.d.filterBy.find((f) => f.disOpt === dis.disOpt);
    if (filter && filter.values.length === 1) continue;
    return dis.disOpt;
  }
  return undefined;
}

function catalogEntry(def: ModuleDefinitionGithub, metric: MetricDefinitionGithub): string[] {
  const lines: string[] = [];
  const label = metric.variantLabel
    ? `${metric.label.en} [${metric.variantLabel.en}]`
    : metric.label.en;
  lines.push(`${metric.id}: ${label} [${metric.formatAs}]`);
  if (metric.aiDescription?.summary) lines.push(`  ${metric.aiDescription.summary.en}`);
  if (metric.importantNotes) lines.push(`  NOTE: ${metric.importantNotes.en}`);
  if (metric.valueProps.length > 0) {
    const props = metric.valueProps.map((p) => {
      const l = metric.valueLabelReplacements[p] || p;
      return l !== p ? `${p} (${l})` : p;
    });
    lines.push(`  Values: ${props.join(", ")}`);
  }
  const avail = availableDisaggregations(def, metric);
  const requiredSet = new Set<string>(metric.requiredDisaggregationOptions);
  const required = avail.options.filter((o) => requiredSet.has(o) && o !== "quarter_id");
  const optional = avail.options.filter((o) => !requiredSet.has(o) && o !== "quarter_id");
  if (required.length > 0) lines.push(`  Auto-disaggregated by: ${required.join(", ")}`);
  if (optional.length > 0) lines.push(`  Optional disaggregations: ${optional.join(", ")}`);
  if (avail.note) lines.push(`  [render note: ${avail.note}]`);
  const missingRequired = metric.requiredDisaggregationOptions.filter((r) =>
    !avail.options.includes(r)
  );
  if (missingRequired.length > 0) {
    lines.push(
      `  [render note: required option(s) not derivable from declared columns: ${missingRequired.join(", ")}]`,
    );
  }
  if (metric.vizPresets.length > 0) {
    lines.push(`  Visualization presets:`);
    for (const preset of metric.vizPresets) {
      const dateFormat = preset.config.d.timeseriesGrouping === "year" ? "YYYY" : "YYYYMM";
      const filterNote = preset.allowedFilters.length > 0
        ? ` — filters: ${preset.allowedFilters.join(", ")}`
        : "";
      const rep = replicateByProp(preset);
      const repNote = rep !== undefined ? ` ** REQUIRES selectedReplicant: one ${rep} value **` : "";
      lines.push(`    - ${preset.id}: ${preset.label.en} (${dateFormat})${filterNote}${repNote}`);
      if (preset.importantNotes) lines.push(`      NOTE: ${preset.importantNotes.en}`);
    }
  }
  lines.push("");
  return lines;
}

function dataHeader(def: ModuleDefinitionGithub, metric: MetricDefinitionGithub): string[] {
  const lines: string[] = [];
  lines.push("# METRIC DATA");
  lines.push("=".repeat(80));
  lines.push("");
  lines.push(`**Metric ID (metricId):** ${metric.id}`);
  lines.push(
    `**Metric Label:** ${metric.label.en}${metric.variantLabel ? ` [${metric.variantLabel.en}]` : ""}`,
  );
  if (metric.formatAs === "indicator") {
    lines.push(
      "**Format:** varies by indicator — each value is in its own indicator's format (percent values are 0-1 fractions; rate_per_10k values are written in the CSV already scaled to counts per 10,000; per-indicator formats are listed in the Dimension Summary below)",
    );
  } else {
    lines.push(`**Format:** ${metric.formatAs}`);
  }
  if (metric.valueProps.length > 0) {
    lines.push("");
    lines.push("**Value properties:**");
    for (const p of metric.valueProps) {
      lines.push(`  - ${p}: ${metric.valueLabelReplacements[p] || p}`);
    }
  }
  if (metric.importantNotes) {
    lines.push("");
    lines.push(`**IMPORTANT:** ${metric.importantNotes.en}`);
  }
  const ai = metric.aiDescription;
  if (ai) {
    lines.push("");
    lines.push(`**Methodology:** ${ai.methodology.en}`);
    lines.push(`**Interpretation:** ${ai.interpretation.en}`);
    lines.push(`**Typical range:** ${ai.typicalRange.en}`);
    if (ai.caveats) lines.push(`**Caveats:** ${ai.caveats.en}`);
    lines.push(`**Disaggregation guidance:** ${ai.disaggregationGuidance.en}`);
  }
  lines.push("");
  lines.push("(The live tool continues with the Dimension Summary, period coverage and CSV, which depend on package data and are not rendered here.)");
  lines.push("");
  const avail = availableDisaggregations(def, metric);
  lines.push(`[render note: disaggregations derivable from declared columns: ${avail.options.join(", ") || "(none)"}${avail.note ? `; ${avail.note}` : ""}]`);
  lines.push("");
  return lines;
}

function computationSummary(metric: MetricDefinitionGithub): string[] {
  const lines: string[] = [];
  lines.push("[ground truth the AI does NOT see]");
  lines.push(`  resultsObjectId: ${metric.resultsObjectId}`);
  lines.push(`  valueFunc: ${metric.valueFunc}; valueProps: ${metric.valueProps.join(", ")}`);
  if (metric.postAggregationExpression) {
    const ing = metric.postAggregationExpression.ingredientValues
      .map((i) => `${i.func}(${i.prop})`).join(", ");
    lines.push(`  postAggregationExpression: ${metric.postAggregationExpression.expression}   over ${ing}`);
  }
  if (metric.catalogExpressionEvaluation) {
    lines.push(`  catalogExpressionEvaluation over: ${metric.catalogExpressionEvaluation.ingredientProps.join(", ")}`);
  }
  lines.push(`  requiredDisaggregationOptions: ${metric.requiredDisaggregationOptions.join(", ") || "(none)"}`);
  lines.push(`  hide: ${metric.hide}`);
  return lines;
}

function translatedFieldsAppendix(metric: MetricDefinitionGithub): string[] {
  const lines: string[] = [];
  lines.push("[instance-language fields: on a fr or pt instance the AI sees these instead of the English]");
  const show = (name: string, v: Lang | null) => {
    if (!v) return;
    lines.push(`  ${name}.fr: ${v.fr}`);
    lines.push(`  ${name}.pt: ${v.pt ?? "(missing, falls back to en)"}`);
  };
  show("label", metric.label);
  show("variantLabel", metric.variantLabel);
  show("importantNotes", metric.importantNotes);
  return lines;
}

// Combined catalog, sorted by id across modules, as get_available_metrics does.
const allMetrics = modules.flatMap(({ def }) =>
  def.metrics.filter((m) => !m.hide).map((m) => ({ def, m }))
);
allMetrics.sort((a, b) => a.m.id.localeCompare(b.m.id));
const catalog: string[] = [
  "AVAILABLE METRICS",
  "=".repeat(80),
  "",
  "Query with get_metric_data for data and detailed context.",
  "Required disaggregations are auto-included.",
  "Period formats: period_id (YYYYMM), year (YYYY), month (1-12 for seasonal).",
  "",
];
for (const { def, m } of allMetrics) catalog.push(...catalogEntry(def, m));
await Deno.writeTextFile(`${outDir}/catalog.md`, catalog.join("\n"));

// Modules list and settings (SPA copilot only).
const modLines: string[] = ["AVAILABLE MODULES", "=".repeat(80), ""];
for (const { dir, def } of modules) {
  modLines.push(`ID: ${dir}`);
  modLines.push(`Name: ${def.label.en}`);
  modLines.push(`Has Parameters: ${def.configRequirements.parameters.length > 0}`);
  modLines.push(`Last Run: (runtime)`);
  modLines.push(`Metrics: ${def.metrics.length}`);
  modLines.push("-".repeat(80));
  modLines.push("");
}
modLines.push("");
for (const { dir, def } of modules) {
  modLines.push(`MODULE SETTINGS: ${dir}`);
  modLines.push("=".repeat(80));
  modLines.push("");
  modLines.push(`Name: ${def.label.en}`);
  modLines.push("");
  modLines.push("CURRENT SETTINGS (rendered with definition defaults)");
  modLines.push("-".repeat(80));
  modLines.push("");
  if (def.configRequirements.parameters.length === 0) {
    modLines.push("No parameters configured");
  } else {
    modLines.push("Parameters:");
    for (const p of def.configRequirements.parameters) {
      modLines.push(`  ${p.description.en}: ${p.input.defaultValue}`);
      modLines.push(`      [replacementString ${p.replacementString}, not shown to the AI]`);
    }
  }
  modLines.push("");
  modLines.push("=".repeat(80));
  modLines.push("");
}
await Deno.writeTextFile(`${outDir}/modules.md`, modLines.join("\n"));

// Per-module file: catalog entry, data header, ground truth, translations.
for (const { dir, def } of modules) {
  const lines: string[] = [];
  lines.push(`# ${dir}: ${def.label.en}`);
  lines.push("");
  lines.push(`Module label fr: ${def.label.fr}`);
  lines.push(`Module label pt: ${def.label.pt ?? "(missing, falls back to en)"}`);
  lines.push(`Prerequisites: ${def.prerequisites.join(", ") || "(none)"}`);
  lines.push("");
  const visible = def.metrics.filter((m) => !m.hide);
  const hidden = def.metrics.filter((m) => m.hide);
  for (const m of visible) {
    lines.push("-".repeat(80));
    lines.push(`## ${m.id}`);
    lines.push("");
    lines.push("### As listed by get_available_metrics");
    lines.push("```");
    lines.push(...catalogEntry(def, m));
    lines.push("```");
    lines.push("### As headed by get_metric_data");
    lines.push("```");
    lines.push(...dataHeader(def, m));
    lines.push("```");
    lines.push("```");
    lines.push(...computationSummary(m));
    lines.push(...translatedFieldsAppendix(m));
    lines.push("```");
    lines.push("");
  }
  if (hidden.length > 0) {
    lines.push("-".repeat(80));
    lines.push("## Hidden metrics (never listed to the AI; still queryable by id)");
    lines.push("");
    for (const m of hidden) {
      lines.push("```");
      lines.push(...catalogEntry(def, m));
      lines.push(...computationSummary(m));
      lines.push("```");
    }
  }
  await Deno.writeTextFile(`${outDir}/${dir}.md`, lines.join("\n"));
}

console.log(`Rendered ${modules.length} modules, ${allMetrics.length} visible metrics, to ${outDir}`);
