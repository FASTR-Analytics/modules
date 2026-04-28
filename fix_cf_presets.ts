const CF_PRESETS: Record<string, { cutoffs: number[]; direction: string }> = {
  CF_90_80: { cutoffs: [0.8, 0.9], direction: "higher-is-better" },
  CF_80_70: { cutoffs: [0.7, 0.8], direction: "higher-is-better" },
  CF_01_03: { cutoffs: [0.01, 0.03], direction: "lower-is-better" },
  CF_05_10: { cutoffs: [0.05, 0.1], direction: "lower-is-better" },
  CF_10_20: { cutoffs: [0.1, 0.2], direction: "lower-is-better" },
  CF_NEG10_POS10: { cutoffs: [-0.1, 0.1], direction: "higher-is-better" },
};

function detectCfPreset(cutoffs: number[], direction: string): string | null {
  for (const [name, preset] of Object.entries(CF_PRESETS)) {
    if (
      preset.cutoffs.length === cutoffs.length &&
      preset.cutoffs.every((v, i) => v === cutoffs[i]) &&
      preset.direction === direction
    ) {
      return name;
    }
  }
  return null;
}

async function fixFile(path: string): Promise<boolean> {
  const content = await Deno.readTextFile(path);

  // Pattern to match the inlined cf fields
  const pattern = /cfMode: "thresholds",\s*cfThresholdCutoffs: \[([\d., -]+)\],\s*cfThresholdBuckets: \[\s*\{\s*color: "[^"]+",?\s*\},\s*\{\s*color: "[^"]+",?\s*\},\s*\{\s*color: "[^"]+",?\s*\},?\s*\],\s*cfThresholdDirection: "([^"]+)",\s*cfThresholdNoDataColor: "[^"]+",\s*cfScalePaletteKind: "[^"]*",\s*cfScalePalettePreset: "[^"]*",\s*cfScaleCustomFrom: "[^"]*",\s*cfScaleCustomMid: "[^"]*",\s*cfScaleCustomTo: "[^"]*",\s*cfScaleReverse: \w+,\s*cfScaleSteps: \d+,\s*cfScaleDomainKind: "[^"]*",\s*cfScaleDomainMin: [\d.]+,\s*cfScaleDomainMax: [\d.]+,\s*cfScaleNoDataColor: "[^"]*",?/g;

  let modified = content;
  let changed = false;

  modified = content.replace(pattern, (match, cutoffsStr, direction) => {
    const cutoffs = cutoffsStr.split(",").map((s: string) => parseFloat(s.trim()));
    const preset = detectCfPreset(cutoffs, direction);

    if (preset) {
      changed = true;
      return `...${preset},`;
    }
    return match;
  });

  if (changed) {
    // Add import if not present
    if (!modified.includes("cf_presets.ts")) {
      modified = modified.replace(
        'import type { MetricDefinitionGithub, VizPreset } from "../../.validation/_module_definition_github.ts";',
        'import type { MetricDefinitionGithub, VizPreset } from "../../.validation/_module_definition_github.ts";\nimport { CF_01_03, CF_05_10, CF_10_20, CF_80_70, CF_90_80, CF_NEG10_POS10 } from "../../.validation/cf_presets.ts";'
      );
    }
    await Deno.writeTextFile(path, modified);
    return true;
  }
  return false;
}

// Process all metric files
for await (const dirEntry of Deno.readDir(".")) {
  if (!dirEntry.isDirectory || !/^m\d+$/.test(dirEntry.name)) continue;

  const metricsDir = `./${dirEntry.name}/_metrics`;
  try {
    for await (const file of Deno.readDir(metricsDir)) {
      if (!file.name.endsWith(".ts")) continue;
      const path = `${metricsDir}/${file.name}`;
      const fixed = await fixFile(path);
      if (fixed) {
        console.log(`fixed ${path}`);
      }
    }
  } catch {
    // No _metrics folder
  }
}

console.log("\nDone.");
