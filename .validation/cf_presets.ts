import { z } from "zod";
import { cfStorageSchema } from "./conditional_formatting_standalone.ts";

type CfStorage = z.infer<typeof cfStorageSchema>;

const CF_LIGHTER_GREEN = "#68C690";
const CF_LIGHTER_YELLOW = "#F6D982";
const CF_LIGHTER_RED = "#F18989";
const CF_NO_DATA = "#ffffff";

// Base with all scale fields defaulted (required by schema)
const CF_BASE_SCALE_FIELDS: Pick<
  CfStorage,
  | "cfScalePaletteKind"
  | "cfScalePalettePreset"
  | "cfScaleCustomFrom"
  | "cfScaleCustomMid"
  | "cfScaleCustomTo"
  | "cfScaleReverse"
  | "cfScaleSteps"
  | "cfScaleDomainKind"
  | "cfScaleDomainMin"
  | "cfScaleDomainMax"
  | "cfScaleNoDataColor"
> = {
  cfScalePaletteKind: "preset",
  cfScalePalettePreset: "",
  cfScaleCustomFrom: "",
  cfScaleCustomMid: "",
  cfScaleCustomTo: "",
  cfScaleReverse: false,
  cfScaleSteps: 0,
  cfScaleDomainKind: "auto",
  cfScaleDomainMin: 0,
  cfScaleDomainMax: 1,
  cfScaleNoDataColor: "",
};

function thresholds(
  cutoffs: number[],
  buckets: Array<{ color: string }>,
  direction: "higher-is-better" | "lower-is-better" = "higher-is-better",
): CfStorage {
  return {
    cfMode: "thresholds",
    cfThresholdCutoffs: cutoffs,
    cfThresholdBuckets: buckets,
    cfThresholdDirection: direction,
    cfThresholdNoDataColor: CF_NO_DATA,
    ...CF_BASE_SCALE_FIELDS,
  };
}

// Higher is better: red < yellow < green
export const CF_90_80 = thresholds(
  [0.8, 0.9],
  [{ color: CF_LIGHTER_RED }, { color: CF_LIGHTER_YELLOW }, { color: CF_LIGHTER_GREEN }],
);

export const CF_80_70 = thresholds(
  [0.7, 0.8],
  [{ color: CF_LIGHTER_RED }, { color: CF_LIGHTER_YELLOW }, { color: CF_LIGHTER_GREEN }],
);

// Lower is better: green < yellow < red
export const CF_01_03 = thresholds(
  [0.01, 0.03],
  [{ color: CF_LIGHTER_GREEN }, { color: CF_LIGHTER_YELLOW }, { color: CF_LIGHTER_RED }],
  "lower-is-better",
);

export const CF_05_10 = thresholds(
  [0.05, 0.1],
  [{ color: CF_LIGHTER_GREEN }, { color: CF_LIGHTER_YELLOW }, { color: CF_LIGHTER_RED }],
  "lower-is-better",
);

export const CF_10_20 = thresholds(
  [0.1, 0.2],
  [{ color: CF_LIGHTER_GREEN }, { color: CF_LIGHTER_YELLOW }, { color: CF_LIGHTER_RED }],
  "lower-is-better",
);

// Diverging: red | neutral | green
export const CF_NEG10_POS10 = thresholds(
  [-0.1, 0.1],
  [{ color: CF_LIGHTER_RED }, { color: "#e0e0e0" }, { color: CF_LIGHTER_GREEN }],
);
