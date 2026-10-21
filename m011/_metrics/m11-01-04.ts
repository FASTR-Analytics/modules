import {
  bandChartMetric,
  childrenOfParentPreset,
  singleAreaPreset,
} from "../_band_chart.ts";

export const metric = bandChartMetric({
  id: "m11-01-04",
  level: 4,
  vizPresets: [
    childrenOfParentPreset(4),
    singleAreaPreset(4),
  ],
});
