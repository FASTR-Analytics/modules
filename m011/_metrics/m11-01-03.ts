import {
  bandChartMetric,
  childrenOfParentPreset,
  singleAreaPreset,
} from "../_band_chart.ts";

export const metric = bandChartMetric({
  id: "m11-01-03",
  level: 3,
  vizPresets: [
    childrenOfParentPreset(3),
    singleAreaPreset(3),
  ],
});
