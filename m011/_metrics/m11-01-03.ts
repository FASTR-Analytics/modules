import {
  bandChartMetric,
  childrenOfParentPreset,
  singleAreaPreset,
} from "../_band_chart.ts";

export const metric = bandChartMetric({
  id: "m11-01-03",
  level: 3,
  vizPresets: [
    childrenOfParentPreset(3, "403b0757-4158-41e7-a9f2-e7c703ff9175"),
    singleAreaPreset(3, "adc9ad91-2449-4d50-bbea-7f4dbf61920b"),
  ],
});
