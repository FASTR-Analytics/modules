import {
  bandChartMetric,
  gridPreset,
  singleAreaPreset,
  singleIndicatorPreset,
} from "../_band_chart.ts";

export const metric = bandChartMetric({
  id: "m11-01-02",
  level: 2,
  vizPresets: [
    gridPreset(2),
    singleAreaPreset(2),
    singleIndicatorPreset(2),
  ],
});
