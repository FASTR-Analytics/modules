import { bandChartMetric, byIndicatorPreset } from "../_band_chart.ts";

export const metric = bandChartMetric({
  id: "m11-01-01",
  level: 1,
  vizPresets: [byIndicatorPreset()],
});
