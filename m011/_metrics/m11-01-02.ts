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
    gridPreset(2, "df19e07e-88dd-41ea-aa40-29b4d82fb0dc"),
    singleAreaPreset(2, "c43c3d27-c1c0-4945-a4fa-21d93c5cf45a"),
    singleIndicatorPreset(2, "19283b6d-4a37-4843-9667-9d66979a5990"),
  ],
});
