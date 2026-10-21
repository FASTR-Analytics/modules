import {
  bandChartMetric,
  childrenOfParentPreset,
  singleAreaPreset,
} from "../_band_chart.ts";

export const metric = bandChartMetric({
  id: "m11-01-04",
  level: 4,
  vizPresets: [
    childrenOfParentPreset(4, "ab26c051-62e8-40fc-a562-2fca058375a3"),
    singleAreaPreset(4, "9040ce32-5daa-4886-814c-c6a1120f3e89"),
  ],
});
