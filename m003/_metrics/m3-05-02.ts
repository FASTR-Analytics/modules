import type { MetricDefinitionGithub, VizPreset } from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "disruption-differences-table-single-admin-area-3-multiple-admin-area-4",
    label: {
      en: "Difference between actual and expected service volume (single Admin Area 3, multiple Admin Area 4, multiple indicators)",
      fr: "Différence entre le volume de services réel et attendu (unique Zone administrative 3, plusieurs Zones administratives 4, plusieurs indicateurs)",
    },
    description: {
      en: "Table showing percentage differences between actual vs expected service volume, with conditional formatting, for a single Admin Area 3, with multiple Admin Areas 4 and multiple indicators",
      fr: "Tableau montrant les différences en pourcentage entre le volume de services réel et attendu, avec mise en forme conditionnelle, pour une unique Zone administrative 3, avec plusieurs Zones administratives 4 et plusieurs indicateurs",
    },
    allowedFilters: ["indicator_common_id", "admin_area_4"],
    config: {
      d: {
        type: "table",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "col",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "col",
          },
          {
            disOpt: "admin_area_3",
            disDisplayOpt: "replicant",
          },
          {
            disOpt: "admin_area_4",
            disDisplayOpt: "row",
          },
        ],
        filterBy: [],
      },
      s: {
        cfMode: "thresholds",
        cfThresholdCutoffs: [-0.1, 0.1],
        cfThresholdBuckets: [
          {
            color: "#F18989",
          },
          {
            color: "#e0e0e0",
          },
          {
            color: "#68C690",
          },
        ],
        cfThresholdDirection: "higher-is-better",
        cfThresholdNoDataColor: "#ffffff",
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
        decimalPlaces: 0,
      },
      t: {
        caption: null,
        captionRelFontSize: null,
        subCaption: null,
        subCaptionRelFontSize: null,
        footnote: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
    createDefaultVisualizationOnInstall: "45735c34-5dd4-43ee-baad-346c74751d4d",
  },
];

export const metric: MetricDefinitionGithub = {
  id: "m3-05-02",
  resultsObjectId: "M3_disruptions_analysis_admin_area_4.csv",
  label: {
    en: "Difference between actual and expected service volume",
    fr: "Différence entre le volume de services réel et attendu",
  },
  variantLabel: {
    en: "Admin area 4",
    fr: "Zone administrative 4",
  },
  valueProps: ["pct_diff"],
  valueFunc: "identity",
  valueLabelReplacements: {
    pct_diff: "Percent difference",
  },
  postAggregationExpression: {
    ingredientValues: [
      {
        prop: "count_sum",
        func: "SUM",
      },
      {
        prop: "count_expect_sum",
        func: "SUM",
      },
    ],
    expression: "pct_diff = (count_sum - count_expect_sum)/count_expect_sum",
  },
  requiredDisaggregationOptions: ["indicator_common_id", "admin_area_4"],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Percentage difference between actual and expected service volumes at sub-district level.",
      fr: "Différence en pourcentage entre volumes réels et attendus au niveau sous-district.",
    },
    methodology: {
      en: "(actual - expected) / expected at sub-district level. Finest geographic granularity for performance assessment.",
      fr: "(réel - attendu) / attendu au niveau sous-district. Granularité géographique la plus fine pour l'évaluation.",
    },
    interpretation: {
      en: "Enables targeted facility-level interventions. Small sample sizes at this level may increase volatility.",
      fr: "Permet des interventions ciblées au niveau de l'établissement. Les petites tailles d'échantillon peuvent augmenter la volatilité.",
    },
    typicalRange: {
      en: "±20-50% variation common; high volatility expected at this granular level.",
      fr: "Variation de ±20-50% commune; forte volatilité attendue à ce niveau granulaire.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id and admin_area_4 (both required). Interpret with caution due to small sample sizes.",
      fr: "Toujours désagréger par indicator_common_id et admin_area_4 (tous deux requis). Interpréter avec prudence en raison de petites tailles d'échantillon.",
    },
    caveats: null,
  },
  importantNotes: null,
  hide: false,
  vizPresets,
};
