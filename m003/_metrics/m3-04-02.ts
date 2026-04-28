import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { CF_NEG10_POS10 } from "../../.validation/cf_presets.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "disruption-differences-table-single-admin-area-2-multiple-admin-area-3",
    label: {
      en: "Difference between actual and expected service volume (single Admin Area 2, multiple Admin Area 3, multiple indicators)",
      fr: "Différence entre le volume de services réel et attendu (unique Zone administrative 2, plusieurs Zones administratives 3, plusieurs indicateurs)",
    },
    description: {
      en: "Table showing percentage differences between actual vs expected service volume, with conditional formatting, for a single Admin Area 2, with multiple Admin Areas 3 and multiple indicators",
      fr: "Tableau montrant les différences en pourcentage entre le volume de services réel et attendu, avec mise en forme conditionnelle, pour une unique Zone administrative 2, avec plusieurs Zones administratives 3 et plusieurs indicateurs",
    },
    allowedFilters: ["indicator_common_id", "admin_area_3"],
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
            disOpt: "admin_area_2",
            disDisplayOpt: "replicant",
          },
          {
            disOpt: "admin_area_3",
            disDisplayOpt: "row",
          },
        ],
        selectedReplicantValue: "",
        filterBy: [],
      },
      s: {
        ...CF_NEG10_POS10,
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
    createDefaultVisualizationOnInstall: "82a91435-b090-4bfe-8477-fe864095155e",
  },
  {
    id: "disruption-differences-map",
    label: {
      en: "Disruption map by admin area",
      fr: "Carte des perturbations par zone administrative",
    },
    description: {
      en: "Map showing disruption levels by admin area for a single indicator and time period",
      fr: "Carte montrant les niveaux de perturbation par zone administrative pour un indicateur et une période unique",
    },
    allowedFilters: [],
    config: {
      d: {
        type: "map",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "cell",
        disaggregateBy: [
          {
            disOpt: "admin_area_3",
            disDisplayOpt: "mapArea",
          },
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "replicant",
          },
        ],
        selectedReplicantValue: "anc1",
        filterBy: [],
      },
      s: {
        cfMode: "scale",
        cfScalePaletteKind: "custom",
        cfScaleCustomFrom: "#fee0d2",
        cfScaleCustomTo: "#de2d26",
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
    createDefaultVisualizationOnInstall: "732d4f78-e53c-4c44-ac9f-234f616baa2c",
  },
];

export const metric: MetricDefinitionGithub = {
  id: "m3-04-02",
  resultsObjectId: "M3_disruptions_analysis_admin_area_3.csv",
  label: {
    en: "Difference between actual and expected service volume",
    fr: "Différence entre le volume de services réel et attendu",
  },
  variantLabel: {
    en: "Admin area 3",
    fr: "Zone administrative 3",
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
  requiredDisaggregationOptions: ["indicator_common_id", "admin_area_3"],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Percentage difference between actual and expected service volumes at district level.",
      fr: "Différence en pourcentage entre volumes réels et attendus au niveau du district.",
    },
    methodology: {
      en: "(actual - expected) / expected at district level. Quantifies local service delivery performance.",
      fr: "(réel - attendu) / attendu au niveau du district. Quantifie la performance locale de prestation.",
    },
    interpretation: {
      en: "Enables district-level performance monitoring. Prioritize districts with largest negative deviations.",
      fr: "Permet la surveillance de la performance au niveau du district. Prioriser les districts avec les plus grands écarts négatifs.",
    },
    typicalRange: {
      en: "±10-40% variation common at district level; >40% deviation warrants investigation.",
      fr: "Variation de ±10-40% commune au niveau du district; écart >40% nécessite investigation.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id and admin_area_3 (both required). Compare with admin area 2 for context.",
      fr: "Toujours désagréger par indicator_common_id et admin_area_3 (tous deux requis). Comparer avec la zone administrative 2 pour le contexte.",
    },
    caveats: null,
  },
  importantNotes: null,
  hide: false,
  vizPresets,
};
