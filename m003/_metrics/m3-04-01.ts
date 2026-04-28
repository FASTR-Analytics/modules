import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "disruption-chart-single-admin-area-2-multiple-admin-area-3",
    label: {
      en: "Disruptions and surpluses (single Admin Area 2, multiple Admin Area 3, multiple indicators)",
      fr: "Perturbations et excédents (unique Zone administrative 2, plusieurs Zones administratives 3, plusieurs indicateurs)",
    },
    description: {
      en: "Area chart showing actual vs expected service volume, for a single Admin Area 2, with multiple Admin Areas 3 and multiple indicators",
      fr: "Graphique en aires montrant le volume de services réel vs attendu, pour une unique Zone administrative 2, avec plusieurs Zones administratives 3 et plusieurs indicateurs",
    },
    allowedFilters: ["indicator_common_id", "admin_area_3"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "series",
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
        filterBy: [],
      },
      s: {
        specialDisruptionsChart: true,
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
    createDefaultVisualizationOnInstall: "68e1b077-eb96-427e-9254-ea6b4a871b5d",
  },
  {
    id: "disruption-chart-single-admin-area-3",
    label: {
      en: "Disruptions and surpluses (single Admin Area 3, multiple indicators)",
      fr: "Perturbations et excédents (unique Zone administrative 3, plusieurs indicateurs)",
    },
    description: {
      en: "Area chart showing actual vs expected service volume, for a single Admin Area 3 and multiple indicators",
      fr: "Graphique en aires montrant le volume de services réel vs attendu, pour une unique Zone administrative 3 et plusieurs indicateurs",
    },
    allowedFilters: ["indicator_common_id"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "series",
        selectedReplicantValue: "anc1",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "cell",
          },
          {
            disOpt: "admin_area_3",
            disDisplayOpt: "replicant",
          },
        ],
        filterBy: [],
      },
      s: {
        specialDisruptionsChart: true,
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
    createDefaultVisualizationOnInstall: "1870d8fa-b7ff-453f-8974-184282643bee",
  },
];

export const metric: MetricDefinitionGithub = {
  id: "m3-04-01",
  resultsObjectId: "M3_disruptions_analysis_admin_area_3.csv",
  label: {
    en: "Actual vs expected service volume",
    fr: "Volume de services réel vs attendu",
  },
  variantLabel: {
    en: "Admin area 3",
    fr: "Zone administrative 3",
  },
  valueProps: ["count_sum", "count_expected_if_above_diff_threshold"],
  valueFunc: "SUM",
  valueLabelReplacements: {
    count_sum: "Actual service volume",
    count_expected_if_above_diff_threshold: "Expected service volume",
  },
  requiredDisaggregationOptions: ["indicator_common_id", "admin_area_3"],
  formatAs: "number",
  aiDescription: {
    summary: {
      en: "Comparison of actual vs expected service volumes at admin area 3 (district) level.",
      fr: "Comparaison des volumes réels vs attendus au niveau de la zone administrative 3 (district).",
    },
    methodology: {
      en: "District-level disruption analysis using robust regression models. Enables fine-grained geographic targeting.",
      fr: "Analyse de perturbation au niveau du district utilisant des modèles de régression robustes. Permet un ciblage géographique précis.",
    },
    interpretation: {
      en: "Identifies district-specific service delivery problems. Use for operational planning and targeted supervision.",
      fr: "Identifie les problèmes de prestation de services spécifiques au district. Utiliser pour la planification opérationnelle et la supervision ciblée.",
    },
    typicalRange: {
      en: "Varies by district size. Expect similar patterns to admin area 2 but with more volatility.",
      fr: "Varie selon la taille du district. S'attendre à des modèles similaires à la zone administrative 2 mais avec plus de volatilité.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id and admin_area_3 (both required). Time series and maps reveal district-level patterns.",
      fr: "Toujours désagréger par indicator_common_id et admin_area_3 (tous deux requis). Les séries temporelles et cartes révèlent les modèles au niveau du district.",
    },
    caveats: null,
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
