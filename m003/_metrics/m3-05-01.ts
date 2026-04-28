import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "disruption-chart-single-admin-area-3-multiple-admin-area-4",
    label: {
      en: "Disruptions and surpluses (single Admin Area 3, multiple Admin Area 4, multiple indicators)",
      fr: "Perturbations et excédents (unique Zone administrative 3, plusieurs Zones administratives 4, plusieurs indicateurs)",
    },
    description: {
      en: "Area chart showing actual vs expected service volume, for a single Admin Area 3, with multiple Admin Areas 4 and multiple indicators",
      fr: "Graphique en aires montrant le volume de services réel vs attendu, pour une unique Zone administrative 3, avec plusieurs Zones administratives 4 et plusieurs indicateurs",
    },
    allowedFilters: ["indicator_common_id", "admin_area_4"],
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
    createDefaultVisualizationOnInstall: "d4e1fad8-d7f4-46b3-8d56-7ca6b42a5cad",
  },
  {
    id: "disruption-chart-single-admin-area-4",
    label: {
      en: "Disruptions and surpluses (single Admin Area 4, multiple indicators)",
      fr: "Perturbations et excédents (unique Zone administrative 4, plusieurs indicateurs)",
    },
    description: {
      en: "Area chart showing actual vs expected service volume, for a single Admin Area 4 and multiple indicators",
      fr: "Graphique en aires montrant le volume de services réel vs attendu, pour une unique Zone administrative 4 et plusieurs indicateurs",
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
            disOpt: "admin_area_4",
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
    createDefaultVisualizationOnInstall: "8eb6e1d5-8bb7-48a1-a91a-9e6573ef8640",
  },
];

export const metric: MetricDefinitionGithub = {
  id: "m3-05-01",
  resultsObjectId: "M3_disruptions_analysis_admin_area_4.csv",
  label: {
    en: "Actual vs expected service volume",
    fr: "Volume de services réel vs attendu",
  },
  variantLabel: {
    en: "Admin area 4",
    fr: "Zone administrative 4",
  },
  valueProps: ["count_sum", "count_expected_if_above_diff_threshold"],
  valueFunc: "SUM",
  valueLabelReplacements: {
    count_sum: "Actual service volume",
    count_expected_if_above_diff_threshold: "Expected service volume",
  },
  requiredDisaggregationOptions: ["indicator_common_id", "admin_area_4"],
  formatAs: "number",
  aiDescription: {
    summary: {
      en: "Comparison of actual vs expected service volumes at admin area 4 (sub-district) level.",
      fr: "Comparaison des volumes réels vs attendus au niveau de la zone administrative 4 (sous-district).",
    },
    methodology: {
      en: "Sub-district level disruption analysis. Only generated when RUN_ADMIN_AREA_4_ANALYSIS parameter is enabled.",
      fr: "Analyse de perturbation au niveau sous-district. Généré uniquement si le paramètre RUN_ADMIN_AREA_4_ANALYSIS est activé.",
    },
    interpretation: {
      en: "Highest geographic resolution for disruption detection. Use for micro-level targeting and facility-level support.",
      fr: "Plus haute résolution géographique pour la détection de perturbations. Utiliser pour le ciblage micro-niveau et le soutien au niveau de l'établissement.",
    },
    typicalRange: {
      en: "Highly variable by sub-district. Expect greater volatility than higher geographic levels.",
      fr: "Très variable selon le sous-district. S'attendre à plus de volatilité que les niveaux géographiques supérieurs.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id and admin_area_4 (both required). Only available when RUN_ADMIN_AREA_4_ANALYSIS enabled.",
      fr: "Toujours désagréger par indicator_common_id et admin_area_4 (tous deux requis). Disponible uniquement si RUN_ADMIN_AREA_4_ANALYSIS activé.",
    },
    caveats: null,
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
