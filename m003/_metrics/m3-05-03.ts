import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [];

export const metric: MetricDefinitionGithub = {
  id: "m3-05-03",
  resultsObjectId: "M3_disruptions_analysis_admin_area_4.csv",
  label: {
    en: "Actual vs expected service volume",
    fr: "Volume de services réel vs attendu",
  },
  variantLabel: {
    en: "Admin area 4",
    fr: "Zone administrative 4",
  },
  valueProps: ["count_sum", "count_expect_sum"],
  valueFunc: "SUM",
  valueLabelReplacements: {
    count_sum: "Actual service volume",
    count_expect_sum: "Expected service volume",
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
