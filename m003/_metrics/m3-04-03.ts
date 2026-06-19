import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [];

export const metric: MetricDefinitionGithub = {
  id: "m3-04-03",
  resultsObjectId: "M3_disruptions_analysis_admin_area_3.csv",
  label: {
    en: "Actual vs expected service volume",
    fr: "Volume de services réel vs attendu",
  },
  variantLabel: {
    en: "Admin area 3",
    fr: "Zone administrative 3",
  },
  valueProps: ["count_sum", "count_expect_sum"],
  valueFunc: "SUM",
  valueLabelReplacements: {
    count_sum: "Actual service volume",
    count_expect_sum: "Expected service volume",
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
