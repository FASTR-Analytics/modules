import type { MetricDefinitionGithub, VizPreset } from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [];

export const metric: MetricDefinitionGithub = {
  id: "m6-02-02",
  resultsObjectId: "M6_coverage_estimation_admin2.csv",
  valueProps: ["coverage_cov"],
  valueFunc: "AVG",
  valueLabelReplacements: {
    coverage_cov: "Coverage calculated from HMIS data",
  },
  label: {
    en: "Coverage (HMIS only)",
    fr: "Couverture (HMIS uniquement)",
  },
  variantLabel: {
    en: "Admin area 2",
    fr: "Zone administrative 2",
  },
  requiredDisaggregationOptions: ["indicator_common_id", "admin_area_2", "year"],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "HMIS-only coverage estimates at admin area 2 level for simplified regional monitoring.",
      fr: "Estimations de couverture HMIS uniquement au niveau de la zone administrative 2 pour une surveillance régionale simplifiée.",
    },
    methodology: {
      en: "AVG of HMIS-derived coverage only (excludes survey and projected values). Useful for operational dashboards focusing solely on routine data.",
      fr: "Moyenne de la couverture dérivée du HMIS uniquement (exclut les valeurs d'enquête et projetées). Utile pour les tableaux de bord opérationnels.",
    },
    interpretation: {
      en: "Simplified metric for routine monitoring without survey complexity. Use for operational decision-making and trend tracking.",
      fr: "Métrique simplifiée pour la surveillance de routine sans complexité d'enquête. Utiliser pour la prise de décision opérationnelle.",
    },
    typicalRange: {
      en: "0-100%. Should align with full coverage metric (m6-02-01) but may differ from survey estimates.",
      fr: "0-100%. Devrait s'aligner avec la métrique de couverture complète (m6-02-01) mais peut différer des estimations d'enquête.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id, admin_area_2, and year (all required). Use for regional operational monitoring.",
      fr: "Toujours désagréger par indicator_common_id, admin_area_2 et year (tous requis). Utiliser pour la surveillance opérationnelle régionale.",
    },
    caveats: null,
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
