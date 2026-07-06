import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [];

export const metric: MetricDefinitionGithub = {
  id: "m6-03-02",
  resultsObjectId: "M6_coverage_estimation_admin3.csv",
  valueProps: ["coverage_cov"],
  valueFunc: "AVG",
  valueLabelReplacements: {
    coverage_cov: "Coverage calculated from HMIS data",
  },
  label: {
    en: "Coverage (HMIS only)",
    fr: "Couverture (HMIS uniquement)",
    pt: "Cobertura (apenas SIS)",
  },
  variantLabel: {
    en: "Admin area 3",
    fr: "Zone administrative 3",
    pt: "Área administrativa 3",
  },
  requiredDisaggregationOptions: [
    "indicator_common_id",
    "admin_area_3",
    "year",
  ],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "HMIS-only coverage estimates at district level for operational monitoring.",
      fr: "Estimations de couverture HMIS uniquement au niveau du district pour la surveillance opérationnelle.",
      pt: "Estimativas de cobertura apenas do SIS ao nível distrital para monitorização operacional.",
    },
    methodology: {
      en: "AVG of HMIS-derived coverage at district level. Excludes survey complexity for simplified operational use.",
      fr: "Moyenne de la couverture dérivée du HMIS au niveau du district. Exclut la complexité de l'enquête pour une utilisation opérationnelle simplifiée.",
      pt: "Média da cobertura derivada do SIS ao nível distrital. Exclui a complexidade dos inquéritos para uma utilização operacional simplificada.",
    },
    interpretation: {
      en: "District-level operational metric. Use for micro-level performance tracking and facility supervision planning.",
      fr: "Métrique opérationnelle au niveau du district. Utiliser pour le suivi de la performance micro-niveau.",
      pt: "Métrica operacional ao nível distrital. Utilizar para o acompanhamento do desempenho ao micro-nível e o planeamento da supervisão das unidades sanitárias.",
    },
    typicalRange: {
      en: "0-100%. Highest variation expected at this level; use for operational monitoring only.",
      fr: "0-100%. Variation la plus élevée attendue à ce niveau; utiliser uniquement pour la surveillance opérationnelle.",
      pt: "0-100%. Espera-se a maior variação a este nível; utilizar apenas para monitorização operacional.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id, admin_area_3, and year (all required). Simplified metric for district operations.",
      fr: "Toujours désagréger par indicator_common_id, admin_area_3 et year (tous requis). Métrique simplifiée pour les opérations de district.",
      pt: "Desagregar sempre por indicator_common_id, admin_area_3 e year (todos obrigatórios). Métrica simplificada para as operações distritais.",
    },
    caveats: null,
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
