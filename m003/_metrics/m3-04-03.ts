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
    pt: "Volume de serviços real vs esperado",
  },
  variantLabel: {
    en: "Admin area 3",
    fr: "Zone administrative 3",
    pt: "Área administrativa 3",
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
      pt: "Comparação dos volumes de serviços reais vs esperados ao nível da área administrativa 3 (distrito).",
    },
    methodology: {
      en: "District-level disruption analysis using robust regression models. Enables fine-grained geographic targeting.",
      fr: "Analyse de perturbation au niveau du district utilisant des modèles de régression robustes. Permet un ciblage géographique précis.",
      pt: "Análise de perturbação ao nível do distrito utilizando modelos de regressão robustos. Permite um direcionamento geográfico preciso.",
    },
    interpretation: {
      en: "Identifies district-specific service delivery problems. Use for operational planning and targeted supervision.",
      fr: "Identifie les problèmes de prestation de services spécifiques au district. Utiliser pour la planification opérationnelle et la supervision ciblée.",
      pt: "Identifica problemas de prestação de serviços específicos do distrito. Utilizar para o planeamento operacional e a supervisão direcionada.",
    },
    typicalRange: {
      en: "Varies by district size. Expect similar patterns to admin area 2 but with more volatility.",
      fr: "Varie selon la taille du district. S'attendre à des modèles similaires à la zone administrative 2 mais avec plus de volatilité.",
      pt: "Varia consoante a dimensão do distrito. Esperar padrões semelhantes aos da área administrativa 2, mas com maior volatilidade.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id and admin_area_3 (both required). Time series and maps reveal district-level patterns.",
      fr: "Toujours désagréger par indicator_common_id et admin_area_3 (tous deux requis). Les séries temporelles et cartes révèlent les modèles au niveau du district.",
      pt: "Desagregar sempre por indicator_common_id e admin_area_3 (ambos obrigatórios). As séries temporais e os mapas revelam os padrões ao nível do distrito.",
    },
    caveats: null,
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
