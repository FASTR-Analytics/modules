import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [];

export const metric: MetricDefinitionGithub = {
  id: "m3-03-03",
  resultsObjectId: "M3_disruptions_analysis_admin_area_2.csv",
  label: {
    en: "Actual vs expected service volume",
    fr: "Volume de services réel vs attendu",
    pt: "Volume de serviços real vs esperado",
  },
  variantLabel: {
    en: "Admin area 2",
    fr: "Zone administrative 2",
    pt: "Área administrativa 2",
  },
  valueProps: ["count_sum", "count_expect_sum"],
  valueFunc: "SUM",
  valueLabelReplacements: {
    count_sum: "Actual service volume",
    count_expect_sum: "Expected service volume",
  },
  requiredDisaggregationOptions: ["indicator_common_id", "admin_area_2"],
  formatAs: "number",
  aiDescription: {
    summary: {
      en: "Comparison of actual reported service volumes against model-predicted expected volumes at admin area 2 (province/state) level.",
      fr: "Comparaison des volumes de services réels contre volumes attendus au niveau de la zone administrative 2 (province/état).",
      pt: "Comparação dos volumes de serviços reais notificados com os volumes esperados previstos pelo modelo ao nível da área administrativa 2 (província/estado).",
    },
    methodology: {
      en: "SUM of actual vs expected counts from area-specific robust regression models. Expected volumes account for local time trends and seasonal patterns. The full expected series is reported for every period (no difference threshold is applied).",
      fr: "Somme des comptes réels vs attendus des modèles de régression robustes spécifiques à la zone. Les volumes attendus tiennent compte des tendances temporelles locales et des modèles saisonniers. La série attendue complète est rapportée pour chaque période (aucun seuil de différence n'est appliqué).",
      pt: "SOMA das contagens reais vs esperadas a partir de modelos de regressão robustos específicos da área. Os volumes esperados têm em conta as tendências temporais locais e os padrões sazonais. A série esperada completa é notificada para cada período (não é aplicado qualquer limiar de diferença).",
    },
    interpretation: {
      en: "Enables subnational identification of service disruptions. Compare across admin areas to identify geographic hotspots of service delivery problems. Areas with persistent negative deviations need targeted support.",
      fr: "Permet l'identification sous-nationale des perturbations de services. Comparer entre les zones administratives pour identifier les points chauds géographiques de problèmes de prestation de services. Les zones avec des écarts négatifs persistants nécessitent un soutien ciblé.",
      pt: "Permite a identificação subnacional das perturbações dos serviços. Comparar entre áreas administrativas para identificar os focos geográficos de problemas na prestação de serviços. As áreas com desvios negativos persistentes necessitam de apoio direcionado.",
    },
    typicalRange: {
      en: "Expected to closely match actual in stable regions. Deviations indicate local disruptions, data issues, or demand changes.",
      fr: "Attendu pour correspondre étroitement au réel dans les régions stables. Les écarts indiquent des perturbations locales, des problèmes de données ou des changements de la demande.",
      pt: "Espera-se que corresponda de perto ao real nas regiões estáveis. Os desvios indicam perturbações locais, problemas de dados ou alterações da procura.",
    },
    caveats: {
      en: "Smaller geographic areas may have more volatile patterns, making model predictions less reliable. Consider aggregating to higher levels if area-specific models perform poorly.",
      fr: "Les zones géographiques plus petites peuvent avoir des modèles plus volatils, rendant les prédictions du modèle moins fiables. Considérer l'agrégation à des niveaux supérieurs si les modèles spécifiques aux zones fonctionnent mal.",
      pt: "As áreas geográficas mais pequenas podem ter padrões mais voláteis, tornando as previsões do modelo menos fiáveis. Considerar a agregação a níveis superiores se os modelos específicos das áreas tiverem um desempenho fraco.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id and admin_area_2 (both required). Time series reveals when and where disruptions occurred. Map visualization effectively shows geographic distribution of service gaps.",
      fr: "Toujours désagréger par indicator_common_id et admin_area_2 (tous deux requis). Les séries temporelles révèlent quand et où les perturbations se sont produites. La visualisation cartographique montre efficacement la distribution géographique des écarts de services.",
      pt: "Desagregar sempre por indicator_common_id e admin_area_2 (ambos obrigatórios). As séries temporais revelam quando e onde ocorreram as perturbações. A visualização cartográfica mostra eficazmente a distribuição geográfica das lacunas de serviços.",
    },
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
