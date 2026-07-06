import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [];

export const metric: MetricDefinitionGithub = {
  id: "m3-02-03",
  resultsObjectId: "M3_disruptions_analysis_admin_area_1.csv",
  label: {
    en: "Actual vs expected service volume",
    fr: "Volume de services réel vs attendu",
    pt: "Volume de serviços real vs esperado",
  },
  variantLabel: {
    en: "National",
    fr: "National",
    pt: "Nacional",
  },
  valueProps: ["count_sum", "count_expect_sum"],
  valueFunc: "SUM",
  valueLabelReplacements: {
    count_sum: "Actual service volume",
    count_expect_sum: "Expected service volume",
  },
  requiredDisaggregationOptions: ["indicator_common_id"],
  formatAs: "number",
  aiDescription: {
    summary: {
      en: "Comparison of actual reported service volumes against model-predicted expected volumes at national level.",
      fr: "Comparaison des volumes de services réels déclarés par rapport aux volumes attendus prédits par modèle au niveau national.",
      pt: "Comparação dos volumes de serviços reais notificados com os volumes esperados previstos pelo modelo a nível nacional.",
    },
    methodology: {
      en: "SUM of actual reported counts vs SUM of expected counts from robust regression models. Expected volumes are calculated using panel regression with time trends, seasonal patterns, and disruption flags. The full expected series is reported for every period (no difference threshold is applied).",
      fr: "Somme des comptes réels déclarés vs somme des comptes attendus des modèles de régression robustes. Les volumes attendus sont calculés en utilisant la régression de panel avec des tendances temporelles, des modèles saisonniers et des indicateurs de perturbation. La série attendue complète est rapportée pour chaque période (aucun seuil de différence n'est appliqué).",
      pt: "SOMA das contagens reais notificadas vs SOMA das contagens esperadas a partir de modelos de regressão robustos. Os volumes esperados são calculados utilizando regressão de painel com tendências temporais, padrões sazonais e indicadores de perturbação. A série esperada completa é notificada para cada período (não é aplicado qualquer limiar de diferença).",
    },
    interpretation: {
      en: "Deviations between actual and expected volumes indicate potential service disruptions or anomalies. Periods where actual falls below expected suggest service delivery problems; actual above expected may indicate data quality issues, campaigns, or genuine increases in demand.",
      fr: "Les écarts entre volumes réels et attendus indiquent des perturbations potentielles de services ou des anomalies. Les périodes où le réel est inférieur à l'attendu suggèrent des problèmes de prestation de services; le réel supérieur à l'attendu peut indiquer des problèmes de qualité des données, des campagnes ou de véritables augmentations de la demande.",
      pt: "Os desvios entre os volumes reais e esperados indicam potenciais perturbações dos serviços ou anomalias. Os períodos em que o real fica abaixo do esperado sugerem problemas na prestação de serviços; o real acima do esperado pode indicar problemas de qualidade dos dados, campanhas ou aumentos genuínos da procura.",
    },
    typicalRange: {
      en: "Expected volumes should closely track actual volumes in stable periods. Deviations during crises, campaigns, or data system changes are normal.",
      fr: "Les volumes attendus devraient suivre de près les volumes réels pendant les périodes stables. Les écarts pendant les crises, les campagnes ou les changements du système de données sont normaux.",
      pt: "Os volumes esperados devem acompanhar de perto os volumes reais nos períodos estáveis. Os desvios durante crises, campanhas ou alterações do sistema de dados são normais.",
    },
    caveats: {
      en: "Model quality depends on sufficient historical data and stable baseline periods. Unlike the disruptions-and-surpluses chart metric, the expected value is reported for every period regardless of the size of the difference.",
      fr: "La qualité du modèle dépend de données historiques suffisantes et de périodes de référence stables. Contrairement à la métrique du graphique des perturbations et excédents, la valeur attendue est rapportée pour chaque période quelle que soit la taille de la différence.",
      pt: "A qualidade do modelo depende de dados históricos suficientes e de períodos de referência estáveis. Ao contrário da métrica do gráfico de perturbações e excedentes, o valor esperado é notificado para cada período independentemente da dimensão da diferença.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id (required) to see service-specific patterns. Time series visualization is essential for identifying disruption periods and recovery trends.",
      fr: "Toujours désagréger par indicator_common_id (requis) pour voir les modèles spécifiques au service. La visualisation en séries temporelles est essentielle pour identifier les périodes de perturbation et les tendances de récupération.",
      pt: "Desagregar sempre por indicator_common_id (obrigatório) para ver os padrões específicos de cada serviço. A visualização em séries temporais é essencial para identificar os períodos de perturbação e as tendências de recuperação.",
    },
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
