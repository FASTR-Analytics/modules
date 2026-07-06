import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "disruption-chart",
    label: {
      en: "Disruptions and surpluses (national)",
      fr: "Perturbations et excédents (national)",
      pt: "Perturbações e excedentes (nacional)",
    },
    description: {
      en: "Area chart showing actual vs expected service volume nationally",
      fr: "Graphique en aires montrant le volume de services réel vs attendu au niveau national",
      pt: "Gráfico de áreas que mostra o volume de serviços real vs esperado a nível nacional",
    },
    createDefaultVisualizationOnInstall: "e51a15fd-acfc-4da9-8797-b462b9626cff",
    allowedFilters: ["indicator_common_id"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "cell",
          },
        ],
        filterBy: [],
      },
      s: {
        specialDisruptionsChart: true,
      },
      t: {
        caption: {
          en: "Disruptions and surpluses in service volume, nationally",
          fr: "Perturbations et excédents du volume de services, au niveau national",
          pt: "Perturbações e excedentes no volume de serviços, a nível nacional",
        },
        subCaption: {
          en: "DATE_RANGE",
          fr: "PLAGE_DE_DATES",
          pt: "INTERVALO_DE_DATAS",
        },
        footnote: {
          en: "This graph quantifies changes in service volume compared to historical trends and accounting for seasonality. These signals should be triangulated to other data and contextual knowledge to determine if the results are an artifact of data quality. Unexpected volume changes are estimated by comparing the observed volume to the expected volume based on historical trends and seasonality. Previous large unexpected changes in the historical data are removed. This analysis is an interrupted time series regression with facility-level fixed effects.",
          fr: "Ce graphique quantifie les changements du volume de services par rapport aux tendances historiques et en tenant compte de la saisonnalité. Ces signaux doivent être triangulés avec d'autres données et connaissances contextuelles pour déterminer si les résultats sont un artefact de la qualité des données. Les changements de volume inattendus sont estimés en comparant le volume observé au volume attendu basé sur les tendances historiques et la saisonnalité. Les grands changements inattendus précédents dans les données historiques sont supprimés. Cette analyse est une régression de séries temporelles interrompues avec effets fixes au niveau de l'établissement.",
          pt: "Este gráfico quantifica as alterações no volume de serviços em comparação com as tendências históricas e tendo em conta a sazonalidade. Estes sinais devem ser triangulados com outros dados e conhecimento contextual para determinar se os resultados são um artefacto da qualidade dos dados. As alterações de volume inesperadas são estimadas comparando o volume observado com o volume esperado com base nas tendências históricas e na sazonalidade. As grandes alterações inesperadas anteriores nos dados históricos são removidas. Esta análise é uma regressão de séries temporais interrompidas com efeitos fixos ao nível do estabelecimento.",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
];

export const metric: MetricDefinitionGithub = {
  id: "m3-02-01",
  resultsObjectId: "M3_disruptions_analysis_admin_area_1.csv",
  label: {
    en: "Disruptions and surpluses",
    fr: "Perturbations et excédents",
    pt: "Perturbações e excedentes",
  },
  variantLabel: {
    en: "National",
    fr: "National",
    pt: "Nacional",
  },
  valueProps: ["count_sum", "count_expected_if_above_diff_threshold"],
  valueFunc: "SUM",
  valueLabelReplacements: {
    count_sum: "Actual service volume",
    count_expected_if_above_diff_threshold: "Expected service volume",
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
      en: "SUM of actual reported counts vs SUM of expected counts from robust regression models. Expected volumes are calculated using panel regression with time trends, seasonal patterns, and disruption flags. Only periods with differences exceeding the configured threshold are included.",
      fr: "Somme des comptes réels déclarés vs somme des comptes attendus des modèles de régression robustes. Les volumes attendus sont calculés en utilisant la régression de panel avec des tendances temporelles, des modèles saisonniers et des indicateurs de perturbation. Seules les périodes avec des différences dépassant le seuil configuré sont incluses.",
      pt: "SOMA das contagens reais notificadas vs SOMA das contagens esperadas a partir de modelos de regressão robustos. Os volumes esperados são calculados utilizando regressão de painel com tendências temporais, padrões sazonais e indicadores de perturbação. São incluídos apenas os períodos com diferenças que excedem o limiar configurado.",
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
      en: "Model quality depends on sufficient historical data and stable baseline periods. Expected values are only shown when difference exceeds the configured DIFFPERCENT threshold (default 10%), so small deviations are filtered out.",
      fr: "La qualité du modèle dépend de données historiques suffisantes et de périodes de référence stables. Les valeurs attendues ne sont affichées que lorsque la différence dépasse le seuil DIFFPERCENT configuré (10% par défaut), donc les petits écarts sont filtrés.",
      pt: "A qualidade do modelo depende de dados históricos suficientes e de períodos de referência estáveis. Os valores esperados só são apresentados quando a diferença excede o limiar DIFFPERCENT configurado (10% por predefinição), pelo que os pequenos desvios são filtrados.",
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
