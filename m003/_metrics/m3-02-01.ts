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
    },
    description: {
      en: "Area chart showing actual vs expected service volume nationally",
      fr: "Graphique en aires montrant le volume de services réel vs attendu au niveau national",
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
        },
        subCaption: {
          en: "DATE_RANGE",
          fr: "DATE_RANGE",
        },
        footnote: {
          en: "This graph quantifies changes in service volume compared to historical trends and accounting for seasonality. These signals should be triangulated to other data and contextual knowledge to determine if the results are an artifact of data quality. Unexpected volume changes are estimated by comparing the observed volume to the expected volume based on historical trends and seasonality. Previous large unexpected changes in the historical data are removed. This analysis is an interrupted time series regression with facility-level fixed effects.",
          fr: "Ce graphique quantifie les changements du volume de services par rapport aux tendances historiques et en tenant compte de la saisonnalité. Ces signaux doivent être triangulés avec d'autres données et connaissances contextuelles pour déterminer si les résultats sont un artefact de la qualité des données. Les changements de volume inattendus sont estimés en comparant le volume observé au volume attendu basé sur les tendances historiques et la saisonnalité. Les grands changements inattendus précédents dans les données historiques sont supprimés. Cette analyse est une régression de séries temporelles interrompues avec effets fixes au niveau de l'établissement.",
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
    en: "Actual vs expected service volume",
    fr: "Volume de services réel vs attendu",
  },
  variantLabel: {
    en: "National",
    fr: "National",
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
    },
    methodology: {
      en: "SUM of actual reported counts vs SUM of expected counts from robust regression models. Expected volumes are calculated using panel regression with time trends, seasonal patterns, and disruption flags. Only periods with differences exceeding the configured threshold are included.",
      fr: "Somme des comptes réels déclarés vs somme des comptes attendus des modèles de régression robustes. Les volumes attendus sont calculés en utilisant la régression de panel avec des tendances temporelles, des modèles saisonniers et des indicateurs de perturbation. Seules les périodes avec des différences dépassant le seuil configuré sont incluses.",
    },
    interpretation: {
      en: "Deviations between actual and expected volumes indicate potential service disruptions or anomalies. Periods where actual falls below expected suggest service delivery problems; actual above expected may indicate data quality issues, campaigns, or genuine increases in demand.",
      fr: "Les écarts entre volumes réels et attendus indiquent des perturbations potentielles de services ou des anomalies. Les périodes où le réel est inférieur à l'attendu suggèrent des problèmes de prestation de services; le réel supérieur à l'attendu peut indiquer des problèmes de qualité des données, des campagnes ou de véritables augmentations de la demande.",
    },
    typicalRange: {
      en: "Expected volumes should closely track actual volumes in stable periods. Deviations during crises, campaigns, or data system changes are normal.",
      fr: "Les volumes attendus devraient suivre de près les volumes réels pendant les périodes stables. Les écarts pendant les crises, les campagnes ou les changements du système de données sont normaux.",
    },
    caveats: {
      en: "Model quality depends on sufficient historical data and stable baseline periods. Expected values are only shown when difference exceeds the configured DIFFPERCENT threshold (default 10%), so small deviations are filtered out.",
      fr: "La qualité du modèle dépend de données historiques suffisantes et de périodes de référence stables. Les valeurs attendues ne sont affichées que lorsque la différence dépasse le seuil DIFFPERCENT configuré (10% par défaut), donc les petits écarts sont filtrés.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id (required) to see service-specific patterns. Time series visualization is essential for identifying disruption periods and recovery trends.",
      fr: "Toujours désagréger par indicator_common_id (requis) pour voir les modèles spécifiques au service. La visualisation en séries temporelles est essentielle pour identifier les périodes de perturbation et les tendances de récupération.",
    },
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
