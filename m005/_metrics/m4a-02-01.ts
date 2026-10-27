import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { DENOMINATOR_VALUE_LABELS } from "../_denominator_labels.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "coverage-timeseries",
    label: {
      en: "Coverage by denominator type",
      fr: "Couverture par type de dénominateur",
      pt: "Cobertura por tipo de denominador",
    },
    description: {
      en: "Timeseries comparing coverage across denominator sources",
      fr:
        "Séries temporelles comparant la couverture entre sources de dénominateur",
      pt:
        "Séries temporais que comparam a cobertura entre fontes de denominador",
    },
    allowedFilters: ["denominator_best_or_survey"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [
          {
            disOpt: "denominator_best_or_survey",
            disDisplayOpt: "series",
          },
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "replicant",
          },
        ],
        filterBy: [],
        selectedReplicantValue: "anc1",
      },
      s: {
        content: "lines",
        showDataLabels: true,
      },
      t: {
        caption: {
          en: "Coverage based on different denominators",
          fr: "Couverture selon différents dénominateurs",
          pt: "Cobertura com base em diferentes denominadores",
        },
        subCaption: {
          en: "DATE_RANGE",
          fr: "PLAGE_DE_DATES",
          pt: "INTERVALO_DE_DATAS",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnote: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
];

export const metric: MetricDefinitionGithub = {
  id: "m4a-02-01",
  resultsObjectId: "M5_combined_results_national.csv",
  valueProps: ["value"],
  valueFunc: "AVG",
  valueLabelReplacements: DENOMINATOR_VALUE_LABELS,
  label: {
    en: "Coverage estimated with different denominators",
    fr: "Couverture estimée avec différents dénominateurs",
    pt: "Cobertura estimada com diferentes denominadores",
  },
  variantLabel: {
    en: "National",
    fr: "National",
    pt: "Nacional",
  },
  requiredDisaggregationOptions: [
    "denominator_best_or_survey",
    "indicator_common_id",
    "year",
  ],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en:
        "Coverage estimates calculated using alternative denominator sources, plus survey benchmarks at national level.",
      fr:
        "Estimations de couverture calculées utilisant des sources de dénominateur alternatives, plus repères d'enquête au niveau national.",
      pt:
        "Estimativas de cobertura calculadas utilizando fontes de denominador alternativas, mais os valores de referência dos inquéritos ao nível nacional.",
    },
    methodology: {
      en:
        "AVG of coverage calculated as HMIS numerators divided by each available denominator type. Includes 'best' denominator (selected to minimize survey error) and 'survey' (original survey estimate). Enables comparison of how denominator choice affects coverage estimates.",
      fr:
        "Moyenne de la couverture calculée comme numérateurs HMIS divisés par chaque type de dénominateur disponible. Inclut le 'meilleur' dénominateur et l'estimé d'enquête.",
      pt:
        "Média (AVG) da cobertura calculada como os numeradores HMIS divididos por cada tipo de denominador disponível. Inclui o 'melhor' denominador e a estimativa de inquérito.",
    },
    interpretation: {
      en:
        "Large variation across denominator types indicates denominator uncertainty. Coverage estimates should be similar to survey values when using appropriate denominators. Coverage >100% suggests denominator underestimation or HMIS over-reporting.",
      fr:
        "Une grande variation entre types de dénominateur indique une incertitude. Les estimations de couverture devraient être similaires aux valeurs d'enquête avec des dénominateurs appropriés.",
      pt:
        "Uma grande variação entre os tipos de denominador indica incerteza. As estimativas de cobertura deveriam ser semelhantes aos valores dos inquéritos quando se utilizam denominadores adequados.",
    },
    typicalRange: {
      en:
        "0-100%. Denominators producing >100% coverage are likely inappropriate or require data quality investigation.",
      fr:
        "0-100%. Les dénominateurs produisant >100% de couverture sont probablement inappropriés ou nécessitent une investigation de qualité.",
      pt:
        "0-100%. Os denominadores que produzem coberturas superiores a 100% são provavelmente inadequados ou exigem uma investigação da qualidade dos dados.",
    },
    caveats: {
      en:
        "Denominator selection is subjective and affects results. The 'best' denominator minimizes squared error against surveys but may not be appropriate for all analytical purposes.",
      fr:
        "La sélection du dénominateur est subjective et affecte les résultats. Le 'meilleur' dénominateur minimise l'erreur quadratique mais peut ne pas être approprié pour tous les objectifs.",
      pt:
        "A seleção do denominador é subjetiva e afeta os resultados. O 'melhor' denominador minimiza o erro quadrático em relação aos inquéritos, mas pode não ser adequado para todos os fins analíticos.",
    },
    disaggregationGuidance: {
      en:
        "Always disaggregate by denominator_best_or_survey, indicator_common_id, and year (all required). Compare 'best' vs 'survey' to validate denominator selection. Visualize all denominator types to assess range of plausible coverage estimates.",
      fr:
        "Toujours désagréger par denominator_best_or_survey, indicator_common_id et year (tous requis). Comparer 'meilleur' vs 'enquête' pour valider la sélection du dénominateur.",
      pt:
        "Desagregar sempre por denominator_best_or_survey, indicator_common_id e year (todos obrigatórios). Comparar 'melhor' vs 'inquérito' para validar a seleção do denominador.",
    },
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  catalogExpressionEvaluation: null,
  vizPresets,
};
