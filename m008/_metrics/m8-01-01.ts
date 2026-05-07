import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [];

export const metric: MetricDefinitionGithub = {
  id: "m8-01-01",
  hide: false,
  resultsObjectId: "M8_output_scorecard.csv",
  valueProps: ["value"],
  valueFunc: "identity",
  valueLabelReplacements: {},
  label: {
    en: "Scorecard",
    fr: "Scorecard",
  },
  requiredDisaggregationOptions: ["indicator_common_id"],
  formatAs: "percent",
  postAggregationExpression: {
    ingredientValues: [
      { prop: "numerator", func: "SUM" },
      { prop: "denominator", func: "SUM" },
    ],
    expression: "value = numerator / denominator",
  },
  aiDescription: {
    summary: {
      en: "Health sector scorecard indicators computed from the calculated indicators catalog.",
      fr: "Indicateurs de scorecard du secteur santé calculés à partir du catalogue des indicateurs calculés.",
    },
    methodology: {
      en: "Each indicator is computed as numerator/denominator. Denominators are either other HMIS indicators or population-based (annual population fraction scaled to monthly).",
      fr: "Chaque indicateur est calculé comme numérateur/dénominateur. Les dénominateurs sont soit d'autres indicateurs HMIS, soit basés sur la population (fraction annuelle de population mise à l'échelle mensuelle).",
    },
    interpretation: {
      en: "Values represent coverage rates or ratios. Thresholds for traffic-light coloring come from the calculated indicators catalog.",
      fr: "Les valeurs représentent des taux de couverture ou des ratios. Les seuils pour la coloration feu tricolore proviennent du catalogue des indicateurs calculés.",
    },
    typicalRange: {
      en: "Varies by indicator - see catalog for thresholds.",
      fr: "Varie selon l'indicateur - voir le catalogue pour les seuils.",
    },
    disaggregationGuidance: {
      en: "Can be disaggregated by admin area (2, 3, or 4) and time period (month, quarter, year). Aggregation uses SUM(numerator)/SUM(denominator).",
      fr: "Peut être désagrégé par zone administrative (2, 3 ou 4) et période de temps (mois, trimestre, année). L'agrégation utilise SUM(numérateur)/SUM(dénominateur).",
    },
    caveats: null,
  },
  variantLabel: null,
  importantNotes: null,
  vizPresets,
};
