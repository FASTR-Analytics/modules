import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { DENOMINATOR_VALUE_LABELS } from "../_denominator_labels.ts";

export const vizPresets: VizPreset[] = [];

export const metric: MetricDefinitionGithub = {
  id: "m4a-02-02",
  resultsObjectId: "M5_combined_results_admin2.csv",
  valueProps: ["value"],
  valueFunc: "AVG",
  valueLabelReplacements: DENOMINATOR_VALUE_LABELS,
  label: {
    en: "Coverage estimated with different denominators",
    fr: "Couverture estimée avec différents dénominateurs",
    pt: "Cobertura estimada com diferentes denominadores",
  },
  variantLabel: {
    en: "Admin area 2",
    fr: "Zone administrative 2",
    pt: "Área administrativa 2",
  },
  requiredDisaggregationOptions: [
    "denominator_best_or_survey",
    "admin_area_2",
    "indicator_common_id",
    "year",
  ],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en:
        "Subnational coverage estimates using alternative denominators at admin area 2 level.",
      fr:
        "Estimations de couverture sous-nationales utilisant des dénominateurs alternatifs au niveau de la zone administrative 2.",
      pt:
        "Estimativas de cobertura subnacionais utilizando denominadores alternativos ao nível da área administrativa 2.",
    },
    methodology: {
      en:
        "AVG of coverage with different denominator types at subnational level. National-only denominators (e.g., UNWPP) replaced with subnational alternatives.",
      fr:
        "Moyenne de la couverture avec différents types de dénominateur au niveau sous-national. Dénominateurs nationaux uniquement remplacés par alternatives sous-nationales.",
      pt:
        "Média (AVG) da cobertura com diferentes tipos de denominador ao nível subnacional. Os denominadores apenas nacionais (por exemplo, UNWPP) são substituídos por alternativas subnacionais.",
    },
    interpretation: {
      en:
        "Enables assessment of denominator uncertainty at regional level. Useful for understanding geographic variation in coverage and denominator quality.",
      fr:
        "Permet l'évaluation de l'incertitude du dénominateur au niveau régional. Utile pour comprendre la variation géographique de la couverture.",
      pt:
        "Permite avaliar a incerteza do denominador ao nível regional. Útil para compreender a variação geográfica da cobertura e a qualidade do denominador.",
    },
    typicalRange: {
      en:
        "0-100%. Variation across denominators indicates uncertainty in regional coverage estimates.",
      fr:
        "0-100%. La variation entre dénominateurs indique une incertitude dans les estimations de couverture régionale.",
      pt:
        "0-100%. A variação entre denominadores indica incerteza nas estimativas de cobertura regional.",
    },
    disaggregationGuidance: {
      en:
        "Always disaggregate by denominator_best_or_survey, admin_area_2, indicator_common_id, and year (all required). Compare denominator types to assess regional estimate uncertainty.",
      fr:
        "Toujours désagréger par denominator_best_or_survey, admin_area_2, indicator_common_id et year (tous requis). Comparer les types de dénominateur pour évaluer l'incertitude des estimations régionales.",
      pt:
        "Desagregar sempre por denominator_best_or_survey, admin_area_2, indicator_common_id e year (todos obrigatórios). Comparar os tipos de denominador para avaliar a incerteza das estimativas regionais.",
    },
    caveats: null,
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  catalogExpressionEvaluation: null,
  vizPresets,
};
