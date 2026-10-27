import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { DENOMINATOR_VALUE_LABELS } from "../_denominator_labels.ts";

export const vizPresets: VizPreset[] = [];

export const metric: MetricDefinitionGithub = {
  id: "m4a-02-03",
  resultsObjectId: "M5_combined_results_admin3.csv",
  valueProps: ["value"],
  valueFunc: "AVG",
  valueLabelReplacements: DENOMINATOR_VALUE_LABELS,
  label: {
    en: "Coverage estimated with different denominators",
    fr: "Couverture estimée avec différents dénominateurs",
    pt: "Cobertura estimada com diferentes denominadores",
  },
  variantLabel: {
    en: "Admin area 3",
    fr: "Zone administrative 3",
    pt: "Área administrativa 3",
  },
  requiredDisaggregationOptions: [
    "denominator_best_or_survey",
    "admin_area_3",
    "indicator_common_id",
    "year",
  ],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en:
        "District-level coverage estimates using alternative denominators at admin area 3 level.",
      fr:
        "Estimations de couverture au niveau du district utilisant des dénominateurs alternatifs au niveau de la zone administrative 3.",
      pt:
        "Estimativas de cobertura ao nível do distrito utilizando denominadores alternativos ao nível da área administrativa 3.",
    },
    methodology: {
      en:
        "AVG of coverage with different denominator types at district level. Highest geographic resolution for coverage estimation.",
      fr:
        "Moyenne de la couverture avec différents types de dénominateur au niveau du district. Plus haute résolution géographique pour l'estimation de couverture.",
      pt:
        "Média (AVG) da cobertura com diferentes tipos de denominador ao nível do distrito. Resolução geográfica mais elevada para a estimativa de cobertura.",
    },
    interpretation: {
      en:
        "District-level denominator uncertainty typically highest. Use to understand micro-level coverage patterns but interpret with caution.",
      fr:
        "L'incertitude du dénominateur au niveau du district est généralement la plus élevée. Utiliser pour comprendre les modèles de couverture micro-niveau avec prudence.",
      pt:
        "A incerteza do denominador ao nível do distrito é geralmente a mais elevada. Utilizar para compreender os padrões de cobertura ao nível micro, mas interpretar com prudência.",
    },
    typicalRange: {
      en:
        "0-100%. Wide variation across denominators common due to small sample sizes and denominator uncertainty.",
      fr:
        "0-100%. Variation large entre dénominateurs commune en raison de petites tailles d'échantillon et incertitude du dénominateur.",
      pt:
        "0-100%. É comum uma ampla variação entre denominadores devido às pequenas dimensões das amostras e à incerteza do denominador.",
    },
    disaggregationGuidance: {
      en:
        "Always disaggregate by denominator_best_or_survey, admin_area_3, indicator_common_id, and year (all required). High uncertainty at this level - consider aggregating to admin area 2.",
      fr:
        "Toujours désagréger par denominator_best_or_survey, admin_area_3, indicator_common_id et year (tous requis). Incertitude élevée à ce niveau - considérer l'agrégation à la zone administrative 2.",
      pt:
        "Desagregar sempre por denominator_best_or_survey, admin_area_3, indicator_common_id e year (todos obrigatórios). Incerteza elevada a este nível - considerar a agregação à área administrativa 2.",
    },
    caveats: null,
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  catalogExpressionEvaluation: null,
  vizPresets,
};
