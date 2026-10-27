import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { DENOMINATOR_VALUE_LABELS } from "../_denominator_labels.ts";

export const vizPresets: VizPreset[] = [];

export const metric: MetricDefinitionGithub = {
  id: "m4a-01-03",
  resultsObjectId: "M5_denominators_admin3.csv",
  valueProps: ["value"],
  valueFunc: "AVG",
  valueLabelReplacements: DENOMINATOR_VALUE_LABELS,
  label: {
    en: "Denominator values",
    fr: "Valeurs de dénominateur",
    pt: "Valores de denominador",
  },
  variantLabel: {
    en: "Admin area 3",
    fr: "Zone administrative 3",
    pt: "Área administrativa 3",
  },
  requiredDisaggregationOptions: ["denominator", "admin_area_3", "year"],
  formatAs: "number",
  aiDescription: {
    summary: {
      en:
        "District-level population denominators for coverage calculation at admin area 3 level.",
      fr:
        "Dénominateurs de population au niveau du district pour le calcul de couverture au niveau de la zone administrative 3.",
      pt:
        "Denominadores populacionais ao nível do distrito para o cálculo da cobertura ao nível da área administrativa 3.",
    },
    methodology: {
      en:
        "AVG of district-level denominators. Finest geographic resolution available. Only generated when ANALYSIS_LEVEL includes admin_area_3.",
      fr:
        "Moyenne des dénominateurs au niveau du district. Résolution géographique la plus fine disponible. Généré uniquement lorsque ANALYSIS_LEVEL inclut admin_area_3.",
      pt:
        "Média (AVG) dos denominadores ao nível do distrito. Resolução geográfica mais fina disponível. Apenas gerado quando ANALYSIS_LEVEL inclui admin_area_3.",
    },
    interpretation: {
      en:
        "District denominators most uncertain due to small sample sizes and population mobility. Use with caution for micro-level targeting.",
      fr:
        "Les dénominateurs de district sont les plus incertains en raison de petites tailles d'échantillon et mobilité de la population.",
      pt:
        "Os denominadores de distrito são os mais incertos devido às pequenas dimensões das amostras e à mobilidade da população.",
    },
    typicalRange: {
      en:
        "Varies by district size. Smaller than admin area 2 denominators, proportional to population.",
      fr:
        "Varie selon la taille du district. Plus petit que les dénominateurs de zone administrative 2, proportionnel à la population.",
      pt:
        "Varia consoante a dimensão do distrito. Menor do que os denominadores da área administrativa 2, proporcional à população.",
    },
    disaggregationGuidance: {
      en:
        "Always disaggregate by denominator, admin_area_3, and year (all required). Interpret with caution due to denominator uncertainty at this level.",
      fr:
        "Toujours désagréger par denominator, admin_area_3 et year (tous requis). Interpréter avec prudence en raison de l'incertitude du dénominateur à ce niveau.",
      pt:
        "Desagregar sempre por denominator, admin_area_3 e year (todos obrigatórios). Interpretar com prudência devido à incerteza do denominador a este nível.",
    },
    caveats: null,
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  catalogExpressionEvaluation: null,
  vizPresets,
};
