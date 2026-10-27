import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { DENOMINATOR_VALUE_LABELS } from "../_denominator_labels.ts";

export const vizPresets: VizPreset[] = [];

export const metric: MetricDefinitionGithub = {
  id: "m4a-01-02",
  resultsObjectId: "M5_denominators_admin2.csv",
  valueProps: ["value"],
  valueFunc: "AVG",
  valueLabelReplacements: DENOMINATOR_VALUE_LABELS,
  label: {
    en: "Denominator values",
    fr: "Valeurs de dénominateur",
    pt: "Valores de denominador",
  },
  variantLabel: {
    en: "Admin area 2",
    fr: "Zone administrative 2",
    pt: "Área administrativa 2",
  },
  requiredDisaggregationOptions: ["denominator", "admin_area_2", "year"],
  formatAs: "number",
  aiDescription: {
    summary: {
      en:
        "Subnational population denominators for coverage calculation at admin area 2 level.",
      fr:
        "Dénominateurs de population sous-nationaux pour le calcul de couverture au niveau de la zone administrative 2.",
      pt:
        "Denominadores populacionais subnacionais para o cálculo da cobertura ao nível da área administrativa 2.",
    },
    methodology: {
      en:
        "AVG of denominators derived from subnational HMIS data and survey estimates. UNWPP denominators not available at this level. Uses national-level methodology adapted to subnational populations.",
      fr:
        "Moyenne des dénominateurs dérivés des données HMIS sous-nationales et estimations d'enquête. Dénominateurs UNWPP non disponibles à ce niveau.",
      pt:
        "Média (AVG) dos denominadores derivados dos dados HMIS subnacionais e das estimativas de inquérito. Os denominadores UNWPP não estão disponíveis a este nível.",
    },
    interpretation: {
      en:
        "Subnational denominators enable regional coverage monitoring. Compare across admin areas to assess population estimate consistency. Larger regional variation suggests denominator uncertainty.",
      fr:
        "Les dénominateurs sous-nationaux permettent la surveillance de la couverture régionale. Comparer entre zones administratives pour évaluer la cohérence des estimations.",
      pt:
        "Os denominadores subnacionais permitem a monitorização da cobertura regional. Comparar entre áreas administrativas para avaliar a coerência das estimativas populacionais.",
    },
    typicalRange: {
      en:
        "Varies by region size and population. Proportional to national denominators.",
      fr:
        "Varie selon la taille de la région et la population. Proportionnel aux dénominateurs nationaux.",
      pt:
        "Varia consoante a dimensão e a população da região. Proporcional aos denominadores nacionais.",
    },
    disaggregationGuidance: {
      en:
        "Always disaggregate by denominator, admin_area_2, and year (all required). Compare denominator sources to assess regional data quality.",
      fr:
        "Toujours désagréger par denominator, admin_area_2 et year (tous requis). Comparer les sources de dénominateur pour évaluer la qualité des données régionales.",
      pt:
        "Desagregar sempre por denominator, admin_area_2 e year (todos obrigatórios). Comparar as fontes de denominador para avaliar a qualidade dos dados regionais.",
    },
    caveats: null,
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  catalogExpressionEvaluation: null,
  vizPresets,
};
