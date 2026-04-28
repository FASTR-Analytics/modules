import type { MetricDefinitionGithub, VizPreset } from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [];

export const metric: MetricDefinitionGithub = {
  id: "m3-02-02",
  resultsObjectId: "M3_disruptions_analysis_admin_area_1.csv",
  label: {
    en: "Difference between actual and expected service volume",
    fr: "Différence entre le volume de services réel et attendu",
  },
  variantLabel: {
    en: "National",
    fr: "National",
  },
  valueProps: ["pct_diff"],
  valueFunc: "identity",
  valueLabelReplacements: {
    pct_diff: "Percent difference",
  },
  postAggregationExpression: {
    ingredientValues: [
      {
        prop: "count_sum",
        func: "SUM",
      },
      {
        prop: "count_expect_sum",
        func: "SUM",
      },
    ],
    expression: "pct_diff = (count_sum - count_expect_sum)/count_expect_sum",
  },
  requiredDisaggregationOptions: ["indicator_common_id"],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Percentage difference between actual and model-predicted service volumes at national level.",
      fr: "Différence en pourcentage entre les volumes de services réels et prédits par modèle au niveau national.",
    },
    methodology: {
      en: "Calculated as (actual - expected) / expected. Positive values indicate volumes above expected; negative values indicate shortfalls. Expected volumes come from robust panel regression models accounting for trends and seasonality.",
      fr: "Calculé comme (réel - attendu) / attendu. Les valeurs positives indiquent des volumes au-dessus de l'attendu; les valeurs négatives indiquent des déficits.",
    },
    interpretation: {
      en: "Negative values indicate service disruptions or underperformance relative to expected levels. Values below -20% warrant investigation. Positive values may reflect data quality issues, special campaigns, increased demand, or genuine service improvements. Use alongside control chart flags for context.",
      fr: "Les valeurs négatives indiquent des perturbations de services ou une sous-performance par rapport aux niveaux attendus. Les valeurs inférieures à -20% nécessitent une investigation. Les valeurs positives peuvent refléter des problèmes de qualité des données, des campagnes spéciales, une demande accrue ou de véritables améliorations de services. Utiliser avec les indicateurs de carte de contrôle pour le contexte.",
    },
    typicalRange: {
      en: "±10% is within normal variation; ±10-30% indicates moderate disruption; >30% deviation suggests major disruption or data issues.",
      fr: "±10% est dans la variation normale; ±10-30% indique une perturbation modérée; >30% d'écart suggère une perturbation majeure ou des problèmes de données.",
    },
    caveats: {
      en: "Percentage differences can be misleading for indicators with small absolute volumes. Model predictions assume stable baseline patterns - structural changes in the health system may appear as disruptions.",
      fr: "Les différences en pourcentage peuvent être trompeuses pour les indicateurs avec de petits volumes absolus. Les prédictions du modèle supposent des modèles de base stables - les changements structurels du système de santé peuvent apparaître comme des perturbations.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id (required). Time series shows disruption evolution and recovery patterns. Consider using absolute differences (m3-02-01) alongside percentages for low-volume indicators.",
      fr: "Toujours désagréger par indicator_common_id (requis). Les séries temporelles montrent l'évolution de la perturbation et les modèles de récupération. Considérer l'utilisation des différences absolues (m3-02-01) aux côtés des pourcentages pour les indicateurs à faible volume.",
    },
  },
  importantNotes: null,
  hide: false,
  vizPresets,
};
