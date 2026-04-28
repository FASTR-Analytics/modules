import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [];

///////////////////////////////////////////////////////////////
//  __       __              __                __            //
// /  \     /  |            /  |              /  |           //
// $$  \   /$$ |  ______   _$$ |_     ______  $$/   _______  //
// $$$  \ /$$$ | /      \ / $$   |   /      \ /  | /       | //
// $$$$  /$$$$ |/$$$$$$  |$$$$$$/   /$$$$$$  |$$ |/$$$$$$$/  //
// $$ $$ $$/$$ |$$    $$ |  $$ | __ $$ |  $$/ $$ |$$ |       //
// $$ |$$$/ $$ |$$$$$$$$/   $$ |/  |$$ |      $$ |$$ \_____  //
// $$ | $/  $$ |$$       |  $$  $$/ $$ |      $$ |$$       | //
// $$/      $$/  $$$$$$$/    $$$$/  $$/       $$/  $$$$$$$/  //
//                                                           //
///////////////////////////////////////////////////////////////

export const metric: MetricDefinitionGithub = {
  id: "m1-05-01",
  hide: true,
  resultsObjectId: "M1_output_outlier_list.csv",
  valueProps: ["count"],
  valueFunc: "SUM",
  valueLabelReplacements: {
    dqa_score: "Indicator outliers",
  },
  label: {
    en: "Indicator outliers",
    fr: "Valeurs aberrantes des indicateurs",
  },
  requiredDisaggregationOptions: [],
  formatAs: "number",
  aiDescription: {
    summary: {
      en: "Total number of outlier data points identified across all facilities.",
      fr: "Nombre total de points de données aberrants identifiés dans tous les établissements.",
    },
    methodology: {
      en: "SUM of outlier counts from the outlier list. Each outlier represents a facility-period-indicator combination flagged by MAD-based detection.",
      fr: "Somme des comptes de valeurs aberrantes de la liste des valeurs aberrantes. Chaque valeur aberrante représente une combinaison établissement-période-indicateur signalée par détection basée sur MAD.",
    },
    interpretation: {
      en: "Higher counts indicate more widespread data quality issues. Use alongside outlier proportion to understand both prevalence and severity of quality problems.",
      fr: "Des comptes plus élevés indiquent des problèmes de qualité des données plus répandus. Utiliser avec la proportion de valeurs aberrantes pour comprendre la prévalence et la gravité.",
    },
    typicalRange: {
      en: "Varies widely by dataset size and quality.",
      fr: "Varie considérablement selon la taille et la qualité du jeu de données.",
    },
    caveats: {
      en: "Absolute counts are less informative than proportions when comparing across different time periods or regions with varying reporting volumes.",
      fr: "Les comptes absolus sont moins informatifs que les proportions lors de comparaisons entre différentes périodes ou régions.",
    },
    disaggregationGuidance: {
      en: "Disaggregate by indicator_common_id to identify problem indicators. Use admin_area to find regions with highest outlier counts. Combine with m1-01-01 for context.",
      fr: "Désagréger par indicator_common_id pour identifier les indicateurs problématiques. Utiliser admin_area pour trouver les régions avec le plus de valeurs aberrantes.",
    },
  },
  variantLabel: null,
  postAggregationExpression: null,
  importantNotes: null,
  vizPresets,
};
