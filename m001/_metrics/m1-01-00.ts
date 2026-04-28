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
  id: "m1-01-00",
  hide: true,
  resultsObjectId: "M1_output_outliers.csv",
  valueProps: ["facility_id"],
  valueFunc: "COUNT",
  valueLabelReplacements: {},
  label: {
    en: "Number of records",
    fr: "Nombre d'enregistrements",
  },
  requiredDisaggregationOptions: [],
  formatAs: "number",
  aiDescription: {
    summary: {
      en: "Count of facility-month-indicator records in the dataset.",
      fr: "Nombre d'enregistrements établissement-mois-indicateur dans le jeu de données.",
    },
    methodology: {
      en: "COUNT of unique facility-indicator-period combinations in the database.",
      fr: "Décompte des combinaisons uniques établissement-indicateur-période dans la base de données.",
    },
    interpretation: {
      en: "Higher counts indicate more complete reporting coverage. Low counts may indicate data gaps or limited facility participation.",
      fr: "Des valeurs plus élevées indiquent une couverture de déclaration plus complète. Des valeurs basses peuvent indiquer des lacunes de données.",
    },
    typicalRange: {
      en: "Varies by country size and time period selected.",
      fr: "Varie selon la taille du pays et la période sélectionnée.",
    },
    disaggregationGuidance: {
      en: "Disaggregate by admin_area to compare regional completeness. Use indicator_common_id to see which services have better reporting.",
      fr: "Désagréger par zone administrative pour comparer la complétude régionale.",
    },
    caveats: null,
  },
  variantLabel: null,
  postAggregationExpression: null,
  importantNotes: null,
  vizPresets,
};
