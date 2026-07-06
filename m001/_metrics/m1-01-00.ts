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
    pt: "Número de registos",
  },
  requiredDisaggregationOptions: [],
  formatAs: "number",
  aiDescription: {
    summary: {
      en: "Count of facility-month-indicator records in the dataset.",
      fr: "Nombre d'enregistrements établissement-mois-indicateur dans le jeu de données.",
      pt: "Contagem de registos unidade sanitária-mês-indicador no conjunto de dados.",
    },
    methodology: {
      en: "COUNT of unique facility-indicator-period combinations in the database.",
      fr: "Décompte des combinaisons uniques établissement-indicateur-période dans la base de données.",
      pt: "CONTAGEM das combinações únicas unidade sanitária-indicador-período na base de dados.",
    },
    interpretation: {
      en: "Higher counts indicate more complete reporting coverage. Low counts may indicate data gaps or limited facility participation.",
      fr: "Des valeurs plus élevées indiquent une couverture de déclaration plus complète. Des valeurs basses peuvent indiquer des lacunes de données.",
      pt: "Valores mais elevados indicam uma cobertura de notificação mais completa. Valores baixos podem indicar lacunas nos dados ou uma participação limitada das unidades sanitárias.",
    },
    typicalRange: {
      en: "Varies by country size and time period selected.",
      fr: "Varie selon la taille du pays et la période sélectionnée.",
      pt: "Varia consoante a dimensão do país e o período de tempo selecionado.",
    },
    disaggregationGuidance: {
      en: "Disaggregate by admin_area to compare regional completeness. Use indicator_common_id to see which services have better reporting.",
      fr: "Désagréger par zone administrative pour comparer la complétude régionale.",
      pt: "Desagregar por admin_area para comparar a completude regional. Utilizar indicator_common_id para ver quais os serviços com melhor notificação.",
    },
    caveats: null,
  },
  variantLabel: null,
  postAggregationExpression: null,
  importantNotes: null,
  vizPresets,
};
