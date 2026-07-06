import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { CF_DIVERGING_5_10_20 } from "../../.validation/cf_presets.ts";

export const vizPresets: VizPreset[] = [
  ///////////////////////////////////////////////
  //  ________         __        __            //
  // /        |       /  |      /  |           //
  // $$$$$$$$/______  $$ |____  $$ |  ______   //
  //    $$ | /      \ $$      \ $$ | /      \  //
  //    $$ | $$$$$$  |$$$$$$$  |$$ |/$$$$$$  | //
  //    $$ | /    $$ |$$ |  $$ |$$ |$$    $$ | //
  //    $$ |/$$$$$$$ |$$ |__$$ |$$ |$$$$$$$$/  //
  //    $$ |$$    $$ |$$    $$/ $$ |$$       | //
  //    $$/  $$$$$$$/ $$$$$$$/  $$/  $$$$$$$/  //
  //                                           //
  ///////////////////////////////////////////////
  {
    id: "disruption-differences-table-single-admin-area-2-multiple-admin-area-3",
    label: {
      en: "Service disruption by Admin Area 3",
      fr: "Perturbation des services par Zone administrative 3",
      pt: "Perturbação dos serviços por Área administrativa 3",
    },
    description: {
      en: "Table showing percentage difference between actual vs expected service volume for a single Admin Area 2, with multiple Admin Areas 3 and multiple indicators",
      fr: "Tableau montrant la différence en pourcentage entre le volume de services réel et attendu pour une unique Zone administrative 2, avec plusieurs Zones administratives 3 et plusieurs indicateurs",
      pt: "Tabela que mostra a diferença percentual entre o volume de serviços real e esperado para uma única Área administrativa 2, com várias Áreas administrativas 3 e vários indicadores",
    },
    allowedFilters: ["indicator_common_id", "admin_area_3"],
    config: {
      d: {
        type: "table",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "col",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "col",
          },
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "replicant",
          },
          {
            disOpt: "admin_area_3",
            disDisplayOpt: "row",
          },
        ],
        selectedReplicantValue: "",
        filterBy: [],
        periodFilter: {
          filterType: "last_n_calendar_quarters",
          nQuarters: 1,
        },
      },
      s: {
        ...CF_DIVERGING_5_10_20,
      },
      t: {
        caption: null,
        captionRelFontSize: null,
        subCaption: {
          en: "DATE_RANGE",
          fr: "PLAGE_DE_DATES",
          pt: "INTERVALO_DE_DATAS",
        },
        subCaptionRelFontSize: null,
        footnote: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
    createDefaultVisualizationOnInstall: "82a91435-b090-4bfe-8477-fe864095155e",
  },
  ///////////////////////////////////////
  //  __       __                      //
  // /  \     /  |                     //
  // $$  \   /$$ |  ______    ______   //
  // $$$  \ /$$$ | /      \  /      \  //
  // $$$$  /$$$$ | $$$$$$  |/$$$$$$  | //
  // $$ $$ $$/$$ | /    $$ |$$ |  $$ | //
  // $$ |$$$/ $$ |/$$$$$$$ |$$ |__$$ | //
  // $$ | $/  $$ |$$    $$ |$$    $$/  //
  // $$/      $$/  $$$$$$$/ $$$$$$$/   //
  //                        $$ |       //
  //                        $$ |       //
  //                        $$/        //
  //                                   //
  ///////////////////////////////////////
  {
    id: "disruption-differences-map-admin-area-3",
    label: {
      en: "Service disruption map by Admin Area 3",
      fr: "Carte des perturbations par Zone administrative 3",
      pt: "Mapa de perturbação dos serviços por Área administrativa 3",
    },
    description: {
      en: "Map showing percentage difference between actual vs expected service volume by Admin Area 3",
      fr: "Carte montrant la différence en pourcentage entre le volume de services réel et attendu par Zone administrative 3",
      pt: "Mapa que mostra a diferença percentual entre o volume de serviços real e esperado por Área administrativa 3",
    },
    allowedFilters: [],
    config: {
      d: {
        type: "map",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "cell",
        disaggregateBy: [
          {
            disOpt: "admin_area_3",
            disDisplayOpt: "mapArea",
          },
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "replicant",
          },
        ],
        selectedReplicantValue: "anc1",
        filterBy: [],
        periodFilter: {
          filterType: "last_n_calendar_quarters",
          nQuarters: 1,
        },
      },
      s: {
        ...CF_DIVERGING_5_10_20,
      },
      t: {
        caption: {
          en: "Service disruption",
          fr: "Perturbation des services",
          pt: "Perturbação dos serviços",
        },
        subCaption: {
          en: "Percentage difference between actual and expected service volume by Admin Area 3, DATE_RANGE",
          fr: "Différence en pourcentage entre le volume de services réel et attendu par Zone administrative 3, PLAGE_DE_DATES",
          pt: "Diferença percentual entre o volume de serviços real e esperado por Área administrativa 3, INTERVALO_DE_DATAS",
        },
        footnote: {
          en: "Negative values indicate service delivery below expected levels. Expected values are based on pre-disruption trends.",
          fr: "Les valeurs négatives indiquent une prestation de services inférieure aux niveaux attendus. Les valeurs attendues sont basées sur les tendances pré-perturbation.",
          pt: "Os valores negativos indicam uma prestação de serviços abaixo dos níveis esperados. Os valores esperados baseiam-se nas tendências anteriores à perturbação.",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
    createDefaultVisualizationOnInstall: "732d4f78-e53c-4c44-ac9f-234f616baa2c",
  },
];

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
  id: "m3-04-02",
  resultsObjectId: "M3_disruptions_analysis_admin_area_3.csv",
  label: {
    en: "Difference between actual and expected service volume (%)",
    fr: "Différence entre le volume de services réel et attendu (%)",
    pt: "Diferença entre o volume de serviços real e esperado (%)",
  },
  variantLabel: {
    en: "Admin area 3",
    fr: "Zone administrative 3",
    pt: "Área administrativa 3",
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
  requiredDisaggregationOptions: ["indicator_common_id", "admin_area_3"],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Percentage difference between actual and expected service volumes at district level.",
      fr: "Différence en pourcentage entre volumes réels et attendus au niveau du district.",
      pt: "Diferença percentual entre os volumes de serviços reais e esperados ao nível do distrito.",
    },
    methodology: {
      en: "(actual - expected) / expected at district level. Quantifies local service delivery performance.",
      fr: "(réel - attendu) / attendu au niveau du district. Quantifie la performance locale de prestation.",
      pt: "(real - esperado) / esperado ao nível do distrito. Quantifica o desempenho local da prestação de serviços.",
    },
    interpretation: {
      en: "Enables district-level performance monitoring. Prioritize districts with largest negative deviations.",
      fr: "Permet la surveillance de la performance au niveau du district. Prioriser les districts avec les plus grands écarts négatifs.",
      pt: "Permite a monitorização do desempenho ao nível do distrito. Priorizar os distritos com os maiores desvios negativos.",
    },
    typicalRange: {
      en: "±10-40% variation common at district level; >40% deviation warrants investigation.",
      fr: "Variation de ±10-40% commune au niveau du district; écart >40% nécessite investigation.",
      pt: "Variação de ±10-40% comum ao nível do distrito; um desvio >40% justifica investigação.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id and admin_area_3 (both required). Compare with admin area 2 for context.",
      fr: "Toujours désagréger par indicator_common_id et admin_area_3 (tous deux requis). Comparer avec la zone administrative 2 pour le contexte.",
      pt: "Desagregar sempre por indicator_common_id e admin_area_3 (ambos obrigatórios). Comparar com a área administrativa 2 para contexto.",
    },
    caveats: null,
  },
  importantNotes: null,
  hide: false,
  vizPresets,
};
