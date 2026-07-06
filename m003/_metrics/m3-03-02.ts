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
    id: "disruption-differences-table",
    label: {
      en: "Service disruption by Admin Area 2",
      fr: "Perturbation des services par Zone administrative 2",
      pt: "Perturbação dos serviços por Área administrativa 2",
    },
    description: {
      en: "Table showing percentage difference between actual vs expected service volume, with multiple Admin Areas 2 and multiple indicators",
      fr: "Tableau montrant la différence en pourcentage entre le volume de services réel et attendu, avec plusieurs Zones administratives 2 et plusieurs indicateurs",
      pt: "Tabela que mostra a diferença percentual entre o volume de serviços real e esperado, com várias Áreas administrativas 2 e vários indicadores",
    },
    allowedFilters: ["indicator_common_id", "admin_area_2"],
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
            disDisplayOpt: "row",
          },
        ],
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
    createDefaultVisualizationOnInstall: "b93b6411-9c6a-410c-93d1-f878f866b5c1",
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
    id: "disruption-differences-map-admin-area-2",
    label: {
      en: "Service disruption map by Admin Area 2",
      fr: "Carte des perturbations par Zone administrative 2",
      pt: "Mapa de perturbação dos serviços por Área administrativa 2",
    },
    description: {
      en: "Map showing percentage difference between actual vs expected service volume by Admin Area 2",
      fr: "Carte montrant la différence en pourcentage entre le volume de services réel et attendu par Zone administrative 2",
      pt: "Mapa que mostra a diferença percentual entre o volume de serviços real e esperado por Área administrativa 2",
    },
    allowedFilters: [],
    config: {
      d: {
        type: "map",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "cell",
        disaggregateBy: [
          {
            disOpt: "admin_area_2",
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
          en: "Percentage difference between actual and expected service volume by Admin Area 2, DATE_RANGE",
          fr: "Différence en pourcentage entre le volume de services réel et attendu par Zone administrative 2, PLAGE_DE_DATES",
          pt: "Diferença percentual entre o volume de serviços real e esperado por Área administrativa 2, INTERVALO_DE_DATAS",
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
    createDefaultVisualizationOnInstall: "5780650b-22e7-4d92-8239-8a1c51dbf0a5",
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
  id: "m3-03-02",
  resultsObjectId: "M3_disruptions_analysis_admin_area_2.csv",
  label: {
    en: "Difference between actual and expected service volume (%)",
    fr: "Différence entre le volume de services réel et attendu (%)",
    pt: "Diferença entre o volume de serviços real e esperado (%)",
  },
  variantLabel: {
    en: "Admin area 2",
    fr: "Zone administrative 2",
    pt: "Área administrativa 2",
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
  requiredDisaggregationOptions: ["indicator_common_id", "admin_area_2"],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Percentage difference between actual and expected service volumes at admin area 2 level.",
      fr: "Différence en pourcentage entre volumes réels et attendus au niveau de la zone administrative 2.",
      pt: "Diferença percentual entre os volumes de serviços reais e esperados ao nível da área administrativa 2.",
    },
    methodology: {
      en: "(actual - expected) / expected at subnational level. Quantifies service delivery gaps or surpluses for each province/state.",
      fr: "(réel - attendu) / attendu au niveau sous-national. Quantifie les écarts ou excédents de prestation de services.",
      pt: "(real - esperado) / esperado ao nível subnacional. Quantifica as lacunas ou excedentes na prestação de serviços para cada província/estado.",
    },
    interpretation: {
      en: "Negative values indicate regional shortfalls; positive values may indicate improved service delivery or data issues. Compare across regions to identify areas needing support.",
      fr: "Les valeurs négatives indiquent des déficits régionaux; les valeurs positives peuvent indiquer une amélioration de la prestation de services ou des problèmes de données. Comparer entre les régions pour identifier les zones nécessitant un soutien.",
      pt: "Os valores negativos indicam défices regionais; os valores positivos podem indicar uma melhoria na prestação de serviços ou problemas de dados. Comparar entre regiões para identificar as áreas que necessitam de apoio.",
    },
    typicalRange: {
      en: "±10-30% variation common; >30% deviation warrants investigation.",
      fr: "Variation de ±10-30% commune; écart >30% nécessite investigation.",
      pt: "Variação de ±10-30% comum; um desvio >30% justifica investigação.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id and admin_area_2 (both required). Time series reveals disruption patterns over time.",
      fr: "Toujours désagréger par indicator_common_id et admin_area_2 (tous deux requis). Les séries temporelles révèlent les modèles de perturbation au fil du temps.",
      pt: "Desagregar sempre por indicator_common_id e admin_area_2 (ambos obrigatórios). As séries temporais revelam os padrões de perturbação ao longo do tempo.",
    },
    caveats: null,
  },
  importantNotes: null,
  hide: false,
  vizPresets,
};
