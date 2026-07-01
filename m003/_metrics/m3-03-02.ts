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
    },
    description: {
      en: "Table showing percentage difference between actual vs expected service volume, with multiple Admin Areas 2 and multiple indicators",
      fr: "Tableau montrant la différence en pourcentage entre le volume de services réel et attendu, avec plusieurs Zones administratives 2 et plusieurs indicateurs",
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
    },
    description: {
      en: "Map showing percentage difference between actual vs expected service volume by Admin Area 2",
      fr: "Carte montrant la différence en pourcentage entre le volume de services réel et attendu par Zone administrative 2",
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
        },
        subCaption: {
          en: "Percentage difference between actual and expected service volume by Admin Area 2, DATE_RANGE",
          fr: "Différence en pourcentage entre le volume de services réel et attendu par Zone administrative 2, PLAGE_DE_DATES",
        },
        footnote: {
          en: "Negative values indicate service delivery below expected levels. Expected values are based on pre-disruption trends.",
          fr: "Les valeurs négatives indiquent une prestation de services inférieure aux niveaux attendus. Les valeurs attendues sont basées sur les tendances pré-perturbation.",
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
  },
  variantLabel: {
    en: "Admin area 2",
    fr: "Zone administrative 2",
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
    },
    methodology: {
      en: "(actual - expected) / expected at subnational level. Quantifies service delivery gaps or surpluses for each province/state.",
      fr: "(réel - attendu) / attendu au niveau sous-national. Quantifie les écarts ou excédents de prestation de services.",
    },
    interpretation: {
      en: "Negative values indicate regional shortfalls; positive values may indicate improved service delivery or data issues. Compare across regions to identify areas needing support.",
      fr: "Les valeurs négatives indiquent des déficits régionaux; les valeurs positives peuvent indiquer une amélioration de la prestation de services ou des problèmes de données. Comparer entre les régions pour identifier les zones nécessitant un soutien.",
    },
    typicalRange: {
      en: "±10-30% variation common; >30% deviation warrants investigation.",
      fr: "Variation de ±10-30% commune; écart >30% nécessite investigation.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id and admin_area_2 (both required). Time series reveals disruption patterns over time.",
      fr: "Toujours désagréger par indicator_common_id et admin_area_2 (tous deux requis). Les séries temporelles révèlent les modèles de perturbation au fil du temps.",
    },
    caveats: null,
  },
  importantNotes: null,
  hide: false,
  vizPresets,
};
