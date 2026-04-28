import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { CF_90_80 } from "../../.validation/cf_presets.ts";

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
    id: "consistency-table",
    label: {
      en: "Internal consistency table",
      fr: "Tableau de cohérence interne",
    },
    description: {
      en: "Table showing consistency by ratio type and region",
      fr: "Tableau montrant la cohérence par type de ratio et région",
    },
    createDefaultVisualizationOnInstall: "cf5b8649-93c2-4bbe-8f2d-773f42ce8ec3",
    allowedFilters: ["ratio_type", "admin_area_2"],
    config: {
      d: {
        type: "table",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "col",
        disaggregateBy: [
          {
            disOpt: "ratio_type",
            disDisplayOpt: "col",
          },
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "row",
          },
        ],
        filterBy: [],
        periodFilter: {
          filterType: "last_n_months",
          nMonths: 12,
        },
      },
      s: {
        content: "lines",
        ...CF_90_80,
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Internal consistency",
          fr: "Cohérence interne",
        },
        subCaption: {
          en: "Percentage of sub-national areas meeting consistency benchmarks, DATE_RANGE",
          fr: "Pourcentage de zones sous-nationales atteignant les critères de cohérence, DATE_RANGE",
        },
        footnote: {
          en: "Internal consistency assesses the plausibility of reported data based on related indicators. Consistency metrics are approximate - depending on timing and seasonality, indicator definitions, and the nature of service delivery and reporting, values may be expected to sit outside plausible ranges. Indicators which are similar are expected to have roughy the same volume over the year (within a 30% margin). The data in this analysis is adjusted for outliers.",
          fr: "La cohérence interne évalue la plausibilité des données déclarées sur la base d'indicateurs liés. Les mesures de cohérence sont approximatives - selon le calendrier et la saisonnalité, les définitions des indicateurs, et la nature de la prestation de services et de la déclaration, les valeurs peuvent se situer en dehors des plages plausibles. Les indicateurs similaires sont censés avoir approximativement le même volume sur l'année (avec une marge de 30%). Les données de cette analyse sont ajustées pour les valeurs aberrantes.",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
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
    id: "consistency-map",
    label: {
      en: "Internal consistency map",
      fr: "Carte de cohérence interne",
    },
    description: {
      en: "Map showing consistency by ratio type and region",
      fr: "Carte montrant la cohérence par type de ratio et région",
    },
    createDefaultVisualizationOnInstall: "a6c091c0-8df5-4dc9-aa5f-857c697310d1",
    allowedFilters: [],
    config: {
      d: {
        type: "map",
        valuesDisDisplayOpt: "mapArea",
        disaggregateBy: [
          {
            disOpt: "ratio_type",
            disDisplayOpt: "replicant",
          },
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "mapArea",
          },
        ],
        selectedReplicantValue: "anc1_anc4",
        filterBy: [],
        periodFilter: {
          filterType: "last_n_months",
          nMonths: 12,
        },
      },
      s: {
        ...CF_90_80,
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Internal consistency",
          fr: "Cohérence interne",
        },
        subCaption: {
          en: "Percentage of sub-national areas meeting consistency benchmarks, DATE_RANGE",
          fr: "Pourcentage de zones sous-nationales atteignant les critères de cohérence, DATE_RANGE",
        },
        footnote: {
          en: "Internal consistency assesses the plausibility of reported data based on related indicators. Consistency metrics are approximate - depending on timing and seasonality, indicator definitions, and the nature of service delivery and reporting, values may be expected to sit outside plausible ranges. Indicators which are similar are expected to have roughy the same volume over the year (within a 30% margin). The data in this analysis is adjusted for outliers.",
          fr: "La cohérence interne évalue la plausibilité des données déclarées sur la base d'indicateurs liés. Les mesures de cohérence sont approximatives - selon le calendrier et la saisonnalité, les définitions des indicateurs, et la nature de la prestation de services et de la déclaration, les valeurs peuvent se situer en dehors des plages plausibles. Les indicateurs similaires sont censés avoir approximativement le même volume sur l'année (avec une marge de 30%). Les données de cette analyse sont ajustées pour les valeurs aberrantes.",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
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
  id: "m1-03-01",
  resultsObjectId: "M1_output_consistency_geo.csv",
  valueProps: ["sconsistency"],
  valueFunc: "AVG",
  valueLabelReplacements: {
    ratio_type: "Type of ratio being assessed",
    pair_anc: "ANC1 is larger than ANC4",
    pair_delivery: "Delivery is approximately equal to BCG",
    pair_pnc: "Delivery is larger than PNC1",
    pair_penta: "Penta 1 is larger than Penta 3",
  },
  label: {
    en: "Proportion of sub-national areas meeting consistency criteria",
    fr: "Proportion de zones sous-nationales répondant aux critères de cohérence",
  },
  requiredDisaggregationOptions: ["ratio_type"],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Proportion of sub-national areas where related indicators show logical consistency.",
      fr: "Proportion de zones sous-nationales où les indicateurs liés montrent une cohérence logique.",
    },
    methodology: {
      en: "AVG of sconsistency flag. Checks logical relationships between indicator pairs (e.g., ANC1 > ANC4, Penta1 > Penta3).",
      fr: "Moyenne du drapeau de cohérence. Vérifie les relations logiques entre paires d'indicateurs.",
    },
    interpretation: {
      en: "Higher values indicate better data quality. Low consistency suggests data entry errors or aggregation issues. Required disaggregation by ratio_type to see which consistency checks fail most often.",
      fr: "Des valeurs plus élevées indiquent une meilleure qualité des données. Une faible cohérence suggère des erreurs de saisie.",
    },
    typicalRange: {
      en: "90-100% is good; 70-90% acceptable; <70% needs investigation.",
      fr: "90-100% est bon; 70-90% acceptable; <70% nécessite investigation.",
    },
    caveats: {
      en: "Different ratio types have different expected pass rates. Some inconsistency may be clinically valid (e.g., vaccine stock-outs).",
      fr: "Différents types de ratios ont différents taux de réussite attendus.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by ratio_type as each consistency check has different implications. Use admin_area to find regions with systematic issues.",
      fr: "Toujours désagréger par ratio_type car chaque contrôle de cohérence a des implications différentes.",
    },
  },
  variantLabel: null,
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
