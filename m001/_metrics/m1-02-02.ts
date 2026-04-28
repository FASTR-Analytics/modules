import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

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
    id: "completeness-table",
    label: {
      en: "Completeness table by region",
      fr: "Tableau de complétude par région",
    },
    description: {
      en: "Table showing completeness by indicator and region",
      fr: "Tableau montrant la complétude par indicateur et région",
    },
    createDefaultVisualizationOnInstall: "c20f1672-edfc-4140-ae2c-09a30b50443a",
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
          filterType: "last_n_months",
          nMonths: 12,
        },
      },
      s: {
        content: "lines",
        cfMode: "thresholds",
        cfThresholdCutoffs: [0.8, 0.9],
        cfThresholdBuckets: [
          {
            color: "#F18989",
          },
          {
            color: "#F6D982",
          },
          {
            color: "#68C690",
          },
        ],
        cfThresholdDirection: "higher-is-better",
        cfThresholdNoDataColor: "#ffffff",
        cfScalePaletteKind: "preset",
        cfScalePalettePreset: "",
        cfScaleCustomFrom: "",
        cfScaleCustomMid: "",
        cfScaleCustomTo: "",
        cfScaleReverse: false,
        cfScaleSteps: 0,
        cfScaleDomainKind: "auto",
        cfScaleDomainMin: 0,
        cfScaleDomainMax: 1,
        cfScaleNoDataColor: "",
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Indicator Completeness",
          fr: "Complétude des indicateurs",
        },
        subCaption: {
          en: "Percentage of facility-months with complete data, DATE_RANGE",
          fr: "Pourcentage de mois-établissements avec des données complètes, DATE_RANGE",
        },
        footnote: {
          en: "Higher completeness improves the reliability of the data, especially when completeness is stable over time. Completeness is defined as the percentage of reporting facilities each month out of the total number of facilities expected to report. A facility is expected to report if it has reported any volume for each indicator anytime within a year. A high completeness does not indicate that the HMIS is representative of all service delivery in the country, as some services may not be delivered in facilities, or some facilities may not report.",
          fr: "Une complétude élevée améliore la fiabilité des données, surtout lorsqu'elle est stable dans le temps. La complétude est définie comme le pourcentage d'établissements déclarants chaque mois par rapport au nombre total d'établissements censés déclarer. Un établissement est censé déclarer s'il a déclaré un volume pour chaque indicateur à tout moment au cours de l'année. Une complétude élevée n'indique pas que le HMIS est représentatif de toute la prestation de services dans le pays, car certains services peuvent ne pas être fournis dans les établissements, ou certains établissements peuvent ne pas déclarer.",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //  ________  __                                                        __                      //
  // /        |/  |                                                      /  |                     //
  // $$$$$$$$/ $$/  _____  ____    ______    _______   ______    ______  $$/   ______    _______  //
  //    $$ |   /  |/     \/    \  /      \  /       | /      \  /      \ /  | /      \  /       | //
  //    $$ |   $$ |$$$$$$ $$$$  |/$$$$$$  |/$$$$$$$/ /$$$$$$  |/$$$$$$  |$$ |/$$$$$$  |/$$$$$$$/  //
  //    $$ |   $$ |$$ | $$ | $$ |$$    $$ |$$      \ $$    $$ |$$ |  $$/ $$ |$$    $$ |$$      \  //
  //    $$ |   $$ |$$ | $$ | $$ |$$$$$$$$/  $$$$$$  |$$$$$$$$/ $$ |      $$ |$$$$$$$$/  $$$$$$  | //
  //    $$ |   $$ |$$ | $$ | $$ |$$       |/     $$/ $$       |$$ |      $$ |$$       |/     $$/  //
  //    $$/    $$/ $$/  $$/  $$/  $$$$$$$/ $$$$$$$/   $$$$$$$/ $$/       $$/  $$$$$$$/ $$$$$$$/   //
  //                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  {
    id: "completeness-timeseries",
    label: {
      en: "Completeness over time",
      fr: "Complétude dans le temps",
    },
    description: {
      en: "Area chart showing completeness trends over time by indicator",
      fr: "Graphique en aires montrant les tendances de complétude dans le temps par indicateur",
    },
    createDefaultVisualizationOnInstall: "26dedd7c-4577-4022-928c-69e0ee790a71",
    allowedFilters: ["indicator_common_id"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "row",
          },
        ],
        filterBy: [],
      },
      s: {
        content: "lines-area",
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Indicator completeness over time",
          fr: "Complétude des indicateurs dans le temps",
        },
        subCaption: {
          en: "Percentage of facility-months with complete data DATE_RANGE",
          fr: "Pourcentage de mois-établissements avec des données complètes DATE_RANGE",
        },
        footnote: {
          en: "Higher completeness improves the reliability of the data, especially when completeness is stable over time. Completeness is defined as the percentage of reporting facilities each month out of the total number of facilities expected to report. A facility is expected to report if it has reported any volume for each indicator anytime within a year. A high completeness does not indicate that the HMIS is representative of all service delivery in the country, as some services may not be delivered in facilities, or some facilities may not report.",
          fr: "Une complétude élevée améliore la fiabilité des données, surtout lorsqu'elle est stable dans le temps. La complétude est définie comme le pourcentage d'établissements déclarants chaque mois par rapport au nombre total d'établissements censés déclarer. Un établissement est censé déclarer s'il a déclaré un volume pour chaque indicateur à tout moment au cours de l'année. Une complétude élevée n'indique pas que le HMIS est représentatif de toute la prestation de services dans le pays, car certains services peuvent ne pas être fournis dans les établissements, ou certains établissements peuvent ne pas déclarer.",
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
  id: "m1-02-02",
  resultsObjectId: "M1_output_completeness.csv",
  valueProps: ["completeness_flag"],
  valueFunc: "AVG",
  valueLabelReplacements: {
    completeness_flag:
      "Binary variable indicating whether the facility meets criteria",
  },
  label: {
    en: "Proportion of completed records",
    fr: "Proportion d'enregistrements complets",
  },
  requiredDisaggregationOptions: [],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Proportion of facility-indicator-period combinations meeting completeness criteria.",
      fr: "Proportion de combinaisons établissement-indicateur-période répondant aux critères de complétude.",
    },
    methodology: {
      en: "AVG of binary completeness_flag. Facilities must report consistently across expected periods to be flagged as complete.",
      fr: "Moyenne du drapeau binaire de complétude. Les établissements doivent déclarer régulièrement pour être considérés comme complets.",
    },
    interpretation: {
      en: "Higher values indicate better reporting consistency. Values below 80% suggest significant reporting gaps that may bias analysis.",
      fr: "Des valeurs plus élevées indiquent une meilleure cohérence des déclarations. Les valeurs inférieures à 80% suggèrent des lacunes significatives.",
    },
    typicalRange: {
      en: "80-100% is good; 60-80% moderate; <60% indicates major gaps.",
      fr: "80-100% est bon; 60-80% modéré; <60% indique des lacunes majeures.",
    },
    caveats: {
      en: "Definition of completeness can vary. Check module parameters for specific criteria used.",
      fr: "La définition de complétude peut varier. Vérifiez les paramètres du module.",
    },
    disaggregationGuidance: {
      en: "Disaggregate by admin_area to identify regions with reporting challenges. Use indicator_common_id to see if specific services have lower compliance.",
      fr: "Désagréger par zone administrative pour identifier les régions avec des défis de déclaration.",
    },
  },
  variantLabel: null,
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
