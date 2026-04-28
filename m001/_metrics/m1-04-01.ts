import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { CF_80_70 } from "../../.validation/cf_presets.ts";

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
    id: "dqa-score-table",
    label: {
      en: "Overall DQA score table",
      fr: "Tableau du score EQD global",
    },
    description: {
      en: "Table showing DQA scores by region and year",
      fr: "Tableau montrant les scores EQD par région et année",
    },
    createDefaultVisualizationOnInstall: "MAKE THIS FOR ME",
    allowedFilters: ["admin_area_2"],
    config: {
      d: {
        type: "table",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "col",
        disaggregateBy: [
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "row",
          },
          {
            disOpt: "year",
            disDisplayOpt: "col",
          },
        ],
        filterBy: [],
      },
      s: {
        content: "lines",
        ...CF_80_70,
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Overall DQA score",
          fr: "Score EQD global",
        },
        subCaption: {
          en: "Percentage of facility-months with adequate data quality over time",
          fr: "Pourcentage de mois-établissements avec une qualité des données adéquate dans le temps",
        },
        footnote: {
          en: "Adequate data quality is defined as: 1) No missing data or outliers for OPD, Penta1, and ANC1, where available 2) Consistent reporting between Penta1/Penta3 and ANC1/ANC4.",
          fr: "La qualité adéquate des données est définie comme : 1) Pas de données manquantes ou de valeurs aberrantes pour OPD, Penta1 et ANC1, lorsque disponibles 2) Déclaration cohérente entre Penta1/Penta3 et ANC1/ANC4.",
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
    id: "dqa-score-map",
    label: {
      en: "Overall DQA score map",
      fr: "Carte du score EQD global",
    },
    description: {
      en: "Map showing DQA scores by region",
      fr: "Carte montrant les scores EQD par région",
    },
    createDefaultVisualizationOnInstall: "7f08fbde-48e6-44d3-b1db-8771d96120c6",
    allowedFilters: [],
    config: {
      d: {
        type: "map",
        valuesDisDisplayOpt: "mapArea",
        disaggregateBy: [
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "mapArea",
          },
        ],
        filterBy: [],
        periodFilter: {
          filterType: "last_n_months",
          nMonths: 12,
        },
      },
      s: {
        ...CF_80_70,
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Overall DQA score",
          fr: "Score EQD global",
        },
        subCaption: {
          en: "Percentage of facility-months with adequate data quality, DATE_RANGE",
          fr: "Pourcentage de mois-établissements avec une qualité des données adéquate, DATE_RANGE",
        },
        footnote: {
          en: "Adequate data quality is defined as: 1) No missing data or outliers for OPD, Penta1, and ANC1, where available 2) Consistent reporting between Penta1/Penta3 and ANC1/ANC4.",
          fr: "La qualité adéquate des données est définie comme : 1) Pas de données manquantes ou de valeurs aberrantes pour OPD, Penta1 et ANC1, lorsque disponibles 2) Déclaration cohérente entre Penta1/Penta3 et ANC1/ANC4.",
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
  id: "m1-04-01",
  resultsObjectId: "M1_output_dqa.csv",
  valueProps: ["dqa_score"],
  valueFunc: "AVG",
  valueLabelReplacements: {
    dqa_score: "Binary variable indicating adequate data quality",
  },
  label: {
    en: "Proportion of facilities with adequate data quality",
    fr: "Proportion d'établissements avec une qualité des données adéquate",
  },
  requiredDisaggregationOptions: [],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Proportion of facilities meeting the composite data quality assessment threshold.",
      fr: "Proportion d'établissements atteignant le seuil d'évaluation composite de la qualité des données.",
    },
    methodology: {
      en: "AVG of binary dqa_score based on composite assessment of completeness, outliers, and consistency.",
      fr: "Moyenne du score DQA binaire basé sur une évaluation composite de la complétude, des valeurs aberrantes et de la cohérence.",
    },
    interpretation: {
      en: "Higher values indicate more facilities with trustworthy data. Facilities below threshold may need additional support or data verification.",
      fr: "Des valeurs plus élevées indiquent plus d'établissements avec des données fiables.",
    },
    typicalRange: {
      en: "70-100% is acceptable; <70% indicates widespread quality issues.",
      fr: "70-100% est acceptable; <70% indique des problèmes de qualité généralisés.",
    },
    caveats: {
      en: "Composite score weights may need adjustment based on local context. Consider individual components for detailed diagnosis.",
      fr: "Les poids du score composite peuvent nécessiter un ajustement selon le contexte local.",
    },
    disaggregationGuidance: {
      en: "Disaggregate by admin_area to identify regions needing quality improvement support. Use facility_type to see if certain facility levels have more challenges.",
      fr: "Désagréger par zone administrative pour identifier les régions nécessitant un soutien.",
    },
  },
  variantLabel: null,
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
