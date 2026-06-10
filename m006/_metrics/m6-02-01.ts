import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import {
  FOOTNOTE_COVERAGE_WITH_CURRENT_YEAR,
  SUBCAPTION_DISCLAIMER,
} from "../../_shared/text_presets.ts";

export const vizPresets: VizPreset[] = [
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
    id: "coverage-timeseries",
    label: {
      en: "Coverage timeseries by Admin Area 2",
      fr: "Séries temporelles de couverture par Zone administrative 2",
    },
    description: {
      en: "Coverage trends over time by admin area 2",
      fr: "Tendances de couverture dans le temps par zone administrative 2",
    },
    createDefaultVisualizationOnInstall: "e5f8740b-a690-4a84-a0cd-05d529676f26",
    allowedFilters: ["admin_area_2"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "cell",
          },
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "replicant",
          },
        ],
        filterBy: [],
        selectedReplicantValue: "anc1",
      },
      s: {
        content: "lines",
        showDataLabels: true,
        specialCoverageChart: true,
      },
      t: {
        caption: {
          en: "Subnational coverage estimates for REPLICANT",
          fr: "Estimations de couverture sous-nationales pour REPLICANT",
        },
        subCaption: SUBCAPTION_DISCLAIMER,
        footnote: FOOTNOTE_COVERAGE_WITH_CURRENT_YEAR,
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
  ////////////////////////////////////
  //  _______                       //
  // /       \                      //
  // $$$$$$$  |  ______    ______   //
  // $$ |__$$ | /      \  /      \  //
  // $$    $$<  $$$$$$  |/$$$$$$  | //
  // $$$$$$$  | /    $$ |$$ |  $$/  //
  // $$ |__$$ |/$$$$$$$ |$$ |       //
  // $$    $$/ $$    $$ |$$ |       //
  // $$$$$$$/   $$$$$$$/ $$/        //
  //                                //
  ////////////////////////////////////
  {
    id: "coverage-bar",
    label: {
      en: "Coverage bar chart by Admin Area 2",
      fr: "Diagramme à barres de couverture par Zone administrative 2",
    },
    description: {
      en: "Bar chart comparing coverage across Admin Areas 2",
      fr: "Diagramme à barres comparant la couverture entre Zones administratives 2",
    },
    createDefaultVisualizationOnInstall: "9d4977b4-0d87-44e1-b2bd-3eddcba623f4",
    allowedFilters: ["admin_area_2"],
    config: {
      d: {
        type: "chart",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "indicator",
          },
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "replicant",
          },
        ],
        filterBy: [],
        periodFilter: {
          filterType: "last_n_months",
          nMonths: 12,
        },
        selectedReplicantValue: "anc4",
        valuesFilter: ["coverage_cov"],
      },
      s: {
        showDataLabels: true,
        colorScale: "single-grey",
        sortIndicatorValues: "descending",
      },
      t: {
        caption: {
          en: "Sub-national level coverage estimates, REPLICANT",
          fr: "Estimations de couverture au niveau sous-national, REPLICANT",
        },
        subCaption: {
          en: "DATE_RANGE",
          fr: "PLAGE_DE_DATES",
        },
        footnote: FOOTNOTE_COVERAGE_WITH_CURRENT_YEAR,
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
  id: "m6-02-01",
  resultsObjectId: "M6_coverage_estimation_admin2.csv",
  valueProps: [
    "coverage_original_estimate",
    "coverage_avgsurveyprojection",
    "coverage_cov",
  ],
  valueFunc: "AVG",
  valueLabelReplacements: {
    coverage_original_estimate: "Survey-based estimate (when available)",
    coverage_avgsurveyprojection:
      "Projected survey estimate (when survey data is missing)",
    coverage_cov: "Coverage calculated from HMIS data",
  },
  label: {
    en: "Coverage (all estimation types)",
    fr: "Couverture (tous types d'estimation)",
  },
  variantLabel: {
    en: "Admin area 2",
    fr: "Zone administrative 2",
  },
  requiredDisaggregationOptions: [
    "indicator_common_id",
    "admin_area_2",
    "year",
  ],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Subnational coverage estimates at admin area 2 level combining survey, projected, and HMIS-derived values.",
      fr: "Estimations de couverture sous-nationales au niveau de la zone administrative 2 combinant valeurs d'enquête, projetées et dérivées du HMIS.",
    },
    methodology: {
      en: "AVG of survey, projected survey, and HMIS coverage at regional level. Uses user-specified denominators and survey projection methodology.",
      fr: "Moyenne de la couverture d'enquête, d'enquête projetée et HMIS au niveau régional. Utilise des dénominateurs spécifiés par l'utilisateur.",
    },
    interpretation: {
      en: "Enables regional equity analysis and geographic targeting. Compare across admin areas to identify coverage disparities. Three estimation types validate subnational data quality.",
      fr: "Permet l'analyse de l'équité régionale et le ciblage géographique. Comparer entre zones administratives pour identifier les disparités de couverture.",
    },
    typicalRange: {
      en: "0-100%. Regional variation expected; remote areas typically have lower coverage.",
      fr: "0-100%. Variation régionale attendue; zones éloignées ont généralement une couverture plus faible.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id, admin_area_2, and year (all required). Compare three coverage types for regional data quality validation.",
      fr: "Toujours désagréger par indicator_common_id, admin_area_2 et year (tous requis). Comparer trois types de couverture pour la validation de la qualité des données régionales.",
    },
    caveats: null,
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
