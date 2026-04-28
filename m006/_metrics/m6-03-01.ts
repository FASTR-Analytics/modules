import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

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
      en: "Coverage timeseries by district",
      fr: "Séries temporelles de couverture par district",
    },
    description: {
      en: "Coverage trends over time by admin area 3",
      fr: "Tendances de couverture dans le temps par zone administrative 3",
    },
    createDefaultVisualizationOnInstall: "e5f8740b-a690-4a84-a0cd-05d529676f27",
    allowedFilters: ["admin_area_3"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [
          {
            disOpt: "admin_area_3",
            disDisplayOpt: "cell",
          },
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "replicant",
          },
        ],
        filterBy: [],
        includeNationalForAdminArea2: true,
        selectedReplicantValue: "anc1",
      },
      s: {
        scale: 1.7,
        content: "lines",
        showDataLabels: true,
        specialCoverageChart: true,
      },
      t: {
        caption: {
          en: "Subnational coverage estimates for REPLICANT",
          fr: "Estimations de couverture sous-nationales pour REPLICANT",
        },
        subCaption: {
          en: "DATE_RANGE\nDISCLAIMER: These results use routine data to provide rigorous, but not official estimates. They should be interpreted considering any data quality or representation limitations, including data quality findings and any other country specific factors.",
          fr: "DATE_RANGE\nAVERTISSEMENT : Ces résultats utilisent des données de routine pour fournir des estimations rigoureuses mais non officielles. Ils doivent être interprétés en tenant compte des limitations de qualité des données ou de représentativité, y compris les résultats d'évaluation de la qualité des données et tout autre facteur spécifique au pays.",
        },
        footnote: {
          en: "Estimating service coverage from administrative data can provide more timely information on coverage trends, or highlight data quality concerns. Numerators are the volumes reported in HMIS, adjusted for data quality. Denominators are selected from UN projections, survey estimates, or derived from HMIS volume for related indicators. National projections are made by applying HMIS trends to the most recent survey data. Subnational estimates are more sensitive to poor data quality, and projections from surveys are not calculated.\n\nMICS data courtesy of UNICEF. Multiple Indicator Cluster Surveys (various rounds) New York City, New York.\n\nDHS data courtesy of ICF. Demographic and Health Surveys (various rounds). Rockville, Maryland.\n\nData for the current year reflect the period available at the time of analysis; population figures have been scaled to match the corresponding duration.",
          fr: "L'estimation de la couverture des services à partir des données administratives peut fournir des informations plus rapides sur les tendances de couverture, ou mettre en évidence des problèmes de qualité des données. Les numérateurs sont les volumes déclarés dans le HMIS, ajustés pour la qualité des données. Les dénominateurs sont sélectionnés à partir des projections de l'ONU, des estimations d'enquête, ou dérivés du volume HMIS pour les indicateurs liés. Les projections nationales sont réalisées en appliquant les tendances HMIS aux données d'enquête les plus récentes. Les estimations sous-nationales sont plus sensibles à la mauvaise qualité des données, et les projections à partir des enquêtes ne sont pas calculées.\n\nDonnées MICS avec l'aimable autorisation de l'UNICEF. Enquêtes par grappes à indicateurs multiples (divers cycles). New York, New York.\n\nDonnées EDS avec l'aimable autorisation d'ICF. Enquêtes démographiques et de santé (divers cycles). Rockville, Maryland.\n\nLes données de l'année en cours reflètent la période disponible au moment de l'analyse ; les chiffres de population ont été ajustés pour correspondre à la durée correspondante.",
        },
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
      en: "Coverage bar chart by district",
      fr: "Diagramme à barres de couverture par district",
    },
    description: {
      en: "Bar chart comparing coverage across districts",
      fr: "Diagramme à barres comparant la couverture entre districts",
    },
    createDefaultVisualizationOnInstall: "9d4977b4-0d87-44e1-b2bd-3eddcba623f5",
    allowedFilters: ["admin_area_3"],
    config: {
      d: {
        type: "chart",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [
          {
            disOpt: "admin_area_3",
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
          fr: "DATE_RANGE",
        },
        footnote: {
          en: "Estimating service coverage from administrative data can provide more timely information on coverage trends, or highlight data quality concerns. Numerators are the volumes reported in HMIS, adjusted for data quality. Denominators are selected from UN projections, survey estimates, or derived from HMIS volume for related indicators. National projections are made by applying HMIS trends to the most recent survey data. Subnational estimates are more sensitive to poor data quality, and projections from surveys are not calculated.\n\nMICS data courtesy of UNICEF. Multiple Indicator Cluster Surveys (various rounds) New York City, New York.\n\nDHS data courtesy of ICF. Demographic and Health Surveys (various rounds). Rockville, Maryland.\n\nData for the current year reflect the period available at the time of analysis; population figures have been scaled to match the corresponding duration.",
          fr: "L'estimation de la couverture des services à partir des données administratives peut fournir des informations plus rapides sur les tendances de couverture, ou mettre en évidence des problèmes de qualité des données. Les numérateurs sont les volumes déclarés dans le HMIS, ajustés pour la qualité des données. Les dénominateurs sont sélectionnés à partir des projections de l'ONU, des estimations d'enquête, ou dérivés du volume HMIS pour les indicateurs liés. Les projections nationales sont réalisées en appliquant les tendances HMIS aux données d'enquête les plus récentes. Les estimations sous-nationales sont plus sensibles à la mauvaise qualité des données, et les projections à partir des enquêtes ne sont pas calculées.\n\nDonnées MICS avec l'aimable autorisation de l'UNICEF. Enquêtes par grappes à indicateurs multiples (divers cycles). New York, New York.\n\nDonnées EDS avec l'aimable autorisation d'ICF. Enquêtes démographiques et de santé (divers cycles). Rockville, Maryland.\n\nLes données de l'année en cours reflètent la période disponible au moment de l'analyse ; les chiffres de population ont été ajustés pour correspondre à la durée correspondante.",
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
  id: "m6-03-01",
  resultsObjectId: "M6_coverage_estimation_admin3.csv",
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
    en: "Admin area 3",
    fr: "Zone administrative 3",
  },
  requiredDisaggregationOptions: [
    "indicator_common_id",
    "admin_area_3",
    "year",
  ],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "District-level coverage estimates combining survey, projected, and HMIS-derived values.",
      fr: "Estimations de couverture au niveau du district combinant valeurs d'enquête, projetées et dérivées du HMIS.",
    },
    methodology: {
      en: "AVG of survey, projected survey, and HMIS coverage at district level. Finest geographic resolution for coverage estimation.",
      fr: "Moyenne de la couverture d'enquête, d'enquête projetée et HMIS au niveau du district. Résolution géographique la plus fine.",
    },
    interpretation: {
      en: "Enables district-level targeting and micro-planning. Interpret with caution due to smaller sample sizes and denominator uncertainty.",
      fr: "Permet le ciblage au niveau du district et la micro-planification. Interpréter avec prudence en raison de petites tailles d'échantillon.",
    },
    typicalRange: {
      en: "0-100%. Greater variation expected at district level; interpret with caution.",
      fr: "0-100%. Plus grande variation attendue au niveau du district; interpréter avec prudence.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id, admin_area_3, and year (all required). Compare with admin area 2 for context.",
      fr: "Toujours désagréger par indicator_common_id, admin_area_3 et year (tous requis). Comparer avec la zone administrative 2 pour le contexte.",
    },
    caveats: null,
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
