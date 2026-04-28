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
      en: "Coverage timeseries by region",
      fr: "Séries temporelles de couverture par région",
    },
    description: {
      en: "Coverage trends over time by admin area 2",
      fr: "Tendances de couverture dans le temps par zone administrative 2",
    },
    createDefaultVisualizationOnInstall: "a7727717-92d9-4676-b533-9b98be426a81",
    allowedFilters: ["admin_area_2"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "series",
          },
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "cell",
          },
        ],
        filterBy: [],
      },
      s: {
        scale: 1.9,
        content: "lines",
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Coverage estimates",
          fr: "Estimations de couverture",
        },
        subCaption: {
          en: "DATE_RANGE\nDISCLAIMER: These results use routine data to provide rigorous, but not official estimates. They should be interpreted considering any data quality or representation limitations, including data quality findings and any other country specific factors.",
          fr: "DATE_RANGE\nAVERTISSEMENT : Ces résultats utilisent des données de routine pour fournir des estimations rigoureuses mais non officielles. Ils doivent être interprétés en tenant compte des limitations de qualité des données ou de représentativité, y compris les résultats d'évaluation de la qualité des données et tout autre facteur spécifique au pays.",
        },
        footnote: {
          en: "Estimating service coverage from administrative data can provide more timely information on coverage trends, or highlight data quality concerns. Numerators are the volumes reported in HMIS, adjusted for data quality. Denominators are selected from UN projections, survey estimates, or derived from HMIS volume for related indicators. National projections are made by applying HMIS trends to the most recent survey data. Subnational estimates are more sensitive to poor data quality, and projections from surveys are not calculated.\n\nMICS data courtesy of UNICEF. Multiple Indicator Cluster Surveys (various rounds) New York City, New York.\n\nDHS data courtesy of ICF. Demographic and Health Surveys (various rounds). Rockville, Maryland.",
          fr: "L'estimation de la couverture des services à partir des données administratives peut fournir des informations plus rapides sur les tendances de couverture, ou mettre en évidence des problèmes de qualité des données. Les numérateurs sont les volumes déclarés dans le HMIS, ajustés pour la qualité des données. Les dénominateurs sont sélectionnés à partir des projections de l'ONU, des estimations d'enquête, ou dérivés du volume HMIS pour les indicateurs liés. Les projections nationales sont réalisées en appliquant les tendances HMIS aux données d'enquête les plus récentes. Les estimations sous-nationales sont plus sensibles à la mauvaise qualité des données, et les projections à partir des enquêtes ne sont pas calculées.\n\nDonnées MICS avec l'aimable autorisation de l'UNICEF. Enquêtes par grappes à indicateurs multiples (divers cycles). New York, New York.\n\nDonnées EDS avec l'aimable autorisation d'ICF. Enquêtes démographiques et de santé (divers cycles). Rockville, Maryland.",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
  {
    id: "coverage-bar",
    label: {
      en: "Coverage bar chart by region",
      fr: "Diagramme à barres de couverture par région",
    },
    description: {
      en: "Bar chart comparing coverage across regions",
      fr: "Diagramme à barres comparant la couverture entre régions",
    },
    createDefaultVisualizationOnInstall: "d452dfcf-2cc9-4c7f-bfb0-bf5b8ab6433d",
    allowedFilters: ["admin_area_2"],
    config: {
      d: {
        type: "chart",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "indicator",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "replicant",
          },
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "indicator",
          },
          {
            disOpt: "year",
            disDisplayOpt: "cell",
          },
        ],
        filterBy: [],
        selectedReplicantValue: "anc1",
      },
      s: {
        colorScale: "single-grey",
        decimalPlaces: 1,
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
          en: "Estimating service coverage from administrative data can provide more timely information on coverage trends, or highlight data quality concerns. Numerators are the volumes reported in HMIS, adjusted for data quality. Denominators are selected from UN projections, survey estimates, or derived from HMIS volume for related indicators. National projections are made by applying HMIS trends to the most recent survey data. Subnational estimates are more sensitive to poor data quality, and projections from surveys are not calculated.\n\nMICS data courtesy of UNICEF. Multiple Indicator Cluster Surveys (various rounds) New York City, New York.\n\nDHS data courtesy of ICF. Demographic and Health Surveys (various rounds). Rockville, Maryland.",
          fr: "L'estimation de la couverture des services à partir des données administratives peut fournir des informations plus rapides sur les tendances de couverture, ou mettre en évidence des problèmes de qualité des données. Les numérateurs sont les volumes déclarés dans le HMIS, ajustés pour la qualité des données. Les dénominateurs sont sélectionnés à partir des projections de l'ONU, des estimations d'enquête, ou dérivés du volume HMIS pour les indicateurs liés. Les projections nationales sont réalisées en appliquant les tendances HMIS aux données d'enquête les plus récentes. Les estimations sous-nationales sont plus sensibles à la mauvaise qualité des données, et les projections à partir des enquêtes ne sont pas calculées.\n\nDonnées MICS avec l'aimable autorisation de l'UNICEF. Enquêtes par grappes à indicateurs multiples (divers cycles). New York, New York.\n\nDonnées EDS avec l'aimable autorisation d'ICF. Enquêtes démographiques et de santé (divers cycles). Rockville, Maryland.",
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
  id: "m4-02-01",
  resultsObjectId: "M4_coverage_estimation_admin_area_2.csv",
  valueProps: ["coverage_cov"],
  valueFunc: "AVG",
  valueLabelReplacements: {
    coverage_cov: "Coverage calculated from HMIS data",
  },
  label: {
    en: "Coverage calculated from HMIS data",
    fr: "Couverture calculée à partir des données HMIS",
  },
  variantLabel: {
    en: "Admin Area 2",
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
      en: "HMIS-derived health service coverage at admin area 2 (province/state) level.",
      fr: "Couverture des services de santé dérivée du HMIS au niveau de la zone administrative 2 (province/état).",
    },
    methodology: {
      en: "AVG of coverage calculated as HMIS service volumes divided by subnational population denominators. Uses nationally-selected denominator with subnational fallback if national-only denominator chosen.",
      fr: "Moyenne de la couverture calculée comme volumes de services HMIS divisés par dénominateurs de population sous-nationale.",
    },
    interpretation: {
      en: "Enables subnational coverage monitoring and equity analysis. Compare across regions to identify geographic disparities. Coverage >100% indicates denominator or data quality issues.",
      fr: "Permet la surveillance de la couverture sous-nationale et l'analyse de l'équité. Comparer entre régions pour identifier les disparités géographiques.",
    },
    typicalRange: {
      en: "0-100%. Regional variation expected; coverage gaps often larger in remote/underserved areas.",
      fr: "0-100%. Variation régionale attendue; écarts de couverture souvent plus grands dans les zones éloignées.",
    },
    caveats: {
      en: "Subnational denominators may be less reliable than national. Migration and population estimates affect accuracy. Some denominators only available at national level.",
      fr: "Les dénominateurs sous-nationaux peuvent être moins fiables que nationaux. Les estimations de migration et population affectent la précision.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id, admin_area_2, and year (all required). Map visualization effectively shows geographic coverage patterns. Time series reveals regional improvement or deterioration.",
      fr: "Toujours désagréger par indicator_common_id, admin_area_2 et year (tous requis). La visualisation cartographique montre efficacement les modèles de couverture géographique.",
    },
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
