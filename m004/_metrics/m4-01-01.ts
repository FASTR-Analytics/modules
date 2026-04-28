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
      en: "Coverage timeseries (national)",
      fr: "Séries temporelles de couverture (national)",
    },
    description: {
      en: "National coverage timeseries with survey benchmarks",
      fr: "Séries temporelles de couverture nationale avec repères d'enquête",
    },
    createDefaultVisualizationOnInstall: "3e3230cb-ad9e-48b9-b3ce-7bd01255d20b",
    allowedFilters: [],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [
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
        specialCoverageChart: true,
      },
      t: {
        caption: {
          en: "Coverage estimates for REPLICANT",
          fr: "Estimations de couverture pour REPLICANT",
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
  id: "m4-01-01",
  resultsObjectId: "M4_coverage_estimation.csv",
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
    en: "Coverage calculated from HMIS data",
    fr: "Couverture calculée à partir des données HMIS",
  },
  variantLabel: {
    en: "National",
    fr: "National",
  },
  requiredDisaggregationOptions: ["indicator_common_id", "year"],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Health service coverage estimates at national level, comparing HMIS-derived coverage with survey-based benchmarks.",
      fr: "Estimations de couverture des services de santé au niveau national, comparant la couverture dérivée du HMIS avec les repères d'enquête.",
    },
    methodology: {
      en: "AVG of three coverage types: (1) original survey estimates when available, (2) projected survey estimates using HMIS trends, (3) HMIS-derived coverage calculated as service volumes divided by population denominators. Denominators selected based on minimizing error against survey benchmarks.",
      fr: "Moyenne de trois types de couverture: (1) estimations d'enquête originales, (2) estimations d'enquête projetées utilisant les tendances HMIS, (3) couverture dérivée du HMIS.",
    },
    interpretation: {
      en: "Three values provide complementary perspectives: survey estimates are gold standard but sparse; projected estimates fill gaps using HMIS trends; HMIS-derived estimates enable annual monitoring. Large gaps between HMIS and survey coverage suggest data quality issues or denominator problems.",
      fr: "Trois valeurs fournissent des perspectives complémentaires: les estimations d'enquête sont l'étalon-or mais rares; les estimations projetées comblent les lacunes.",
    },
    typicalRange: {
      en: "0-100% for coverage. Maternal services typically 40-80%; vaccination 60-95%; varies by country context.",
      fr: "0-100% pour la couverture. Services maternels généralement 40-80%; vaccination 60-95%; varie selon le contexte.",
    },
    caveats: {
      en: "Denominator selection is critical - inappropriate denominators can produce implausible coverage >100%. Projection assumes HMIS trends reflect true coverage changes. Survey timing and HMIS data quality affect comparability.",
      fr: "La sélection du dénominateur est critique - les dénominateurs inappropriés peuvent produire une couverture >100%. La projection suppose que les tendances HMIS reflètent les vrais changements.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id and year (both required). Compare the three coverage types to assess HMIS-survey concordance. Time series reveals coverage trends and data quality evolution.",
      fr: "Toujours désagréger par indicator_common_id et year (tous deux requis). Comparer les trois types de couverture pour évaluer la concordance HMIS-enquête.",
    },
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
