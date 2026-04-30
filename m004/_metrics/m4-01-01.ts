import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import {
  FOOTNOTE_COVERAGE_BASE,
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
        subCaption: SUBCAPTION_DISCLAIMER,
        footnote: FOOTNOTE_COVERAGE_BASE,
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
