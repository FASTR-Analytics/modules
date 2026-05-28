import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import {
  FOOTNOTE_COVERAGE_WITH_CURRENT_YEAR,
  SUBCAPTION_DISCLAIMER,
} from "../../_shared/text_presets.ts";

export const vizPresets: VizPreset[] = [
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
    createDefaultVisualizationOnInstall: "2a74f737-78e5-41a1-8f6d-7a3f59be2d19",
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
        scale: 2.7,
        content: "lines",
        showDataLabels: true,
        specialCoverageChart: true,
      },
      t: {
        caption: {
          en: "Coverage estimates for REPLICANT",
          fr: "Estimations de couverture pour REPLICANT",
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
];

export const metric: MetricDefinitionGithub = {
  id: "m6-01-01",
  resultsObjectId: "M6_coverage_estimation_national.csv",
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
    en: "National",
    fr: "National",
  },
  requiredDisaggregationOptions: ["indicator_common_id", "year"],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Comprehensive coverage estimates at national level using user-selected denominators, combining survey benchmarks with HMIS-derived trends.",
      fr: "Estimations de couverture complètes au niveau national utilisant des dénominateurs sélectionnés par l'utilisateur, combinant repères d'enquête avec tendances dérivées du HMIS.",
    },
    methodology: {
      en: "AVG of three coverage types: (1) original survey estimates when available, (2) survey estimates projected forward using HMIS year-over-year deltas (additive method), (3) HMIS-derived coverage using user-specified denominators. Denominators selected via module parameters (e.g., DENOM_ANC1, DENOM_PENTA1).",
      fr: "Moyenne de trois types de couverture: (1) estimations d'enquête originales, (2) estimations d'enquête projetées en utilisant les deltas HMIS année par année, (3) couverture dérivée du HMIS avec dénominateurs spécifiés.",
    },
    interpretation: {
      en: "This is the final, policy-relevant coverage metric combining best available data sources. Survey estimates anchor coverage to gold-standard benchmarks; projected estimates fill inter-survey gaps using HMIS momentum; HMIS coverage enables annual monitoring. Concordance between the three types validates data quality and denominator selection.",
      fr: "C'est la métrique de couverture finale pertinente pour les politiques, combinant les meilleures sources de données disponibles. Les estimations d'enquête ancrent la couverture aux repères étalon-or.",
    },
    typicalRange: {
      en: "0-100%. Coverage >100% indicates denominator or data quality problems requiring investigation.",
      fr: "0-100%. Couverture >100% indique des problèmes de dénominateur ou de qualité des données nécessitant investigation.",
    },
    caveats: {
      en: "Projection method assumes HMIS trends accurately reflect true coverage changes. Denominator selection (via module parameters) critically affects results - inappropriate denominators produce implausible coverage. Survey timing and data quality affect baseline accuracy.",
      fr: "La méthode de projection suppose que les tendances HMIS reflètent avec précision les vrais changements de couverture. La sélection du dénominateur affecte de manière critique les résultats.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id and year (both required). Compare the three coverage types to assess data quality - they should show consistent trends if both HMIS and denominators are reliable. Time series reveals coverage evolution and inter-survey projection accuracy.",
      fr: "Toujours désagréger par indicator_common_id et year (tous deux requis). Comparer les trois types de couverture pour évaluer la qualité des données.",
    },
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
