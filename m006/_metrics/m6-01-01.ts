import type { MetricDefinitionGithub, VizPreset } from "../../.validation/_module_definition_github.ts";

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
    coverage_avgsurveyprojection: "Projected survey estimate (when survey data is missing)",
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
