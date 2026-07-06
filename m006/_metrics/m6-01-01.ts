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
      pt: "Séries temporais de cobertura (nacional)",
    },
    description: {
      en: "National coverage timeseries with survey benchmarks",
      fr: "Séries temporelles de couverture nationale avec repères d'enquête",
      pt: "Séries temporais de cobertura nacional com referências de inquéritos",
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
        content: "lines",
        showDataLabels: true,
        specialCoverageChart: true,
      },
      t: {
        caption: {
          en: "Coverage estimates for REPLICANT",
          fr: "Estimations de couverture pour REPLICANT",
          pt: "Estimativas de cobertura para REPLICANT",
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
    pt: "Cobertura (todos os tipos de estimativa)",
  },
  variantLabel: {
    en: "National",
    fr: "National",
    pt: "Nacional",
  },
  requiredDisaggregationOptions: ["indicator_common_id", "year"],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Comprehensive coverage estimates at national level using user-selected denominators, combining survey benchmarks with HMIS-derived trends.",
      fr: "Estimations de couverture complètes au niveau national utilisant des dénominateurs sélectionnés par l'utilisateur, combinant repères d'enquête avec tendances dérivées du HMIS.",
      pt: "Estimativas de cobertura abrangentes ao nível nacional utilizando denominadores selecionados pelo utilizador, combinando referências de inquéritos com tendências derivadas do SIS.",
    },
    methodology: {
      en: "AVG of three coverage types: (1) original survey estimates when available, (2) survey estimates projected forward using HMIS year-over-year deltas (additive method), (3) HMIS-derived coverage using user-specified denominators. Denominators selected via module parameters (e.g., DENOM_ANC1, DENOM_PENTA1).",
      fr: "Moyenne de trois types de couverture: (1) estimations d'enquête originales, (2) estimations d'enquête projetées en utilisant les deltas HMIS année par année, (3) couverture dérivée du HMIS avec dénominateurs spécifiés.",
      pt: "Média de três tipos de cobertura: (1) estimativas de inquéritos originais quando disponíveis, (2) estimativas de inquéritos projetadas utilizando os deltas anuais do SIS (método aditivo), (3) cobertura derivada do SIS utilizando denominadores especificados pelo utilizador. Os denominadores são selecionados através dos parâmetros do módulo (por exemplo, DENOM_ANC1, DENOM_PENTA1).",
    },
    interpretation: {
      en: "This is the final, policy-relevant coverage metric combining best available data sources. Survey estimates anchor coverage to gold-standard benchmarks; projected estimates fill inter-survey gaps using HMIS momentum; HMIS coverage enables annual monitoring. Concordance between the three types validates data quality and denominator selection.",
      fr: "C'est la métrique de couverture finale pertinente pour les politiques, combinant les meilleures sources de données disponibles. Les estimations d'enquête ancrent la couverture aux repères étalon-or.",
      pt: "Esta é a métrica de cobertura final, relevante para as políticas, que combina as melhores fontes de dados disponíveis. As estimativas de inquéritos ancoram a cobertura a referências padrão-ouro; as estimativas projetadas preenchem os intervalos entre inquéritos utilizando o ritmo do SIS; a cobertura do SIS permite a monitorização anual. A concordância entre os três tipos valida a qualidade dos dados e a seleção do denominador.",
    },
    typicalRange: {
      en: "0-100%. Coverage >100% indicates denominator or data quality problems requiring investigation.",
      fr: "0-100%. Couverture >100% indique des problèmes de dénominateur ou de qualité des données nécessitant investigation.",
      pt: "0-100%. Uma cobertura >100% indica problemas de denominador ou de qualidade dos dados que requerem investigação.",
    },
    caveats: {
      en: "Projection method assumes HMIS trends accurately reflect true coverage changes. Denominator selection (via module parameters) critically affects results - inappropriate denominators produce implausible coverage. Survey timing and data quality affect baseline accuracy.",
      fr: "La méthode de projection suppose que les tendances HMIS reflètent avec précision les vrais changements de couverture. La sélection du dénominateur affecte de manière critique les résultats.",
      pt: "O método de projeção pressupõe que as tendências do SIS refletem com precisão as verdadeiras alterações da cobertura. A seleção do denominador (através dos parâmetros do módulo) afeta de forma crítica os resultados - denominadores inadequados produzem coberturas implausíveis. O momento dos inquéritos e a qualidade dos dados afetam a exatidão da linha de base.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id and year (both required). Compare the three coverage types to assess data quality - they should show consistent trends if both HMIS and denominators are reliable. Time series reveals coverage evolution and inter-survey projection accuracy.",
      fr: "Toujours désagréger par indicator_common_id et year (tous deux requis). Comparer les trois types de couverture pour évaluer la qualité des données.",
      pt: "Desagregar sempre por indicator_common_id e year (ambos obrigatórios). Comparar os três tipos de cobertura para avaliar a qualidade dos dados - devem apresentar tendências consistentes se tanto o SIS como os denominadores forem fiáveis. As séries temporais revelam a evolução da cobertura e a exatidão da projeção entre inquéritos.",
    },
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
