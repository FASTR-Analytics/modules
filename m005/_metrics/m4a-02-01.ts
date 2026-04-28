import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "coverage-timeseries",
    label: {
      en: "Coverage by denominator type",
      fr: "Couverture par type de dénominateur",
    },
    description: {
      en: "Timeseries comparing coverage across denominator sources",
      fr: "Séries temporelles comparant la couverture entre sources de dénominateur",
    },
    createDefaultVisualizationOnInstall: "15ca88bd-6183-4e71-bb26-3277dd8eb02f",
    allowedFilters: ["denominator_best_or_survey"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [
          {
            disOpt: "denominator_best_or_survey",
            disDisplayOpt: "series",
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
      },
      t: {
        caption: {
          en: "Coverage based of different denominators, REPLICANT",
          fr: "Couverture selon différents dénominateurs, REPLICANT",
        },
        subCaption: {
          en: "DATE_RANGE",
          fr: "DATE_RANGE",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnote: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
];

export const metric: MetricDefinitionGithub = {
  id: "m4a-02-01",
  resultsObjectId: "M5_combined_results_national.csv",
  valueProps: ["value"],
  valueFunc: "AVG",
  valueLabelReplacements: {
    best: "Best",
    survey: "Survey",
    danc1_birth:
      "Estimated number of total births (live + stillbirths) derived from HMIS data on ANC 1st visits.",
    danc1_delivery:
      "Estimated number of deliveries derived from HMIS data on ANC 1st visits.",
    danc1_dpt:
      "Estimated number of infants eligible for DPT1 derived from HMIS data on ANC 1st visits.",
    danc1_livebirth:
      "Estimated number of live births derived from HMIS data on ANC 1st visits.",
    danc1_measles1:
      "Estimated number of children eligible for measles dose 1 (MCV1) derived from HMIS data on ANC 1st visits.",
    danc1_measles2:
      "Estimated number of children eligible for measles dose 2 (MCV2) derived from HMIS data on ANC 1st visits.",
    danc1_pregnancy:
      "Estimated number of pregnancies derived from HMIS data on ANC 1st visits.",
    dbcg_dpt:
      "Estimated number of infants eligible for DPT1 derived from HMIS data on BCG doses.",
    dbcg_livebirth:
      "Estimated number of live births derived from HMIS data on BCG doses.",
    dbcg_pregnancy:
      "Estimated number of pregnancies derived from HMIS data on BCG doses.",
    ddelivery_birth:
      "Estimated number of total births (live + stillbirths) derived from HMIS data on institutional deliveries.",
    ddelivery_dpt:
      "Estimated number of infants eligible for DPT1 derived from HMIS data on institutional deliveries.",
    ddelivery_livebirth:
      "Estimated number of live births derived from HMIS data on institutional deliveries.",
    ddelivery_measles1:
      "Estimated number of children eligible for measles dose 1 (MCV1) derived from HMIS data on institutional deliveries.",
    ddelivery_measles2:
      "Estimated number of children eligible for measles dose 2 (MCV2) derived from HMIS data on institutional deliveries.",
    ddelivery_pregnancy:
      "Estimated number of pregnancies derived from HMIS data on institutional deliveries.",
    dlivebirths_birth:
      "Estimated number of total births (live + stillbirths) derived from HMIS data on live births.",
    dlivebirths_delivery:
      "Estimated number of deliveries derived from HMIS data on live births.",
    dlivebirths_dpt:
      "Estimated number of infants eligible for DPT1 derived from HMIS data on live births.",
    dlivebirths_livebirth:
      "Estimated number of live births derived from HMIS data on live births.",
    dlivebirths_measles1:
      "Estimated number of children eligible for measles dose 1 (MCV1) derived from HMIS data on live births.",
    dlivebirths_measles2:
      "Estimated number of children eligible for measles dose 2 (MCV2) derived from HMIS data on live births.",
    dlivebirths_pregnancy:
      "Estimated number of pregnancies derived from HMIS data on live births.",
    dpenta1_dpt:
      "Estimated number of infants eligible for DPT1 derived from HMIS data on Penta-1 doses.",
    dpenta1_measles1:
      "Estimated number of children eligible for measles dose 1 (MCV1) derived from HMIS data on Penta-1 doses.",
    dpenta1_measles2:
      "Estimated number of children eligible for measles dose 2 (MCV2) derived from HMIS data on Penta-1 doses.",
    dwpp_dpt:
      "Estimated number of infants eligible for DPT1 based on UN WPP estimates.",
    dwpp_livebirth:
      "Estimated number of live births based on UN WPP estimates.",
    dwpp_measles1:
      "Estimated number of children eligible for measles dose 1 (MCV1) based on UN WPP estimates.",
    dwpp_measles2:
      "Estimated number of children eligible for measles dose 2 (MCV2) based on UN WPP estimates.",
    dwpp_pregnancy:
      "Estimated number of pregnancies based on UN WPP estimates.",
    source_anc1: "HMIS data on ANC 1st visits",
    source_delivery: "HMIS data on institutional deliveries",
    source_bcg: "HMIS data on BCG doses",
    source_penta1: "HMIS data on Penta-1 doses",
    source_wpp: "UN WPP estimates",
    source_livebirths: "HMIS data on live births",
    target_pregnancy: "Pregnancies",
    target_delivery: "Deliveries",
    target_birth: "Total births (live + stillbirths)",
    target_livebirth: "Live births",
    target_dpt: "Infants eligible for DPT1",
    target_measles1: "Children eligible for measles dose 1 (MCV1)",
    target_measles2: "Children eligible for measles dose 2 (MCV2)",
    target_vitaminA:
      "Children aged 6-59 months eligible for Vitamin A supplementation",
    target_fully_immunized:
      "Children under 1 year eligible for full immunization",
  },
  label: {
    en: "Coverage estimated with different denominators",
    fr: "Couverture estimée avec différents dénominateurs",
  },
  variantLabel: {
    en: "National",
    fr: "National",
  },
  requiredDisaggregationOptions: [
    "denominator_best_or_survey",
    "indicator_common_id",
    "year",
  ],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Coverage estimates calculated using alternative denominator sources, plus survey benchmarks at national level.",
      fr: "Estimations de couverture calculées utilisant des sources de dénominateur alternatives, plus repères d'enquête au niveau national.",
    },
    methodology: {
      en: "AVG of coverage calculated as HMIS numerators divided by each available denominator type. Includes 'best' denominator (selected to minimize survey error) and 'survey' (original survey estimate). Enables comparison of how denominator choice affects coverage estimates.",
      fr: "Moyenne de la couverture calculée comme numérateurs HMIS divisés par chaque type de dénominateur disponible. Inclut le 'meilleur' dénominateur et l'estimé d'enquête.",
    },
    interpretation: {
      en: "Large variation across denominator types indicates denominator uncertainty. Coverage estimates should be similar to survey values when using appropriate denominators. Coverage >100% suggests denominator underestimation or HMIS over-reporting.",
      fr: "Une grande variation entre types de dénominateur indique une incertitude. Les estimations de couverture devraient être similaires aux valeurs d'enquête avec des dénominateurs appropriés.",
    },
    typicalRange: {
      en: "0-100%. Denominators producing >100% coverage are likely inappropriate or require data quality investigation.",
      fr: "0-100%. Les dénominateurs produisant >100% de couverture sont probablement inappropriés ou nécessitent une investigation de qualité.",
    },
    caveats: {
      en: "Denominator selection is subjective and affects results. The 'best' denominator minimizes squared error against surveys but may not be appropriate for all analytical purposes.",
      fr: "La sélection du dénominateur est subjective et affecte les résultats. Le 'meilleur' dénominateur minimise l'erreur quadratique mais peut ne pas être approprié pour tous les objectifs.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by denominator_best_or_survey, indicator_common_id, and year (all required). Compare 'best' vs 'survey' to validate denominator selection. Visualize all denominator types to assess range of plausible coverage estimates.",
      fr: "Toujours désagréger par denominator_best_or_survey, indicator_common_id et year (tous requis). Comparer 'meilleur' vs 'enquête' pour valider la sélection du dénominateur.",
    },
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
