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
      en: "Coverage timeseries by Admin Area 3",
      fr: "Séries temporelles de couverture par Zone administrative 3",
      pt: "Séries temporais de cobertura por Área administrativa 3",
    },
    description: {
      en: "Coverage trends over time by admin area 3",
      fr: "Tendances de couverture dans le temps par zone administrative 3",
      pt: "Tendências de cobertura ao longo do tempo por área administrativa 3",
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
          pt: "Estimativas de cobertura subnacionais para REPLICANT",
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
      en: "Coverage bar chart by Admin Area 3",
      fr: "Diagramme à barres de couverture par Zone administrative 3",
      pt: "Gráfico de barras de cobertura por Área administrativa 3",
    },
    description: {
      en: "Bar chart comparing coverage across Admin Areas 3",
      fr: "Diagramme à barres comparant la couverture entre Zones administratives 3",
      pt: "Gráfico de barras que compara a cobertura entre Áreas administrativas 3",
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
          pt: "Estimativas de cobertura ao nível subnacional, REPLICANT",
        },
        subCaption: {
          en: "DATE_RANGE",
          fr: "PLAGE_DE_DATES",
          pt: "INTERVALO_DE_DATAS",
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
    pt: "Cobertura (todos os tipos de estimativa)",
  },
  variantLabel: {
    en: "Admin area 3",
    fr: "Zone administrative 3",
    pt: "Área administrativa 3",
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
      pt: "Estimativas de cobertura ao nível distrital que combinam valores de inquéritos, projetados e derivados do SIS.",
    },
    methodology: {
      en: "AVG of survey, projected survey, and HMIS coverage at district level. Finest geographic resolution for coverage estimation.",
      fr: "Moyenne de la couverture d'enquête, d'enquête projetée et HMIS au niveau du district. Résolution géographique la plus fine.",
      pt: "Média da cobertura de inquéritos, de inquéritos projetados e do SIS ao nível distrital. Resolução geográfica mais fina para a estimativa de cobertura.",
    },
    interpretation: {
      en: "Enables district-level targeting and micro-planning. Interpret with caution due to smaller sample sizes and denominator uncertainty.",
      fr: "Permet le ciblage au niveau du district et la micro-planification. Interpréter avec prudence en raison de petites tailles d'échantillon.",
      pt: "Permite a focalização ao nível distrital e o micro-planeamento. Interpretar com prudência devido às menores dimensões das amostras e à incerteza dos denominadores.",
    },
    typicalRange: {
      en: "0-100%. Greater variation expected at district level; interpret with caution.",
      fr: "0-100%. Plus grande variation attendue au niveau du district; interpréter avec prudence.",
      pt: "0-100%. Espera-se maior variação ao nível distrital; interpretar com prudência.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id, admin_area_3, and year (all required). Compare with admin area 2 for context.",
      fr: "Toujours désagréger par indicator_common_id, admin_area_3 et year (tous requis). Comparer avec la zone administrative 2 pour le contexte.",
      pt: "Desagregar sempre por indicator_common_id, admin_area_3 e year (todos obrigatórios). Comparar com a área administrativa 2 para contextualização.",
    },
    caveats: null,
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
