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
      en: "Coverage timeseries by Admin Area 2",
      fr: "Séries temporelles de couverture par Zone administrative 2",
      pt: "Séries temporais de cobertura por Área Administrativa 2",
    },
    description: {
      en: "Coverage trends over time by admin area 2",
      fr: "Tendances de couverture dans le temps par zone administrative 2",
      pt: "Tendências de cobertura ao longo do tempo por área administrativa 2",
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
        content: "lines",
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Coverage estimates",
          fr: "Estimations de couverture",
          pt: "Estimativas de cobertura",
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
  {
    id: "coverage-bar",
    label: {
      en: "Coverage bar chart by Admin Area 2",
      fr: "Diagramme à barres de couverture par Zone administrative 2",
      pt: "Gráfico de barras de cobertura por Área Administrativa 2",
    },
    description: {
      en: "Bar chart comparing coverage across Admin Areas 2",
      fr: "Diagramme à barres comparant la couverture entre Zones administratives 2",
      pt: "Gráfico de barras comparando a cobertura entre Áreas Administrativas 2",
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
          pt: "Estimativas de cobertura ao nível subnacional, REPLICANT",
        },
        subCaption: {
          en: "DATE_RANGE",
          fr: "PLAGE_DE_DATES",
          pt: "INTERVALO_DE_DATAS",
        },
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
    pt: "Cobertura calculada a partir dos dados do HMIS",
  },
  variantLabel: {
    en: "Admin Area 2",
    fr: "Zone administrative 2",
    pt: "Área Administrativa 2",
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
      pt: "Cobertura dos serviços de saúde derivada do HMIS ao nível da área administrativa 2 (província/estado).",
    },
    methodology: {
      en: "AVG of coverage calculated as HMIS service volumes divided by subnational population denominators. Uses nationally-selected denominator with subnational fallback if national-only denominator chosen.",
      fr: "Moyenne de la couverture calculée comme volumes de services HMIS divisés par dénominateurs de population sous-nationale.",
      pt: "Média da cobertura calculada como volumes de serviços do HMIS divididos pelos denominadores da população subnacional. Utiliza o denominador selecionado a nível nacional, recorrendo ao subnacional caso seja escolhido um denominador apenas nacional.",
    },
    interpretation: {
      en: "Enables subnational coverage monitoring and equity analysis. Compare across regions to identify geographic disparities. Coverage >100% indicates denominator or data quality issues.",
      fr: "Permet la surveillance de la couverture sous-nationale et l'analyse de l'équité. Comparer entre régions pour identifier les disparités géographiques.",
      pt: "Permite a monitorização da cobertura subnacional e a análise da equidade. Comparar entre regiões para identificar disparidades geográficas. Uma cobertura >100% indica problemas no denominador ou na qualidade dos dados.",
    },
    typicalRange: {
      en: "0-100%. Regional variation expected; coverage gaps often larger in remote/underserved areas.",
      fr: "0-100%. Variation régionale attendue; écarts de couverture souvent plus grands dans les zones éloignées.",
      pt: "0-100%. Variação regional esperada; lacunas de cobertura frequentemente maiores em zonas remotas/desfavorecidas.",
    },
    caveats: {
      en: "Subnational denominators may be less reliable than national. Migration and population estimates affect accuracy. Some denominators only available at national level.",
      fr: "Les dénominateurs sous-nationaux peuvent être moins fiables que nationaux. Les estimations de migration et population affectent la précision.",
      pt: "Os denominadores subnacionais podem ser menos fiáveis do que os nacionais. As estimativas de migração e população afetam a precisão. Alguns denominadores apenas estão disponíveis a nível nacional.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id, admin_area_2, and year (all required). Map visualization effectively shows geographic coverage patterns. Time series reveals regional improvement or deterioration.",
      fr: "Toujours désagréger par indicator_common_id, admin_area_2 et year (tous requis). La visualisation cartographique montre efficacement les modèles de couverture géographique.",
      pt: "Desagregar sempre por indicator_common_id, admin_area_2 e year (todos obrigatórios). A visualização cartográfica mostra eficazmente os padrões geográficos de cobertura. As séries temporais revelam a melhoria ou deterioração regional.",
    },
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
