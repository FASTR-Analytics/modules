import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { CF_DIVERGING_5_10_20 } from "../../.validation/cf_presets.ts";

export const vizPresets: VizPreset[] = [
  ///////////////////////////////////////////////
  //  ________         __        __            //
  // /        |       /  |      /  |           //
  // $$$$$$$$/______  $$ |____  $$ |  ______   //
  //    $$ | /      \ $$      \ $$ | /      \  //
  //    $$ | $$$$$$  |$$$$$$$  |$$ |/$$$$$$  | //
  //    $$ | /    $$ |$$ |  $$ |$$ |$$    $$ | //
  //    $$ |/$$$$$$$ |$$ |__$$ |$$ |$$$$$$$$/  //
  //    $$ |$$    $$ |$$    $$/ $$ |$$       | //
  //    $$/  $$$$$$$/ $$$$$$$/  $$/  $$$$$$$/  //
  //                                           //
  ///////////////////////////////////////////////
  {
    id: "disruption-differences-table-single-admin-area-3-multiple-admin-area-4",
    label: {
      en: "Service disruption by Admin Area 4",
      fr: "Perturbation des services par Zone administrative 4",
      pt: "Perturbação dos serviços por Área administrativa 4",
    },
    description: {
      en: "Table showing percentage difference between actual vs expected service volume for a single Admin Area 3, with multiple Admin Areas 4 and multiple indicators",
      fr: "Tableau montrant la différence en pourcentage entre le volume de services réel et attendu pour une unique Zone administrative 3, avec plusieurs Zones administratives 4 et plusieurs indicateurs",
      pt: "Tabela que mostra a diferença percentual entre o volume de serviços real e esperado para uma única Área administrativa 3, com várias Áreas administrativas 4 e vários indicadores",
    },
    allowedFilters: ["indicator_common_id", "admin_area_4"],
    config: {
      d: {
        type: "table",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "col",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "col",
          },
          {
            disOpt: "admin_area_3",
            disDisplayOpt: "replicant",
          },
          {
            disOpt: "admin_area_4",
            disDisplayOpt: "row",
          },
        ],
        selectedReplicantValue: "",
        filterBy: [],
        periodFilter: {
          filterType: "last_n_calendar_quarters",
          nQuarters: 1,
        },
      },
      s: {
        ...CF_DIVERGING_5_10_20,
      },
      t: {
        caption: null,
        captionRelFontSize: null,
        subCaption: {
          en: "DATE_RANGE",
          fr: "PLAGE_DE_DATES",
          pt: "INTERVALO_DE_DATAS",
        },
        subCaptionRelFontSize: null,
        footnote: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
    createDefaultVisualizationOnInstall: "45735c34-5dd4-43ee-baad-346c74751d4d",
  },
];

export const metric: MetricDefinitionGithub = {
  id: "m3-05-02",
  resultsObjectId: "M3_disruptions_analysis_admin_area_4.csv",
  label: {
    en: "Difference between actual and expected service volume (%)",
    fr: "Différence entre le volume de services réel et attendu (%)",
    pt: "Diferença entre o volume de serviços real e esperado (%)",
  },
  variantLabel: {
    en: "Admin area 4",
    fr: "Zone administrative 4",
    pt: "Área administrativa 4",
  },
  valueProps: ["pct_diff"],
  valueFunc: "identity",
  valueLabelReplacements: {
    pct_diff: "Percent difference",
  },
  postAggregationExpression: {
    ingredientValues: [
      {
        prop: "count_sum",
        func: "SUM",
      },
      {
        prop: "count_expect_sum",
        func: "SUM",
      },
    ],
    expression: "pct_diff = (count_sum - count_expect_sum)/count_expect_sum",
  },
  requiredDisaggregationOptions: ["indicator_common_id", "admin_area_4"],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Percentage difference between actual and expected service volumes at sub-district level.",
      fr: "Différence en pourcentage entre volumes réels et attendus au niveau sous-district.",
      pt: "Diferença percentual entre os volumes de serviços reais e esperados ao nível do subdistrito.",
    },
    methodology: {
      en: "(actual - expected) / expected at sub-district level. Finest geographic granularity for performance assessment.",
      fr: "(réel - attendu) / attendu au niveau sous-district. Granularité géographique la plus fine pour l'évaluation.",
      pt: "(real - esperado) / esperado ao nível do subdistrito. Granularidade geográfica mais fina para a avaliação do desempenho.",
    },
    interpretation: {
      en: "Enables targeted facility-level interventions. Small sample sizes at this level may increase volatility.",
      fr: "Permet des interventions ciblées au niveau de l'établissement. Les petites tailles d'échantillon peuvent augmenter la volatilité.",
      pt: "Permite intervenções direcionadas ao nível do estabelecimento. As pequenas dimensões da amostra a este nível podem aumentar a volatilidade.",
    },
    typicalRange: {
      en: "±20-50% variation common; high volatility expected at this granular level.",
      fr: "Variation de ±20-50% commune; forte volatilité attendue à ce niveau granulaire.",
      pt: "Variação de ±20-50% comum; alta volatilidade esperada a este nível granular.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id and admin_area_4 (both required). Interpret with caution due to small sample sizes.",
      fr: "Toujours désagréger par indicator_common_id et admin_area_4 (tous deux requis). Interpréter avec prudence en raison de petites tailles d'échantillon.",
      pt: "Desagregar sempre por indicator_common_id e admin_area_4 (ambos obrigatórios). Interpretar com cautela devido às pequenas dimensões da amostra.",
    },
    caveats: null,
  },
  importantNotes: null,
  hide: false,
  vizPresets,
};
