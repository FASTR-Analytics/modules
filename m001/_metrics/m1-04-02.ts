import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { CF_80_70 } from "../../.validation/cf_presets.ts";

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
    id: "mean-dqa-table",
    label: {
      en: "Mean DQA score table",
      fr: "Tableau du score EQD moyen",
      pt: "Tabela da pontuação média de AQD",
    },
    description: {
      en: "Table showing mean DQA scores by Admin Area 2 and year",
      fr: "Tableau montrant les scores EQD moyens par Zone administrative 2 et année",
      pt: "Tabela que mostra as pontuações médias de AQD por Área Administrativa 2 e ano",
    },
    createDefaultVisualizationOnInstall: "4dc02c21-29da-4a01-9812-469deedaaac8",
    allowedFilters: ["admin_area_2"],
    config: {
      d: {
        type: "table",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "col",
        disaggregateBy: [
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "row",
          },
          {
            disOpt: "year",
            disDisplayOpt: "col",
          },
        ],
        filterBy: [],
      },
      s: {
        content: "lines",
        ...CF_80_70,
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Mean DQA score",
          fr: "Score EQD moyen",
          pt: "Pontuação média de AQD",
        },
        subCaption: {
          en: "Average data quality score across facility-months, DATE_RANGE",
          fr: "Score moyen de qualité des données à travers les mois-établissements, PLAGE_DE_DATES",
          pt: "Pontuação média da qualidade dos dados ao longo dos meses-unidade sanitária, INTERVALO_DE_DATAS",
        },
        footnote: {
          en: "Items included in the DQA score include: No missing data for 1) OPD, 2) Penta1, and 3) ANC1, where available; No outliers for 4) OPD, 5) Penta1, and 6) ANC1, where available; Consistent reporting between 7) Penta1/Penta3, 8) ANC1/ANC4, 9)BCG/Delivery, where available.",
          fr: "Les éléments inclus dans le score EQD comprennent : Pas de données manquantes pour 1) OPD, 2) Penta1 et 3) ANC1, lorsque disponibles ; Pas de valeurs aberrantes pour 4) OPD, 5) Penta1 et 6) ANC1, lorsque disponibles ; Déclaration cohérente entre 7) Penta1/Penta3, 8) ANC1/ANC4, 9) BCG/Accouchement, lorsque disponibles.",
          pt: "Os elementos incluídos na pontuação de AQD são: Ausência de dados em falta para 1) OPD, 2) Penta1 e 3) ANC1, quando disponíveis; Ausência de valores atípicos para 4) OPD, 5) Penta1 e 6) ANC1, quando disponíveis; Notificação coerente entre 7) Penta1/Penta3, 8) ANC1/ANC4, 9) BCG/Parto, quando disponíveis.",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
  ///////////////////////////////////////
  //  __       __                      //
  // /  \     /  |                     //
  // $$  \   /$$ |  ______    ______   //
  // $$$  \ /$$$ | /      \  /      \  //
  // $$$$  /$$$$ | $$$$$$  |/$$$$$$  | //
  // $$ $$ $$/$$ | /    $$ |$$ |  $$ | //
  // $$ |$$$/ $$ |/$$$$$$$ |$$ |__$$ | //
  // $$ | $/  $$ |$$    $$ |$$    $$/  //
  // $$/      $$/  $$$$$$$/ $$$$$$$/   //
  //                        $$ |       //
  //                        $$ |       //
  //                        $$/        //
  //                                   //
  ///////////////////////////////////////
  {
    id: "mean-dqa-map",
    label: {
      en: "Mean DQA score map",
      fr: "Carte du score EQD moyen",
      pt: "Mapa da pontuação média de AQD",
    },
    description: {
      en: "Map showing mean DQA scores by Admin Area 2",
      fr: "Carte montrant les scores EQD moyens par Zone administrative 2",
      pt: "Mapa que mostra as pontuações médias de AQD por Área Administrativa 2",
    },
    createDefaultVisualizationOnInstall: "0730bd53-f951-4286-83ef-5b299395912c",
    allowedFilters: [],
    config: {
      d: {
        type: "map",
        valuesDisDisplayOpt: "mapArea",
        disaggregateBy: [
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "mapArea",
          },
        ],
        filterBy: [],
        periodFilter: {
          filterType: "last_n_months",
          nMonths: 12,
        },
      },
      s: {
        ...CF_80_70,
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Mean DQA score",
          fr: "Score EQD moyen",
          pt: "Pontuação média de AQD",
        },
        subCaption: {
          en: "Average data quality score across facility-months, DATE_RANGE",
          fr: "Score moyen de qualité des données à travers les mois-établissements, PLAGE_DE_DATES",
          pt: "Pontuação média da qualidade dos dados ao longo dos meses-unidade sanitária, INTERVALO_DE_DATAS",
        },
        footnote: {
          en: "Items included in the DQA score include: No missing data for 1) OPD, 2) Penta1, and 3) ANC1, where available; No outliers for 4) OPD, 5) Penta1, and 6) ANC1, where available; Consistent reporting between 7) Penta1/Penta3, 8) ANC1/ANC4, 9)BCG/Delivery, where available.",
          fr: "Les éléments inclus dans le score EQD comprennent : Pas de données manquantes pour 1) OPD, 2) Penta1 et 3) ANC1, lorsque disponibles ; Pas de valeurs aberrantes pour 4) OPD, 5) Penta1 et 6) ANC1, lorsque disponibles ; Déclaration cohérente entre 7) Penta1/Penta3, 8) ANC1/ANC4, 9) BCG/Accouchement, lorsque disponibles.",
          pt: "Os elementos incluídos na pontuação de AQD são: Ausência de dados em falta para 1) OPD, 2) Penta1 e 3) ANC1, quando disponíveis; Ausência de valores atípicos para 4) OPD, 5) Penta1 e 6) ANC1, quando disponíveis; Notificação coerente entre 7) Penta1/Penta3, 8) ANC1/ANC4, 9) BCG/Parto, quando disponíveis.",
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
  id: "m1-04-02",
  resultsObjectId: "M1_output_dqa.csv",
  valueProps: ["dqa_mean"],
  valueFunc: "AVG",
  valueLabelReplacements: {
    dqa_mean: "Data quality score across facilities",
  },
  label: {
    en: "Average data quality score across facilities",
    fr: "Score moyen de qualité des données à travers les établissements",
    pt: "Pontuação média da qualidade dos dados nas unidades sanitárias",
  },
  requiredDisaggregationOptions: [],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Average composite data quality score across all facilities.",
      fr: "Score moyen composite de qualité des données à travers tous les établissements.",
      pt: "Pontuação média composta da qualidade dos dados em todas as unidades sanitárias.",
    },
    methodology: {
      en: "AVG of continuous dqa_mean score. Combines completeness, outlier, and consistency assessments.",
      fr: "Moyenne du score dqa_mean continu. Combine les évaluations de complétude, valeurs aberrantes et cohérence.",
      pt: "MÉDIA da pontuação contínua dqa_mean. Combina as avaliações de completude, de valores atípicos e de coerência.",
    },
    interpretation: {
      en: "Higher values indicate better overall data quality. Use alongside m1-04-01 to understand both average performance and threshold compliance.",
      fr: "Des valeurs plus élevées indiquent une meilleure qualité globale des données.",
      pt: "Valores mais elevados indicam uma melhor qualidade global dos dados. Utilizar em conjunto com m1-04-01 para compreender tanto o desempenho médio como a conformidade com o limiar.",
    },
    typicalRange: {
      en: "0.7-1.0 is good; 0.5-0.7 moderate; <0.5 indicates significant issues.",
      fr: "0.7-1.0 est bon; 0.5-0.7 modéré; <0.5 indique des problèmes significatifs.",
      pt: "0,7-1,0 é bom; 0,5-0,7 moderado; <0,5 indica problemas significativos.",
    },
    disaggregationGuidance: {
      en: "Disaggregate by admin_area for regional comparison. Use time series to track improvement over time.",
      fr: "Désagréger par zone administrative pour comparaison régionale.",
      pt: "Desagregar por admin_area para comparação regional. Utilizar séries temporais para acompanhar a melhoria ao longo do tempo.",
    },
    caveats: null,
  },
  variantLabel: null,
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
