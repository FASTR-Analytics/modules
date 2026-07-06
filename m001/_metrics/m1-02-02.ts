import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { CF_90_80 } from "../../.validation/cf_presets.ts";

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
    id: "completeness-table",
    label: {
      en: "Completeness table by Admin Area 2",
      fr: "Tableau de complétude par Zone administrative 2",
      pt: "Tabela de completude por Área Administrativa 2",
    },
    description: {
      en: "Table showing completeness by indicator and region",
      fr: "Tableau montrant la complétude par indicateur et région",
      pt: "Tabela que mostra a completude por indicador e região",
    },
    createDefaultVisualizationOnInstall: "c20f1672-edfc-4140-ae2c-09a30b50443a",
    allowedFilters: ["indicator_common_id", "admin_area_2"],
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
            disOpt: "admin_area_2",
            disDisplayOpt: "row",
          },
        ],
        filterBy: [],
        periodFilter: {
          filterType: "last_n_months",
          nMonths: 12,
        },
      },
      s: {
        content: "lines",
        ...CF_90_80,
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Indicator Completeness",
          fr: "Complétude des indicateurs",
          pt: "Completude dos indicadores",
        },
        subCaption: {
          en: "Percentage of facility-months with complete data, DATE_RANGE",
          fr: "Pourcentage de mois-établissements avec des données complètes, PLAGE_DE_DATES",
          pt: "Percentagem de meses-unidade sanitária com dados completos, INTERVALO_DE_DATAS",
        },
        footnote: {
          en: "Higher completeness improves the reliability of the data, especially when completeness is stable over time. Completeness is defined as the percentage of reporting facilities each month out of the total number of facilities expected to report. A facility is expected to report if it has reported any volume for each indicator anytime within a year. A high completeness does not indicate that the HMIS is representative of all service delivery in the country, as some services may not be delivered in facilities, or some facilities may not report.",
          fr: "Une complétude élevée améliore la fiabilité des données, surtout lorsqu'elle est stable dans le temps. La complétude est définie comme le pourcentage d'établissements déclarants chaque mois par rapport au nombre total d'établissements censés déclarer. Un établissement est censé déclarer s'il a déclaré un volume pour chaque indicateur à tout moment au cours de l'année. Une complétude élevée n'indique pas que le HMIS est représentatif de toute la prestation de services dans le pays, car certains services peuvent ne pas être fournis dans les établissements, ou certains établissements peuvent ne pas déclarer.",
          pt: "Uma completude mais elevada melhora a fiabilidade dos dados, sobretudo quando a completude é estável ao longo do tempo. A completude é definida como a percentagem de unidades sanitárias notificadoras em cada mês face ao número total de unidades sanitárias que se espera que notifiquem. Espera-se que uma unidade sanitária notifique se tiver notificado algum volume para cada indicador em qualquer momento ao longo de um ano. Uma completude elevada não indica que o SIS seja representativo de toda a prestação de serviços no país, pois alguns serviços podem não ser prestados nas unidades sanitárias, ou algumas unidades sanitárias podem não notificar.",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
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
    id: "completeness-timeseries",
    label: {
      en: "Completeness over time",
      fr: "Complétude dans le temps",
      pt: "Completude ao longo do tempo",
    },
    description: {
      en: "Area chart showing completeness trends over time by indicator",
      fr: "Graphique en aires montrant les tendances de complétude dans le temps par indicateur",
      pt: "Gráfico de áreas que mostra as tendências de completude ao longo do tempo por indicador",
    },
    createDefaultVisualizationOnInstall: "26dedd7c-4577-4022-928c-69e0ee790a71",
    allowedFilters: ["indicator_common_id"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "row",
          },
        ],
        filterBy: [],
      },
      s: {
        content: "lines-area",
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Indicator completeness over time",
          fr: "Complétude des indicateurs dans le temps",
          pt: "Completude dos indicadores ao longo do tempo",
        },
        subCaption: {
          en: "Percentage of facility-months with complete data DATE_RANGE",
          fr: "Pourcentage de mois-établissements avec des données complètes PLAGE_DE_DATES",
          pt: "Percentagem de meses-unidade sanitária com dados completos INTERVALO_DE_DATAS",
        },
        footnote: {
          en: "Higher completeness improves the reliability of the data, especially when completeness is stable over time. Completeness is defined as the percentage of reporting facilities each month out of the total number of facilities expected to report. A facility is expected to report if it has reported any volume for each indicator anytime within a year. A high completeness does not indicate that the HMIS is representative of all service delivery in the country, as some services may not be delivered in facilities, or some facilities may not report.",
          fr: "Une complétude élevée améliore la fiabilité des données, surtout lorsqu'elle est stable dans le temps. La complétude est définie comme le pourcentage d'établissements déclarants chaque mois par rapport au nombre total d'établissements censés déclarer. Un établissement est censé déclarer s'il a déclaré un volume pour chaque indicateur à tout moment au cours de l'année. Une complétude élevée n'indique pas que le HMIS est représentatif de toute la prestation de services dans le pays, car certains services peuvent ne pas être fournis dans les établissements, ou certains établissements peuvent ne pas déclarer.",
          pt: "Uma completude mais elevada melhora a fiabilidade dos dados, sobretudo quando a completude é estável ao longo do tempo. A completude é definida como a percentagem de unidades sanitárias notificadoras em cada mês face ao número total de unidades sanitárias que se espera que notifiquem. Espera-se que uma unidade sanitária notifique se tiver notificado algum volume para cada indicador em qualquer momento ao longo de um ano. Uma completude elevada não indica que o SIS seja representativo de toda a prestação de serviços no país, pois alguns serviços podem não ser prestados nas unidades sanitárias, ou algumas unidades sanitárias podem não notificar.",
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
    id: "completeness-map",
    label: {
      en: "Completeness map by Admin Area 2",
      fr: "Carte de complétude par Zone administrative 2",
      pt: "Mapa de completude por Área Administrativa 2",
    },
    description: {
      en: "Map showing completeness by indicator and region",
      fr: "Carte montrant la complétude par indicateur et région",
      pt: "Mapa que mostra a completude por indicador e região",
    },
    createDefaultVisualizationOnInstall: "f9d926c9-e9d9-41b8-bca9-c3ae252ab259",
    allowedFilters: [],
    config: {
      d: {
        type: "map",
        valuesDisDisplayOpt: "mapArea",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "replicant",
          },
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "mapArea",
          },
        ],
        selectedReplicantValue: "anc1",
        filterBy: [],
        periodFilter: {
          filterType: "last_n_months",
          nMonths: 12,
        },
      },
      s: {
        ...CF_90_80,
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Indicator Completeness",
          fr: "Complétude des indicateurs",
          pt: "Completude dos indicadores",
        },
        subCaption: {
          en: "Percentage of facility-months with complete data, DATE_RANGE",
          fr: "Pourcentage de mois-établissements avec des données complètes, PLAGE_DE_DATES",
          pt: "Percentagem de meses-unidade sanitária com dados completos, INTERVALO_DE_DATAS",
        },
        footnote: {
          en: "Higher completeness improves the reliability of the data, especially when completeness is stable over time. Completeness is defined as the percentage of reporting facilities each month out of the total number of facilities expected to report. A facility is expected to report if it has reported any volume for each indicator anytime within a year. A high completeness does not indicate that the HMIS is representative of all service delivery in the country, as some services may not be delivered in facilities, or some facilities may not report.",
          fr: "Une complétude élevée améliore la fiabilité des données, surtout lorsqu'elle est stable dans le temps. La complétude est définie comme le pourcentage d'établissements déclarants chaque mois par rapport au nombre total d'établissements censés déclarer. Un établissement est censé déclarer s'il a déclaré un volume pour chaque indicateur à tout moment au cours de l'année. Une complétude élevée n'indique pas que le HMIS est représentatif de toute la prestation de services dans le pays, car certains services peuvent ne pas être fournis dans les établissements, ou certains établissements peuvent ne pas déclarer.",
          pt: "Uma completude mais elevada melhora a fiabilidade dos dados, sobretudo quando a completude é estável ao longo do tempo. A completude é definida como a percentagem de unidades sanitárias notificadoras em cada mês face ao número total de unidades sanitárias que se espera que notifiquem. Espera-se que uma unidade sanitária notifique se tiver notificado algum volume para cada indicador em qualquer momento ao longo de um ano. Uma completude elevada não indica que o SIS seja representativo de toda a prestação de serviços no país, pois alguns serviços podem não ser prestados nas unidades sanitárias, ou algumas unidades sanitárias podem não notificar.",
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
  id: "m1-02-02",
  resultsObjectId: "M1_output_completeness.csv",
  valueProps: ["completeness_flag"],
  valueFunc: "AVG",
  valueLabelReplacements: {
    completeness_flag:
      "Binary variable indicating whether the facility meets criteria",
  },
  label: {
    en: "Proportion of completed records",
    fr: "Proportion d'enregistrements complets",
    pt: "Proporção de registos completos",
  },
  requiredDisaggregationOptions: [],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Proportion of facility-indicator-period combinations meeting completeness criteria.",
      fr: "Proportion de combinaisons établissement-indicateur-période répondant aux critères de complétude.",
      pt: "Proporção de combinações unidade sanitária-indicador-período que cumprem os critérios de completude.",
    },
    methodology: {
      en: "AVG of binary completeness_flag. Facilities must report consistently across expected periods to be flagged as complete.",
      fr: "Moyenne du drapeau binaire de complétude. Les établissements doivent déclarer régulièrement pour être considérés comme complets.",
      pt: "MÉDIA do indicador binário completeness_flag. As unidades sanitárias têm de notificar de forma coerente ao longo dos períodos esperados para serem assinaladas como completas.",
    },
    interpretation: {
      en: "Higher values indicate better reporting consistency. Values below 80% suggest significant reporting gaps that may bias analysis.",
      fr: "Des valeurs plus élevées indiquent une meilleure cohérence des déclarations. Les valeurs inférieures à 80% suggèrent des lacunes significatives.",
      pt: "Valores mais elevados indicam uma melhor coerência da notificação. Valores abaixo de 80% sugerem lacunas significativas na notificação que podem enviesar a análise.",
    },
    typicalRange: {
      en: "80-100% is good; 60-80% moderate; <60% indicates major gaps.",
      fr: "80-100% est bon; 60-80% modéré; <60% indique des lacunes majeures.",
      pt: "80-100% é bom; 60-80% moderado; <60% indica lacunas importantes.",
    },
    caveats: {
      en: "Definition of completeness can vary. Check module parameters for specific criteria used.",
      fr: "La définition de complétude peut varier. Vérifiez les paramètres du module.",
      pt: "A definição de completude pode variar. Verifique os parâmetros do módulo para conhecer os critérios específicos utilizados.",
    },
    disaggregationGuidance: {
      en: "Disaggregate by admin_area to identify regions with reporting challenges. Use indicator_common_id to see if specific services have lower compliance.",
      fr: "Désagréger par zone administrative pour identifier les régions avec des défis de déclaration.",
      pt: "Desagregar por admin_area para identificar regiões com dificuldades de notificação. Utilizar indicator_common_id para verificar se serviços específicos têm menor conformidade.",
    },
  },
  variantLabel: null,
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
