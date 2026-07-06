import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { CF_01_03 } from "../../.validation/cf_presets.ts";

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
    id: "outlier-table",
    label: {
      en: "Outlier proportion table",
      fr: "Tableau de proportion de valeurs aberrantes",
      pt: "Tabela de proporção de valores atípicos",
    },
    description: {
      en: "Table showing proportion of outliers by indicator and Admin Area 2",
      fr: "Tableau montrant la proportion de valeurs aberrantes par indicateur et Zone administrative 2",
      pt: "Tabela que mostra a proporção de valores atípicos por indicador e Área Administrativa 2",
    },
    createDefaultVisualizationOnInstall: "c3cb0cc9-4352-4b27-8532-f18e465faec8",
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
        ...CF_01_03,
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Outliers",
          fr: "Valeurs aberrantes",
          pt: "Valores atípicos",
        },
        subCaption: {
          en: "Percentage of facility-months that are outliers, DATE_RANGE",
          fr: "Pourcentage de mois-établissements qui sont des valeurs aberrantes, PLAGE_DE_DATES",
          pt: "Percentagem de meses-unidade sanitária que são valores atípicos, INTERVALO_DE_DATAS",
        },
        footnote: {
          en: "Outliers are reports which are suspiciously high compared to the usual volume reported by the facility in other months. Outliers are identified by assessing the within-facility variation in monthly reporting for each indicator. Outliers are defined observations which are greater than 10 times the median absolute deviation (MAD) from the monthly median value for the indicator in each time period, OR a value for which the proportional contribution in volume for a facility, indicator, and time period  is greater than 80%. Outliers are only identified for indicators where the volume is greater than or equal to the median, the volume is not missing, and the average volume is greater than 100.",
          fr: "Les valeurs aberrantes sont des rapports anormalement élevés par rapport au volume habituel déclaré par l'établissement au cours des autres mois. Elles sont identifiées en évaluant la variation intra-établissement des déclarations mensuelles pour chaque indicateur. Les valeurs aberrantes sont définies comme des observations supérieures à 10 fois l'écart absolu médian (MAD) par rapport à la valeur médiane mensuelle de l'indicateur pour chaque période, OU une valeur dont la contribution proportionnelle au volume pour un établissement, indicateur et période est supérieure à 80%. Les valeurs aberrantes ne sont identifiées que pour les indicateurs dont le volume est supérieur ou égal à la médiane, le volume n'est pas manquant, et le volume moyen est supérieur à 100.",
          pt: "Os valores atípicos são notificações suspeitosamente elevadas em comparação com o volume habitual notificado pela unidade sanitária nos outros meses. Os valores atípicos são identificados avaliando a variação intra-unidade na notificação mensal de cada indicador. Os valores atípicos são definidos como observações superiores a 10 vezes o desvio absoluto da mediana (MAD) relativamente ao valor mediano mensal do indicador em cada período de tempo, OU um valor cuja contribuição proporcional em volume para uma unidade sanitária, indicador e período de tempo é superior a 80%. Os valores atípicos só são identificados para indicadores cujo volume é maior ou igual à mediana, cujo volume não está em falta e cujo volume médio é superior a 100.",
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
    id: "outlier-map",
    label: {
      en: "Outlier proportion map",
      fr: "Carte de proportion de valeurs aberrantes",
      pt: "Mapa de proporção de valores atípicos",
    },
    description: {
      en: "Map showing proportion of outliers by indicator and Admin Area 2",
      fr: "Carte montrant la proportion de valeurs aberrantes par indicateur et Zone administrative 2",
      pt: "Mapa que mostra a proporção de valores atípicos por indicador e Área Administrativa 2",
    },
    createDefaultVisualizationOnInstall: "efc60afb-0b25-4e5e-95e4-e97754d3af83",
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
        ...CF_01_03,
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Outliers",
          fr: "Valeurs aberrantes",
          pt: "Valores atípicos",
        },
        subCaption: {
          en: "Percentage of facility-months that are outliers, DATE_RANGE",
          fr: "Pourcentage de mois-établissements qui sont des valeurs aberrantes, PLAGE_DE_DATES",
          pt: "Percentagem de meses-unidade sanitária que são valores atípicos, INTERVALO_DE_DATAS",
        },
        footnote: {
          en: "Outliers are reports which are suspiciously high compared to the usual volume reported by the facility in other months. Outliers are identified by assessing the within-facility variation in monthly reporting for each indicator. Outliers are defined observations which are greater than 10 times the median absolute deviation (MAD) from the monthly median value for the indicator in each time period, OR a value for which the proportional contribution in volume for a facility, indicator, and time period  is greater than 80%. Outliers are only identified for indicators where the volume is greater than or equal to the median, the volume is not missing, and the average volume is greater than 100.",
          fr: "Les valeurs aberrantes sont des rapports anormalement élevés par rapport au volume habituel déclaré par l'établissement au cours des autres mois. Elles sont identifiées en évaluant la variation intra-établissement des déclarations mensuelles pour chaque indicateur. Les valeurs aberrantes sont définies comme des observations supérieures à 10 fois l'écart absolu médian (MAD) par rapport à la valeur médiane mensuelle de l'indicateur pour chaque période, OU une valeur dont la contribution proportionnelle au volume pour un établissement, indicateur et période est supérieure à 80%. Les valeurs aberrantes ne sont identifiées que pour les indicateurs dont le volume est supérieur ou égal à la médiane, le volume n'est pas manquant, et le volume moyen est supérieur à 100.",
          pt: "Os valores atípicos são notificações suspeitosamente elevadas em comparação com o volume habitual notificado pela unidade sanitária nos outros meses. Os valores atípicos são identificados avaliando a variação intra-unidade na notificação mensal de cada indicador. Os valores atípicos são definidos como observações superiores a 10 vezes o desvio absoluto da mediana (MAD) relativamente ao valor mediano mensal do indicador em cada período de tempo, OU um valor cuja contribuição proporcional em volume para uma unidade sanitária, indicador e período de tempo é superior a 80%. Os valores atípicos só são identificados para indicadores cujo volume é maior ou igual à mediana, cujo volume não está em falta e cujo volume médio é superior a 100.",
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
  id: "m1-01-01",
  resultsObjectId: "M1_output_outliers.csv",
  valueProps: ["outlier_flag"],
  valueFunc: "AVG",
  valueLabelReplacements: {
    outlier_flag: "Binary variable indicating whether this an outlier",
  },
  label: {
    en: "Proportion of outliers",
    fr: "Proportion de valeurs aberrantes",
    pt: "Proporção de valores atípicos",
  },
  requiredDisaggregationOptions: [],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Proportion of data points flagged as statistical outliers in the dataset.",
      fr: "Proportion de points de données signalés comme valeurs aberrantes statistiques.",
      pt: "Proporção de pontos de dados assinalados como valores atípicos estatísticos no conjunto de dados.",
    },
    methodology: {
      en: "AVG of binary outlier_flag column. Outliers identified using Median Absolute Deviation (MAD) with configurable threshold (default: 10 MADs).",
      fr: "Moyenne de la colonne binaire outlier_flag. Valeurs aberrantes identifiées par l'écart absolu médian (MAD).",
      pt: "MÉDIA da coluna binária outlier_flag. Valores atípicos identificados através do desvio absoluto da mediana (MAD).",
    },
    interpretation: {
      en: "Higher values indicate more data quality issues. Values above 5% typically warrant investigation. Compare across indicators and regions to identify systematic problems.",
      fr: "Des valeurs plus élevées indiquent davantage de problèmes de qualité des données. Les valeurs supérieures à 5% nécessitent généralement une investigation.",
      pt: "Valores mais elevados indicam mais problemas de qualidade dos dados. Valores acima de 5% justificam normalmente uma investigação.",
    },
    typicalRange: {
      en: "0-5% for good quality data; 5-10% acceptable; >10% indicates significant issues.",
      fr: "0-5% pour des données de bonne qualité; 5-10% acceptable; >10% indique des problèmes significatifs.",
      pt: "0-5% para dados de boa qualidade; 5-10% aceitável; >10% indica problemas significativos.",
    },
    caveats: {
      en: "Threshold is configurable in module parameters. Compare results using consistent thresholds. Low reporting may mask outliers.",
      fr: "Le seuil est configurable. Comparez les résultats avec des seuils cohérents.",
      pt: "O limiar é configurável nos parâmetros do módulo. Compare os resultados utilizando limiares coerentes. Uma notificação baixa pode mascarar valores atípicos.",
    },
    disaggregationGuidance: {
      en: "Disaggregate by indicator_common_id to identify problem indicators. Use admin_area_2 for regional patterns. Combine with facility_type to see if certain facility types have more issues.",
      fr: "Désagréger par indicator_common_id pour identifier les indicateurs problématiques. Utiliser admin_area_2 pour les tendances régionales.",
      pt: "Desagregar por indicator_common_id para identificar os indicadores problemáticos. Utilizar admin_area_2 para os padrões regionais. Combinar com facility_type para verificar se certos tipos de unidade sanitária têm mais problemas.",
    },
  },
  variantLabel: null,
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
