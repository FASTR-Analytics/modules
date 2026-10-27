import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "adjustment-comparison",
    label: {
      en: "Data quality adjustment comparison",
      fr: "Comparaison des ajustements de qualité des données",
      pt: "Comparação dos ajustes de qualidade dos dados",
    },
    description: {
      en:
        "Line chart comparing annual service volume under each adjustment, one panel per indicator",
      fr:
        "Graphique linéaire comparant le volume annuel de services selon chaque ajustement, un panneau par indicateur",
      pt:
        "Gráfico de linhas que compara o volume anual de serviços segundo cada ajuste, um painel por indicador",
    },
    allowedFilters: ["indicator_common_id", "admin_area_2"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [{
          disOpt: "indicator_common_id",
          disDisplayOpt: "col",
        }],
        filterBy: [],
      },
      s: {
        colorScale: "custom",
        customSeriesStyles: [
          { color: "#00897b", lineStyle: "solid", strokeWidth: 5 },
          { color: "#757575", lineStyle: "solid", strokeWidth: 5 },
          { color: "#8e24aa", lineStyle: "solid", strokeWidth: 5 },
          { color: "#7cb342", lineStyle: "solid", strokeWidth: 5 },
        ],
      },
      t: {
        caption: {
          en: "Change in volume due to data quality adjustments",
          fr: "Variation du volume due aux ajustements de qualité des données",
          pt: "Variação do volume devida aos ajustes de qualidade dos dados",
        },
        subCaption: {
          en: "DATE_RANGE",
          fr: "PLAGE_DE_DATES",
          pt: "INTERVALO_DE_DATAS",
        },
        footnote: null,
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
];

// The four adjustment bases exist side by side only here, at facility grain:
// every downstream module runs on the one basis its package chose.
export const metric: MetricDefinitionGithub = {
  id: "m2-02-01",
  resultsObjectId: "M2_adjusted_data.csv",
  valueProps: [
    "count_final_none",
    "count_final_outliers",
    "count_final_completeness",
    "count_final_both",
  ],
  valueFunc: "SUM",
  formatAs: "number",
  label: {
    en: "Service volume by adjustment type",
    fr: "Volume de services par type d'ajustement",
    pt: "Volume de serviços por tipo de ajuste",
  },
  variantLabel: null,
  requiredDisaggregationOptions: ["indicator_common_id"],
  valueLabelReplacements: {
    count_final_none: "Reported volume",
    count_final_outliers: "Adjusted for outliers",
    count_final_completeness: "Adjusted for completeness",
    count_final_both: "Adjusted for outliers and completeness",
  },
  postAggregationExpression: null,
  catalogExpressionEvaluation: null,
  aiDescription: {
    summary: {
      en:
        "Service volume as reported and after each data quality adjustment, side by side.",
      fr:
        "Volume de services tel que déclaré et après chaque ajustement de qualité des données, côte à côte.",
      pt:
        "Volume de serviços tal como notificado e após cada ajuste de qualidade dos dados, lado a lado.",
    },
    methodology: {
      en:
        "Sums the four count variants of the adjusted facility data: as reported, with outliers replaced, with missing reports imputed, and with both. The adjustments are the module's own outlier and completeness corrections.",
      fr:
        "Somme les quatre variantes de comptage des données ajustées par établissement : telles que déclarées, avec les valeurs aberrantes remplacées, avec les rapports manquants imputés, et avec les deux. Les ajustements sont les corrections des valeurs aberrantes et de la complétude du module.",
      pt:
        "Soma as quatro variantes de contagem dos dados ajustados por unidade sanitária: tal como notificados, com valores atípicos substituídos, com relatórios em falta imputados, e com ambos. Os ajustes são as correções de valores atípicos e de completude do módulo.",
    },
    interpretation: {
      en:
        "The gap between the reported line and an adjusted line is the volume that adjustment adds or removes. Lines that sit close together mean the adjustments change little.",
      fr:
        "L'écart entre la ligne déclarée et une ligne ajustée est le volume que cet ajustement ajoute ou retire. Des lignes proches signifient que les ajustements changent peu.",
      pt:
        "A diferença entre a linha notificada e uma linha ajustada é o volume que esse ajuste acrescenta ou retira. Linhas próximas significam que os ajustes alteram pouco.",
    },
    typicalRange: {
      en:
        "Adjusted volumes usually lie within a few percent of the reported volume; completeness adjustment raises the volume, outlier adjustment lowers it.",
      fr:
        "Les volumes ajustés se situent généralement à quelques pour cent du volume déclaré ; l'ajustement de la complétude augmente le volume, celui des valeurs aberrantes le diminue.",
      pt:
        "Os volumes ajustados situam-se normalmente a poucos pontos percentuais do volume notificado; o ajuste de completude aumenta o volume, o de valores atípicos diminui-o.",
    },
    caveats: {
      en:
        "Facility-level rows are large; filter by indicator before disaggregating by facility.",
      fr:
        "Les lignes au niveau des établissements sont nombreuses ; filtrer par indicateur avant de désagréger par établissement.",
      pt:
        "As linhas ao nível das unidades sanitárias são numerosas; filtrar por indicador antes de desagregar por unidade sanitária.",
    },
    disaggregationGuidance: {
      en:
        "Always by indicator_common_id (required). By year for the comparison chart; by admin area or facility to find where adjustments matter most.",
      fr:
        "Toujours par indicator_common_id (requis). Par année pour le graphique de comparaison ; par zone administrative ou établissement pour voir où les ajustements comptent le plus.",
      pt:
        "Sempre por indicator_common_id (obrigatório). Por ano para o gráfico de comparação; por área administrativa ou unidade sanitária para ver onde os ajustes mais importam.",
    },
  },
  importantNotes: null,
  hide: false,
  vizPresets,
};
