import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "scorecard-table",
    label: {
      en: "Scorecard table",
      fr: "Tableau de bord",
      pt: "Tabela de painel de indicadores",
    },
    description: {
      en: "Table showing calculated indicators with threshold-based coloring",
      fr: "Tableau montrant les indicateurs calculés avec coloration basée sur les seuils",
      pt: "Tabela que mostra os indicadores calculados com coloração baseada em limiares",
    },
    createDefaultVisualizationOnInstall: "b8c128c6-31e9-4481-8253-8cda7bc2bd74",
    allowedFilters: ["indicator_common_id", "admin_area_2", "admin_area_3"],
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
        specialScorecardTable: true,
      },
      t: {
        caption: {
          en: "Health Sector Scorecard",
          fr: "Tableau de bord du secteur santé",
          pt: "Painel de Indicadores do Sector da Saúde",
        },
        subCaption: {
          en: "Calculated indicators by area, DATE_RANGE",
          fr: "Indicateurs calculés par zone, PLAGE_DE_DATES",
          pt: "Indicadores calculados por área, INTERVALO_DE_DATAS",
        },
        footnote: {
          en: "",
          fr: "",
          pt: "",
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
  id: "m8-01-01",
  hide: false,
  resultsObjectId: "M8_output_scorecard.csv",
  valueProps: ["value"],
  valueFunc: "identity",
  valueLabelReplacements: {},
  label: {
    en: "Scorecard",
    fr: "Scorecard",
    pt: "Scorecard",
  },
  requiredDisaggregationOptions: ["indicator_common_id"],
  formatAs: "percent",
  postAggregationExpression: {
    ingredientValues: [
      { prop: "numerator", func: "SUM" },
      { prop: "denominator", func: "SUM" },
    ],
    expression: "value = numerator / denominator",
  },
  aiDescription: {
    summary: {
      en: "Health sector scorecard indicators computed from the calculated indicators catalog.",
      fr: "Indicateurs de scorecard du secteur santé calculés à partir du catalogue des indicateurs calculés.",
      pt: "Indicadores do painel do sector da saúde calculados a partir do catálogo de indicadores calculados.",
    },
    methodology: {
      en: "Each indicator is computed as numerator/denominator. Denominators are either other HMIS indicators or population-based (annual population fraction scaled to monthly).",
      fr: "Chaque indicateur est calculé comme numérateur/dénominateur. Les dénominateurs sont soit d'autres indicateurs HMIS, soit basés sur la population (fraction annuelle de population mise à l'échelle mensuelle).",
      pt: "Cada indicador é calculado como numerador/denominador. Os denominadores são outros indicadores do SIS ou baseados na população (fracção anual da população ajustada à escala mensal).",
    },
    interpretation: {
      en: "Values represent coverage rates or ratios. Thresholds for traffic-light coloring come from the calculated indicators catalog.",
      fr: "Les valeurs représentent des taux de couverture ou des ratios. Les seuils pour la coloration feu tricolore proviennent du catalogue des indicateurs calculés.",
      pt: "Os valores representam taxas de cobertura ou rácios. Os limiares para a coloração tipo semáforo provêm do catálogo de indicadores calculados.",
    },
    typicalRange: {
      en: "Varies by indicator - see catalog for thresholds.",
      fr: "Varie selon l'indicateur - voir le catalogue pour les seuils.",
      pt: "Varia consoante o indicador - ver o catálogo para os limiares.",
    },
    disaggregationGuidance: {
      en: "Can be disaggregated by admin area (2, 3, or 4) and time period (month, quarter, year). Aggregation uses SUM(numerator)/SUM(denominator).",
      fr: "Peut être désagrégé par zone administrative (2, 3 ou 4) et période de temps (mois, trimestre, année). L'agrégation utilise SUM(numérateur)/SUM(dénominateur).",
      pt: "Pode ser desagregado por área administrativa (2, 3 ou 4) e período de tempo (mês, trimestre, ano). A agregação utiliza SUM(numerador)/SUM(denominador).",
    },
    caveats: null,
  },
  variantLabel: null,
  importantNotes: null,
  vizPresets,
};
