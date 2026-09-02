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
      en: "Table of indicators by area with threshold-based colouring",
      fr: "Tableau des indicateurs par zone avec coloration basée sur les seuils",
      pt: "Tabela de indicadores por área com coloração baseada em limiares",
    },
    createDefaultVisualizationOnInstall: "4d1f6a02-9c31-4c8a-9f2c-1b7a5e0d3c44",
    allowedFilters: ["indicator_common_id", "admin_area_2", "admin_area_3"],
    config: {
      d: {
        type: "table",
        valuesDisDisplayOpt: "col",
        disaggregateBy: [
          { disOpt: "indicator_common_id", disDisplayOpt: "col" },
          { disOpt: "admin_area_2", disDisplayOpt: "row" },
        ],
        filterBy: [],
        periodFilter: { filterType: "last_n_months", nMonths: 12 },
      },
      s: { specialScorecardTable: true },
      t: {
        caption: {
          en: "Health Sector Scorecard",
          fr: "Tableau de bord du secteur santé",
          pt: "Painel de Indicadores do Sector da Saúde",
        },
        subCaption: {
          en: "Indicators by area, DATE_RANGE",
          fr: "Indicateurs par zone, PLAGE_DE_DATES",
          pt: "Indicadores por área, INTERVALO_DE_DATAS",
        },
        footnote: { en: "", fr: "", pt: "" },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
];

export const metric: MetricDefinitionGithub = {
  id: "m12-01-01",
  hide: false,
  resultsObjectId: "M12_indicator_values.csv",
  // One displayed value. The ingredient columns are wire-only: the server sums
  // them and applies the indicator's catalog expression, so they are never a
  // user-facing prop and no value filter applies to them.
  valueProps: ["value"],
  valueFunc: "identity",
  valueLabelReplacements: {},
  label: {
    en: "Indicator value",
    fr: "Valeur de l'indicateur",
    pt: "Valor do indicador",
  },
  // REQUIRED, and it is what makes cross-indicator pooling impossible: the
  // server enforces required options as GROUP BYs, so every aggregated row the
  // evaluator sees is keyed by exactly one indicator. AUTHORING INVARIANT: the
  // guard is the INTERSECTION across all metrics sharing a results object, so
  // every metric ever added to M12_indicator_values.csv must require this too,
  // or it dissolves for all of them.
  requiredDisaggregationOptions: ["indicator_common_id"],
  formatAs: "indicator",
  postAggregationExpression: null,
  // The declared evaluation: these props travel on the wire, SUM-aggregated,
  // and the server replaces them with a single `value` computed from each
  // indicator's own catalog expression.
  catalogExpressionEvaluation: {
    ingredientProps: [
      "ing1",
      "ing2",
      "ing3",
      "ing4",
      "ing5",
      "ing6",
      "ing7",
      "ing8",
    ],
  },
  aiDescription: {
    summary: {
      en: "Every common indicator's value, computed from the instance's indicator dictionary.",
      fr: "La valeur de chaque indicateur commun, calculée à partir du dictionnaire d'indicateurs de l'instance.",
      pt: "O valor de cada indicador comum, calculado a partir do dicionário de indicadores da instância.",
    },
    methodology: {
      en: "The module materialises each indicator's additive ingredients at admin area x month grain. A chart re-sums those ingredients at whatever grouping it shows and then applies the indicator's formula once, so a regional or annual figure is the formula over the summed parts, never an average of ratios.",
      fr: "Le module matérialise les composantes additives de chaque indicateur au niveau zone administrative x mois. Un graphique resomme ces composantes au niveau qu'il affiche puis applique la formule une seule fois : un chiffre régional ou annuel est donc la formule appliquée aux sommes, jamais une moyenne de ratios.",
      pt: "O módulo materializa as componentes aditivas de cada indicador ao nível área administrativa x mês. Um gráfico volta a somar essas componentes no agrupamento que mostra e aplica a fórmula uma só vez, pelo que um valor regional ou anual é a fórmula sobre as somas e nunca uma média de rácios.",
    },
    interpretation: {
      en: "Values carry each indicator's own format (count, percent, or rate). Thresholds for traffic-light colouring come from the indicator dictionary. Population rates are ANNUALISED: the numerator is divided by person-years of the population (annual population × months / 12) and scaled by the indicator's multiplier, so a monthly or quarterly value reads as a rate per year and is directly comparable with an annual one.",
      fr: "Les valeurs portent le format propre à chaque indicateur (effectif, pourcentage ou taux). Les seuils de coloration feu tricolore proviennent du dictionnaire d'indicateurs. Les taux de population sont ANNUALISÉS : le numérateur est divisé par les personnes-années de la population (population annuelle × mois / 12) puis multiplié par le multiplicateur de l'indicateur, de sorte qu'une valeur mensuelle ou trimestrielle se lit comme un taux annuel, directement comparable à une valeur annuelle.",
      pt: "Os valores têm o formato próprio de cada indicador (contagem, percentagem ou taxa). Os limiares para a coloração tipo semáforo provêm do dicionário de indicadores. As taxas populacionais são ANUALIZADAS: o numerador é dividido pelas pessoas-ano da população (população anual × meses / 12) e multiplicado pelo multiplicador do indicador, pelo que um valor mensal ou trimestral se lê como uma taxa anual, diretamente comparável com um valor anual.",
    },
    typicalRange: {
      en: "Varies by indicator — see the indicator dictionary.",
      fr: "Varie selon l'indicateur — voir le dictionnaire d'indicateurs.",
      pt: "Varia consoante o indicador — ver o dicionário de indicadores.",
    },
    disaggregationGuidance: {
      en: "Always grouped by indicator. Can also be disaggregated by admin area (2, 3 or 4) and time period (month, quarter, year). Facility-level analysis is not available here — use m3-01-01.",
      fr: "Toujours groupé par indicateur. Peut aussi être désagrégé par zone administrative (2, 3 ou 4) et par période (mois, trimestre, année). L'analyse au niveau des établissements n'est pas disponible ici — utilisez m3-01-01.",
      pt: "Sempre agrupado por indicador. Pode também ser desagregado por área administrativa (2, 3 ou 4) e por período (mês, trimestre, ano). A análise ao nível dos estabelecimentos não está disponível aqui — utilize m3-01-01.",
    },
    caveats: {
      en: "Population denominators are interpolated between the instance's annual figures (mid-year anchors) and extrapolated at most one year beyond them, so a population rate's precision is that of the population estimates, not of the facility counts.",
      fr: "Les dénominateurs de population sont interpolés entre les chiffres annuels de l'instance (ancrés en milieu d'année) et extrapolés d'au plus un an au-delà ; la précision d'un taux de population est donc celle des estimations de population, pas celle des effectifs des établissements.",
      pt: "Os denominadores de população são interpolados entre os valores anuais da instância (âncoras a meio do ano) e extrapolados no máximo um ano para além deles, pelo que a precisão de uma taxa populacional é a das estimativas de população, não a das contagens dos estabelecimentos.",
    },
  },
  variantLabel: null,
  importantNotes: null,
  vizPresets,
};
