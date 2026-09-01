import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { CF_80_70 } from "../../.validation/cf_presets.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "hfa-variant-indicator-item-table-carried",
    label: {
      en: "Indicators by variant item, with carry-forward (per category)",
      fr: "Indicateurs par élément de variante, avec valeurs reportées (par catégorie)",
      pt: "Indicadores por item de variante, com valores transportados (por categoria)",
    },
    description: {
      en: "Table of per-item variant values for a selected category, with indicators as rows and variant items as columns, one column group per time point. Switch the category using the replicant selector. Rounds where a variant item was not measured show the nearest measured round's value.",
      fr: "Tableau des valeurs par élément de variante pour une catégorie sélectionnée, avec les indicateurs en lignes et les éléments de variante en colonnes, un groupe de colonnes par point temporel. Changez de catégorie à l'aide du sélecteur de réplicant. Les points temporels où un élément de variante n'a pas été mesuré affichent la valeur du point mesuré le plus proche.",
      pt: "Tabela dos valores por item de variante para uma categoria selecionada, com os indicadores em linhas e os itens de variante em colunas, um grupo de colunas por momento temporal. Mude de categoria através do seletor de replicante. Os momentos temporais em que um item de variante não foi medido mostram o valor do momento medido mais próximo.",
    },
    allowedFilters: [
      "hfa_category",
      "hfa_sub_category",
      "hfa_indicator",
      "hfa_variant_item",
      "time_point",
    ],
    createDefaultVisualizationOnInstall: null,
    config: {
      d: {
        type: "table",
        valuesDisDisplayOpt: "col",
        disaggregateBy: [
          { disOpt: "hfa_category", disDisplayOpt: "replicant" },
          { disOpt: "hfa_indicator", disDisplayOpt: "row" },
          { disOpt: "hfa_variant_item", disDisplayOpt: "col" },
          { disOpt: "time_point", disDisplayOpt: "colGroup" },
        ],
        // "" per DOC_VIZPRESET_STANDARDS (tables: user selects) — the app
        // auto-selects the first valid category until they do.
        selectedReplicantValue: "",
        filterBy: [],
      },
      s: {
        ...CF_80_70,
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "HFA indicators by variant item",
          fr: "Indicateurs HFA par élément de variante",
          pt: "Indicadores HFA por item de variante",
        },
        captionRelFontSize: null,
        subCaption: null,
        subCaptionRelFontSize: null,
        footnote: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
];

export const metric: MetricDefinitionGithub = {
  id: "m10-03-02",
  hide: false,
  resultsObjectId: "M10_hfa_results_variants_carried.csv",
  valueProps: ["value"],
  valueFunc: "identity",
  valueLabelReplacements: {},
  label: {
    en: "HFA indicators by variant",
    fr: "Indicateurs HFA par variante",
    pt: "Indicadores HFA por variante",
  },
  requiredDisaggregationOptions: [
    "hfa_indicator",
    "hfa_variant_item",
    "time_point",
  ],
  formatAs: "indicator",
  postAggregationExpression: {
    ingredientValues: [
      { prop: "sum_val", func: "SUM" },
      { prop: "avg_num", func: "SUM" },
      { prop: "avg_weight", func: "SUM" },
    ],
    // Bare division only: the PO query evaluator auto-wraps "/column" in a
    // NULLIF guard; writing NULLIF here would be double-wrapped into bad SQL
    expression: "value = COALESCE(sum_val, avg_num / avg_weight)",
  },
  aiDescription: {
    summary: {
      en: "Same as the HFA indicators by variant metric, but with gaps filled: when a variant item was not measured in a survey round, values are carried over from the nearest round where it was measured (previous round preferred; the next round only when no earlier one exists). Each value is the parent indicator's aggregation applied to one variant item's numerator.",
      fr: "Identique à la métrique des indicateurs HFA par variante, mais avec comblement des lacunes : lorsqu'un élément de variante n'a pas été mesuré lors d'un cycle d'enquête, les valeurs sont reportées depuis le cycle le plus proche où il a été mesuré (cycle précédent en priorité ; cycle suivant uniquement s'il n'existe pas de cycle antérieur). Chaque valeur est l'agrégation de l'indicateur parent appliquée au numérateur d'un élément de variante.",
      pt: "Igual à métrica dos indicadores HFA por variante, mas com preenchimento de lacunas: quando um item de variante não foi medido numa ronda do inquérito, os valores são transportados da ronda mais próxima em que foi medido (ronda anterior de preferência; a ronda seguinte apenas quando não existe uma anterior). Cada valor é a agregação do indicador principal aplicada ao numerador de um item de variante.",
    },
    methodology: {
      en: "Each (indicator, variant item) pair has its own authored numerator R code sharing the parent's filter, type, and aggregation. Rounds where a pair has no observed values at all are filled by duplicating the facility rows of the donor round, so filled rounds exactly reproduce the donor round's aggregates at every administrative level. Values do not sum to the parent's overall value: items may overlap or be non-exhaustive.",
      fr: "Chaque paire (indicateur, élément de variante) a son propre code R de numérateur partageant le filtre, le type et l'agrégation du parent. Les cycles où une paire n'a aucune valeur observée sont comblés en dupliquant les lignes d'établissements du cycle donneur ; les cycles comblés reproduisent donc exactement les agrégats du cycle donneur à tous les niveaux administratifs. Les valeurs ne s'additionnent pas à la valeur globale du parent : les éléments peuvent se chevaucher ou être non exhaustifs.",
      pt: "Cada par (indicador, item de variante) tem o seu próprio código R de numerador partilhando o filtro, o tipo e a agregação do indicador principal. As rondas em que um par não tem quaisquer valores observados são preenchidas duplicando as linhas de unidades sanitárias da ronda dadora, pelo que as rondas preenchidas reproduzem exatamente os agregados da ronda dadora em todos os níveis administrativos. Os valores não somam para o valor global do indicador principal: os itens podem sobrepor-se ou ser não exaustivos.",
    },
    interpretation: {
      en: "Carried values are not new measurements — a flat segment across rounds may simply mean the item's question was not asked in between. Compare with the 'Observed only' variant of this metric to see which rounds contain real data. hfa_indicator carries the PARENT indicator id; hfa_variant_item carries the item id.",
      fr: "Les valeurs reportées ne sont pas de nouvelles mesures — un segment plat entre cycles peut simplement signifier que la question de l'élément n'a pas été posée. Comparez avec la variante « Valeurs observées uniquement » de cette métrique pour identifier les cycles contenant des données réelles. hfa_indicator porte l'identifiant de l'indicateur PARENT ; hfa_variant_item porte l'identifiant de l'élément.",
      pt: "Os valores transportados não são novas medições — um segmento plano entre rondas pode significar simplesmente que a pergunta do item não foi feita. Compare com a variante «Apenas valores observados» desta métrica para identificar as rondas com dados reais. hfa_indicator contém o identificador do indicador PRINCIPAL; hfa_variant_item contém o identificador do item.",
    },
    typicalRange: {
      en: "Mirrors the parent indicator: binary avg indicators range 0-1 (percentage), others vary based on what is being measured.",
      fr: "Reflète l'indicateur parent: les indicateurs binaires moyens varient de 0 à 1 (pourcentage), les autres varient selon ce qui est mesuré.",
      pt: "Reflete o indicador principal: os indicadores binários por média variam de 0 a 1 (percentagem), os restantes variam consoante o que está a ser medido.",
    },
    caveats: {
      en: "Indicators from different variant groups have disjoint item sets, so an unscoped indicator × item cross is block-sparse — scope by category or filter to one variant group. Values in rounds where the item was not measured are copies of the nearest measured round: do not interpret them as evidence of stability, and do not compute change between a filled round and its donor round (the difference is zero by construction).",
      fr: "Les indicateurs de groupes de variantes différents ont des ensembles d'éléments disjoints; un croisement non restreint indicateur × élément est donc épars par blocs — restreindre par catégorie ou filtrer sur un seul groupe de variantes. Les valeurs des cycles où l'élément n'a pas été mesuré sont des copies du cycle mesuré le plus proche : ne les interprétez pas comme une preuve de stabilité et ne calculez pas d'évolution entre un cycle comblé et son cycle donneur (la différence est nulle par construction).",
      pt: "Os indicadores de grupos de variantes diferentes têm conjuntos de itens disjuntos, pelo que um cruzamento não restringido indicador × item é esparso por blocos — restringir por categoria ou filtrar para um único grupo de variantes. Os valores das rondas em que o item não foi medido são cópias da ronda medida mais próxima: não os interprete como prova de estabilidade e não calcule variações entre uma ronda preenchida e a sua ronda dadora (a diferença é zero por construção).",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by hfa_indicator, hfa_variant_item, and time_point (all required). The headline view is indicators as rows × variant items as columns, scoped to one category or variant group.",
      fr: "Toujours désagréger par hfa_indicator, hfa_variant_item et time_point (tous requis). La vue principale est indicateurs en lignes × éléments de variante en colonnes, restreinte à une catégorie ou un groupe de variantes.",
      pt: "Desagregar sempre por hfa_indicator, hfa_variant_item e time_point (todos obrigatórios). A vista principal é indicadores em linhas × itens de variante em colunas, restringida a uma categoria ou grupo de variantes.",
    },
  },
  variantLabel: {
    en: "With carry-forward",
    fr: "Avec valeurs reportées",
    pt: "Com valores transportados",
  },
  importantNotes: null,
  catalogExpressionEvaluation: null,
  vizPresets,
};
