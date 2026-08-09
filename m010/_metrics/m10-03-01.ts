import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { CF_80_70 } from "../../.validation/cf_presets.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "hfa-variant-indicator-item-table",
    label: {
      en: "Indicators by variant item (per category)",
      fr: "Indicateurs par élément de variante (par catégorie)",
      pt: "Indicadores por item de variante (por categoria)",
    },
    description: {
      en: "Table of per-item variant values for a selected category, with indicators as rows and variant items as columns, one column group per time point. Switch the category using the replicant selector.",
      fr: "Tableau des valeurs par élément de variante pour une catégorie sélectionnée, avec les indicateurs en lignes et les éléments de variante en colonnes, un groupe de colonnes par point temporel. Changez de catégorie à l'aide du sélecteur de réplicant.",
      pt: "Tabela dos valores por item de variante para uma categoria selecionada, com os indicadores em linhas e os itens de variante em colunas, um grupo de colunas por momento temporal. Mude de categoria através do seletor de replicante.",
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
  id: "m10-03-01",
  hide: false,
  resultsObjectId: "M10_hfa_results_variants.csv",
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
      en: "Per-item breakdown of Health Facility Assessment indicators authored with a variant group: each value is the parent indicator's aggregation applied to one variant item's numerator (e.g. vaccination through campaigns vs routine services).",
      fr: "Ventilation par élément des indicateurs d'évaluation des établissements de santé définis avec un groupe de variantes: chaque valeur est l'agrégation de l'indicateur parent appliquée au numérateur d'un élément de variante (par exemple vaccination par campagnes vs services de routine).",
      pt: "Desagregação por item dos indicadores de Avaliação de Unidades Sanitárias definidos com um grupo de variantes: cada valor é a agregação do indicador principal aplicada ao numerador de um item de variante (por exemplo vacinação por campanhas vs serviços de rotina).",
    },
    methodology: {
      en: "Each (indicator, variant item) pair has its own authored numerator R code sharing the parent's filter, type, and aggregation. Per-item denominators can legitimately differ because each item's expression gates on its own missingness. Values do not sum to the parent's overall value: the overall is authored separately and items may overlap or be non-exhaustive.",
      fr: "Chaque paire (indicateur, élément de variante) a son propre code R de numérateur partageant le filtre, le type et l'agrégation du parent. Les dénominateurs par élément peuvent légitimement différer. Les valeurs ne s'additionnent pas à la valeur globale du parent: celle-ci est définie séparément et les éléments peuvent se chevaucher ou être non exhaustifs.",
      pt: "Cada par (indicador, item de variante) tem o seu próprio código R de numerador partilhando o filtro, o tipo e a agregação do indicador principal. Os denominadores por item podem legitimamente diferir. Os valores não somam para o valor global do indicador principal: este é definido separadamente e os itens podem sobrepor-se ou ser não exaustivos.",
    },
    interpretation: {
      en: "Use to compare how an indicator's result splits across its variant items. hfa_indicator carries the PARENT indicator id; hfa_variant_item carries the item id. Only indicators assigned a variant group appear here.",
      fr: "À utiliser pour comparer la répartition du résultat d'un indicateur entre ses éléments de variante. hfa_indicator porte l'identifiant de l'indicateur PARENT; hfa_variant_item porte l'identifiant de l'élément. Seuls les indicateurs assignés à un groupe de variantes apparaissent ici.",
      pt: "Utilizar para comparar como o resultado de um indicador se reparte pelos seus itens de variante. hfa_indicator contém o identificador do indicador PRINCIPAL; hfa_variant_item contém o identificador do item. Apenas os indicadores atribuídos a um grupo de variantes aparecem aqui.",
    },
    typicalRange: {
      en: "Mirrors the parent indicator: binary avg indicators range 0-1 (percentage), others vary based on what is being measured.",
      fr: "Reflète l'indicateur parent: les indicateurs binaires moyens varient de 0 à 1 (pourcentage), les autres varient selon ce qui est mesuré.",
      pt: "Reflete o indicador principal: os indicadores binários por média variam de 0 a 1 (percentagem), os restantes variam consoante o que está a ser medido.",
    },
    caveats: {
      en: "Indicators from different variant groups have disjoint item sets, so an unscoped indicator × item cross is block-sparse — scope by category or filter to one variant group. Item values are not comparable to the parent's overall value row-for-row.",
      fr: "Les indicateurs de groupes de variantes différents ont des ensembles d'éléments disjoints; un croisement non restreint indicateur × élément est donc épars par blocs — restreindre par catégorie ou filtrer sur un seul groupe de variantes. Les valeurs des éléments ne sont pas comparables ligne à ligne à la valeur globale du parent.",
      pt: "Os indicadores de grupos de variantes diferentes têm conjuntos de itens disjuntos, pelo que um cruzamento não restringido indicador × item é esparso por blocos — restringir por categoria ou filtrar para um único grupo de variantes. Os valores dos itens não são comparáveis linha a linha ao valor global do indicador principal.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by hfa_indicator, hfa_variant_item, and time_point (all required). The headline view is indicators as rows × variant items as columns, scoped to one category or variant group.",
      fr: "Toujours désagréger par hfa_indicator, hfa_variant_item et time_point (tous requis). La vue principale est indicateurs en lignes × éléments de variante en colonnes, restreinte à une catégorie ou un groupe de variantes.",
      pt: "Desagregar sempre por hfa_indicator, hfa_variant_item e time_point (todos obrigatórios). A vista principal é indicadores em linhas × itens de variante em colunas, restringida a uma categoria ou grupo de variantes.",
    },
  },
  variantLabel: {
    en: "Observed only",
    fr: "Valeurs observées uniquement",
    pt: "Apenas valores observados",
  },
  importantNotes: null,
  vizPresets,
};
