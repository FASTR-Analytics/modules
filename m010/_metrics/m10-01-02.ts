import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { CF_80_70 } from "../../.validation/cf_presets.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "hfa-category-indicator-table-carried",
    label: {
      en: "Indicators by time point, with carry-forward (per service category)",
      fr: "Indicateurs par point temporel, avec valeurs reportées (par catégorie de service)",
      pt: "Indicadores por momento temporal, com valores transportados (por categoria de serviço)",
    },
    description: {
      en: "Table of HFA indicator values for a selected service category, with indicators as rows and time points as columns. Switch the category using the replicant selector. Rounds where an indicator was not measured show the nearest measured round's value.",
      fr: "Tableau des valeurs des indicateurs HFA pour une catégorie de service sélectionnée, avec les indicateurs en lignes et les points temporels en colonnes. Changez de catégorie à l'aide du sélecteur de réplicant. Les points temporels où un indicateur n'a pas été mesuré affichent la valeur du point mesuré le plus proche.",
      pt: "Tabela dos valores dos indicadores HFA para uma categoria de serviço selecionada, com os indicadores em linhas e os momentos temporais em colunas. Mude de categoria através do seletor de replicante. Os momentos temporais em que um indicador não foi medido mostram o valor do momento medido mais próximo.",
    },
    allowedFilters: [
      "hfa_category",
      "hfa_sub_category",
      "hfa_indicator",
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
          { disOpt: "time_point", disDisplayOpt: "col" },
        ],
        selectedReplicantValue: "financing",
        filterBy: [],
      },
      s: {
        ...CF_80_70,
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "HFA indicators",
          fr: "Indicateurs HFA",
          pt: "Indicadores HFA",
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
  id: "m10-01-02",
  hide: false,
  resultsObjectId: "M10_hfa_results_carried.csv",
  valueProps: ["value"],
  valueFunc: "identity",
  valueLabelReplacements: {},
  label: {
    en: "HFA indicators",
    fr: "Indicateurs HFA",
    pt: "Indicadores HFA",
  },
  requiredDisaggregationOptions: ["hfa_indicator", "time_point"],
  formatAs: "number",
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
      en: "Same as the HFA indicators metric, but with gaps filled: when an indicator was not measured in a survey round, values are carried over from the nearest round where it was measured (previous round preferred; the next round only when no earlier one exists). Aggregation follows each indicator's defined method (sum or average) and type (binary or numeric).",
      fr: "Identique à la métrique des indicateurs HFA, mais avec comblement des lacunes : lorsqu'un indicateur n'a pas été mesuré lors d'un cycle d'enquête, les valeurs sont reportées depuis le cycle le plus proche où il a été mesuré (cycle précédent en priorité ; cycle suivant uniquement s'il n'existe pas de cycle antérieur). L'agrégation suit la méthode définie pour chaque indicateur (somme ou moyenne) et son type (binaire ou numérique).",
      pt: "Igual à métrica dos indicadores HFA, mas com preenchimento de lacunas: quando um indicador não foi medido numa ronda do inquérito, os valores são transportados da ronda mais próxima em que foi medido (ronda anterior de preferência; a ronda seguinte apenas quando não existe uma anterior). A agregação segue o método definido para cada indicador (soma ou média) e o seu tipo (binário ou numérico).",
    },
    methodology: {
      en: "Each indicator is automatically aggregated using its configured method: binary indicators with avg aggregation yield percentages (0-1), binary with sum yield facility counts, numeric with avg yield mean values, numeric with sum yield totals. Rounds where an indicator has no observed values at all are filled by duplicating the facility rows of the donor round, so filled rounds exactly reproduce the donor round's aggregates at every administrative level.",
      fr: "Chaque indicateur est automatiquement agrégé selon sa méthode configurée: les indicateurs binaires avec agrégation moyenne donnent des pourcentages (0-1), binaires avec somme donnent des comptes d'établissements, numériques avec moyenne donnent des valeurs moyennes, numériques avec somme donnent des totaux. Les cycles où un indicateur n'a aucune valeur observée sont comblés en dupliquant les lignes d'établissements du cycle donneur ; les cycles comblés reproduisent donc exactement les agrégats du cycle donneur à tous les niveaux administratifs.",
      pt: "Cada indicador é automaticamente agregado segundo o seu método configurado: os indicadores binários com agregação por média dão percentagens (0-1), binários com soma dão contagens de unidades sanitárias, numéricos com média dão valores médios, numéricos com soma dão totais. As rondas em que um indicador não tem quaisquer valores observados são preenchidas duplicando as linhas de unidades sanitárias da ronda dadora, pelo que as rondas preenchidas reproduzem exatamente os agregados da ronda dadora em todos os níveis administrativos.",
    },
    interpretation: {
      en: "Carried values are not new measurements — a flat segment across rounds may simply mean the question was not asked in between. Compare with the 'Observed only' variant of this metric to see which rounds contain real data. Otherwise interpretation follows the indicator's type and aggregation method.",
      fr: "Les valeurs reportées ne sont pas de nouvelles mesures — un segment plat entre cycles peut simplement signifier que la question n'a pas été posée. Comparez avec la variante « Valeurs observées uniquement » de cette métrique pour identifier les cycles contenant des données réelles. Sinon, l'interprétation suit le type d'indicateur et la méthode d'agrégation.",
      pt: "Os valores transportados não são novas medições — um segmento plano entre rondas pode significar simplesmente que a pergunta não foi feita. Compare com a variante «Apenas valores observados» desta métrica para identificar as rondas com dados reais. De resto, a interpretação segue o tipo de indicador e o método de agregação.",
    },
    typicalRange: {
      en: "Varies by indicator type: binary avg indicators range 0-1 (percentage), others vary based on what is being measured.",
      fr: "Varie selon le type d'indicateur: les indicateurs binaires moyens varient de 0 à 1 (pourcentage), les autres varient selon ce qui est mesuré.",
      pt: "Varia consoante o tipo de indicador: os indicadores binários por média variam de 0 a 1 (percentagem), os restantes variam consoante o que está a ser medido.",
    },
    caveats: {
      en: "HFA data is typically cross-sectional and may not reflect temporal trends. Survey timing and facility sampling affect comparability across assessments. When sampling weights are enabled, sums are weighted population-total estimates (non-integer) and may be unreliable at fine admin-area disaggregation if the weights were designed for national or stratum-level estimates. Values in rounds where the indicator was not measured are copies of the nearest measured round: do not interpret them as evidence of stability, and do not compute change between a filled round and its donor round (the difference is zero by construction).",
      fr: "Les données HFA sont généralement transversales et peuvent ne pas refléter les tendances temporelles. Lorsque les pondérations d'échantillonnage sont activées, les sommes sont des estimations pondérées de totaux de population (non entières) et peuvent être peu fiables à une désagrégation administrative fine si les pondérations ont été conçues pour des estimations nationales ou par strate. Les valeurs des cycles où l'indicateur n'a pas été mesuré sont des copies du cycle mesuré le plus proche : ne les interprétez pas comme une preuve de stabilité et ne calculez pas d'évolution entre un cycle comblé et son cycle donneur (la différence est nulle par construction).",
      pt: "Os dados HFA são geralmente transversais e podem não refletir tendências temporais. O momento do inquérito e a amostragem das unidades sanitárias afetam a comparabilidade entre avaliações. Quando as ponderações da amostra estão ativadas, as somas são estimativas ponderadas de totais populacionais (não inteiras) e podem ser pouco fiáveis numa desagregação administrativa fina se as ponderações tiverem sido concebidas para estimativas nacionais ou ao nível do estrato. Os valores das rondas em que o indicador não foi medido são cópias da ronda medida mais próxima: não os interprete como prova de estabilidade e não calcule variações entre uma ronda preenchida e a sua ronda dadora (a diferença é zero por construção).",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by hfa_indicator and time_point (both required). Use admin_area to compare regional results. Disaggregate by facility_type or facility_ownership to identify disparities.",
      fr: "Toujours désagréger par hfa_indicator et time_point (tous deux requis). Utiliser admin_area pour comparer les résultats régionaux.",
      pt: "Desagregar sempre por hfa_indicator e time_point (ambos obrigatórios). Utilizar admin_area para comparar resultados regionais. Desagregar por facility_type ou facility_ownership para identificar disparidades.",
    },
  },
  variantLabel: {
    en: "With carry-forward",
    fr: "Avec valeurs reportées",
    pt: "Com valores transportados",
  },
  importantNotes: null,
  vizPresets,
};
