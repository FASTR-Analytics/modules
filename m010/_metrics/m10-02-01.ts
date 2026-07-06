import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { CF_80_70 } from "../../.validation/cf_presets.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "hfa-dont-know-rate-table",
    label: {
      en: "Don't-know rate by indicator and time point",
      fr: "Taux de « Ne sait pas » par indicateur et point temporel",
      pt: "Taxa de \"Não sabe\" por indicador e momento temporal",
    },
    description: {
      en: "Table of the share of facilities answering \"Don't know\" for each HFA indicator, with indicators as rows and time points as columns.",
      fr: "Tableau de la part des établissements ayant répondu « Ne sait pas » pour chaque indicateur HFA, avec les indicateurs en lignes et les points temporels en colonnes.",
      pt: "Tabela da percentagem de unidades sanitárias que responderam \"Não sabe\" para cada indicador HFA, com os indicadores em linhas e os momentos temporais em colunas.",
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
          { disOpt: "hfa_indicator", disDisplayOpt: "row" },
          { disOpt: "time_point", disDisplayOpt: "col" },
        ],
        filterBy: [],
      },
      s: {
        ...CF_80_70,
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "HFA don't-know rates",
          fr: "Taux de « Ne sait pas » HFA",
          pt: "Taxas de \"Não sabe\" HFA",
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
  id: "m10-02-01",
  hide: false,
  resultsObjectId: "M10_hfa_response_status.csv",
  valueProps: ["value"],
  valueFunc: "identity",
  valueLabelReplacements: {},
  label: {
    en: "HFA don't-know rate",
    fr: "Taux de « Ne sait pas » HFA",
    pt: "Taxa de \"Não sabe\" HFA",
  },
  requiredDisaggregationOptions: ["hfa_indicator", "time_point"],
  formatAs: "percent",
  postAggregationExpression: {
    ingredientValues: [
      { prop: "dk_num", func: "SUM" },
      { prop: "resp_weight", func: "SUM" },
    ],
    // Bare division only: the PO query evaluator auto-wraps "/column" in a
    // NULLIF guard; writing NULLIF here would be double-wrapped into bad SQL
    expression: "value = dk_num / resp_weight",
  },
  aiDescription: {
    summary: {
      en: "Share of facilities that answered \"Don't know\" (or the numeric don't-know code) on the survey questions behind each HFA indicator, among facilities the indicator applies to.",
      fr: "Part des établissements ayant répondu « Ne sait pas » (ou le code numérique correspondant) aux questions d'enquête derrière chaque indicateur HFA, parmi les établissements concernés par l'indicateur.",
      pt: "Percentagem de unidades sanitárias que responderam \"Não sabe\" (ou o código numérico correspondente) às perguntas do inquérito subjacentes a cada indicador HFA, entre as unidades às quais o indicador se aplica.",
    },
    methodology: {
      en: "A facility counts as don't-know for an indicator when any survey variable the indicator's formula depends on holds a don't-know code (-99 for select questions, -999999 for numeric questions). The denominator is all facilities the indicator applies to (facilities excluded by the indicator's filter are not counted). Weighted when sampling weights are enabled.",
      fr: "Un établissement compte comme « Ne sait pas » pour un indicateur lorsqu'une variable d'enquête dont dépend la formule de l'indicateur contient un code « Ne sait pas » (-99 pour les questions à choix, -999999 pour les questions numériques). Le dénominateur est l'ensemble des établissements concernés par l'indicateur (les établissements exclus par le filtre de l'indicateur ne sont pas comptés). Pondéré lorsque les pondérations d'échantillonnage sont activées.",
      pt: "Uma unidade sanitária conta como \"Não sabe\" para um indicador quando qualquer variável do inquérito da qual a fórmula do indicador depende contém um código de \"Não sabe\" (-99 para perguntas de escolha, -999999 para perguntas numéricas). O denominador é o conjunto das unidades às quais o indicador se aplica (as unidades excluídas pelo filtro do indicador não são contadas). Ponderado quando as ponderações da amostra estão ativadas.",
    },
    interpretation: {
      en: "High don't-know rates flag questions respondents could not answer reliably — a data-quality and supervision signal. They also show how sensitive the indicator's value is to the don't-know treatment parameter.",
      fr: "Des taux élevés de « Ne sait pas » signalent des questions auxquelles les répondants n'ont pas pu répondre de manière fiable — un signal de qualité des données et de supervision. Ils montrent aussi la sensibilité de la valeur de l'indicateur au paramètre de traitement des « Ne sait pas ».",
      pt: "Taxas elevadas de \"Não sabe\" assinalam perguntas às quais os inquiridos não conseguiram responder de forma fiável — um sinal de qualidade dos dados e de supervisão. Mostram também a sensibilidade do valor do indicador ao parâmetro de tratamento dos \"Não sabe\".",
    },
    typicalRange: {
      en: "0-1 (share of facilities). Usually well below 0.1; higher values warrant scrutiny of the underlying question.",
      fr: "0-1 (part des établissements). Généralement bien en dessous de 0,1 ; des valeurs plus élevées justifient un examen de la question sous-jacente.",
      pt: "0-1 (percentagem de unidades). Geralmente bem abaixo de 0,1; valores mais elevados justificam a análise da pergunta subjacente.",
    },
    caveats: {
      en: "Rates are attributed to indicators, not questions: one don't-know answer can raise the rate of every indicator that uses that question. Data staged before the select_multiple don't-know fix reads as \"No\" rather than don't-know until re-imported.",
      fr: "Les taux sont attribués aux indicateurs, pas aux questions : une seule réponse « Ne sait pas » peut augmenter le taux de tous les indicateurs qui utilisent cette question. Les données importées avant la correction des « Ne sait pas » des questions à choix multiples apparaissent comme « Non » tant qu'elles ne sont pas réimportées.",
      pt: "As taxas são atribuídas aos indicadores, não às perguntas: uma única resposta \"Não sabe\" pode aumentar a taxa de todos os indicadores que utilizam essa pergunta. Os dados carregados antes da correção dos \"Não sabe\" das perguntas de escolha múltipla aparecem como \"Não\" até serem reimportados.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by hfa_indicator and time_point (both required). Disaggregate by admin_area or facility_type to locate where don't-know answers concentrate.",
      fr: "Toujours désagréger par hfa_indicator et time_point (tous deux requis). Désagréger par admin_area ou facility_type pour localiser où se concentrent les réponses « Ne sait pas ».",
      pt: "Desagregar sempre por hfa_indicator e time_point (ambos obrigatórios). Desagregar por admin_area ou facility_type para localizar onde se concentram as respostas \"Não sabe\".",
    },
  },
  variantLabel: null,
  importantNotes: null,
  vizPresets,
};
