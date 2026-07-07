import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [];

export const metric: MetricDefinitionGithub = {
  id: "m10-02-02",
  hide: false,
  resultsObjectId: "M10_hfa_response_status.csv",
  valueProps: ["value"],
  valueFunc: "identity",
  valueLabelReplacements: {},
  label: {
    en: "HFA missing-data rate",
    fr: "Taux de données manquantes HFA",
    pt: "Taxa de dados em falta HFA",
  },
  requiredDisaggregationOptions: ["hfa_indicator", "time_point"],
  formatAs: "percent",
  postAggregationExpression: {
    ingredientValues: [
      { prop: "missing_num", func: "SUM" },
      { prop: "resp_weight", func: "SUM" },
    ],
    // Bare division only: the PO query evaluator auto-wraps "/column" in a
    // NULLIF guard; writing NULLIF here would be double-wrapped into bad SQL
    expression: "value = missing_num / resp_weight",
  },
  aiDescription: {
    summary: {
      en: "Share of facilities with no answer recorded on the survey questions behind each HFA indicator, among facilities the indicator applies to.",
      fr: "Part des établissements sans réponse enregistrée aux questions d'enquête derrière chaque indicateur HFA, parmi les établissements concernés par l'indicateur.",
      pt: "Percentagem de unidades sanitárias sem resposta registada às perguntas do inquérito subjacentes a cada indicador HFA, entre as unidades às quais o indicador se aplica.",
    },
    methodology: {
      en: "A facility counts as missing for an indicator when any survey variable the indicator's formula depends on is empty (skipped by form logic or not recorded) and no don't-know code is present. The denominator is all facilities the indicator applies to. Weighted when sampling weights are enabled.",
      fr: "Un établissement compte comme manquant pour un indicateur lorsqu'une variable d'enquête dont dépend la formule de l'indicateur est vide (sautée par la logique du formulaire ou non enregistrée) et qu'aucun code « Ne sait pas » n'est présent. Le dénominateur est l'ensemble des établissements concernés par l'indicateur. Pondéré lorsque les pondérations d'échantillonnage sont activées.",
      pt: 'Uma unidade sanitária conta como em falta para um indicador quando qualquer variável do inquérito da qual a fórmula do indicador depende está vazia (saltada pela lógica do formulário ou não registada) e não há código de "Não sabe". O denominador é o conjunto das unidades às quais o indicador se aplica. Ponderado quando as ponderações da amostra estão ativadas.',
    },
    interpretation: {
      en: "High missing rates shrink the effective denominator of the indicator itself and can bias it if missingness is not random. Compare with the don't-know rate to separate skipped questions from uncertain respondents.",
      fr: "Des taux élevés de données manquantes réduisent le dénominateur effectif de l'indicateur lui-même et peuvent le biaiser si les données manquantes ne sont pas aléatoires. Comparer avec le taux de « Ne sait pas » pour distinguer les questions sautées des répondants incertains.",
      pt: 'Taxas elevadas de dados em falta reduzem o denominador efetivo do próprio indicador e podem enviesá-lo se a falta de dados não for aleatória. Comparar com a taxa de "Não sabe" para separar perguntas saltadas de inquiridos incertos.',
    },
    typicalRange: {
      en: "0-1 (share of facilities). Values above 0.2 suggest form-logic or fieldwork problems for the underlying questions.",
      fr: "0-1 (part des établissements). Des valeurs supérieures à 0,2 suggèrent des problèmes de logique de formulaire ou de travail de terrain pour les questions sous-jacentes.",
      pt: "0-1 (percentagem de unidades). Valores acima de 0,2 sugerem problemas de lógica do formulário ou de trabalho de campo nas perguntas subjacentes.",
    },
    caveats: {
      en: 'Rates are attributed to indicators, not questions: one unanswered question can raise the rate of every indicator that uses it. Data staged before the select_multiple missingness fix reads unanswered multi-select questions as "No" (never missing) until re-imported.',
      fr: "Les taux sont attribués aux indicateurs, pas aux questions : une seule question sans réponse peut augmenter le taux de tous les indicateurs qui l'utilisent. Les données importées avant la correction des questions à choix multiples sans réponse apparaissent comme « Non » (jamais manquantes) tant qu'elles ne sont pas réimportées.",
      pt: 'As taxas são atribuídas aos indicadores, não às perguntas: uma única pergunta sem resposta pode aumentar a taxa de todos os indicadores que a utilizam. Os dados carregados antes da correção das perguntas de escolha múltipla sem resposta aparecem como "Não" (nunca em falta) até serem reimportados.',
    },
    disaggregationGuidance: {
      en: "Always disaggregate by hfa_indicator and time_point (both required). Disaggregate by admin_area or facility_type to locate where missing data concentrates.",
      fr: "Toujours désagréger par hfa_indicator et time_point (tous deux requis). Désagréger par admin_area ou facility_type pour localiser où se concentrent les données manquantes.",
      pt: "Desagregar sempre por hfa_indicator e time_point (ambos obrigatórios). Desagregar por admin_area ou facility_type para localizar onde se concentram os dados em falta.",
    },
  },
  variantLabel: null,
  importantNotes: null,
  vizPresets,
};
