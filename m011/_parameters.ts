import type { ModuleParameter } from "../.validation/_module_definition_github.ts";

export const parameters: ModuleParameter[] = [
  {
    replacementString: "SELECTEDCOUNT",
    description: {
      en: "Count variable to model. Completeness-imputed counts are deliberately not offered: this module models reporting gaps itself, so it must see them in the data",
      fr: "Variable de comptage à modéliser. Les comptages imputés pour complétude ne sont volontairement pas proposés : ce module modélise lui-même les lacunes de rapportage et doit donc les voir dans les données",
      pt: "Variável de contagem a modelar. As contagens imputadas por completude não são propositadamente oferecidas: este módulo modela ele próprio as lacunas de reporte, pelo que deve vê-las nos dados",
    },
    input: {
      inputType: "select",
      valueType: "string",
      options: [
        {
          value: "count_final_none",
          label: "Count (unadjusted)",
        },
        {
          value: "count_final_outliers",
          label: "Count (adjusted for outliers)",
        },
      ],
      defaultValue: "count_final_outliers",
    },
  },
  {
    replacementString: "TRIM_FROM_PERIOD",
    description: {
      en: "Trim data from this period (YYYYMM, leave blank for all data)",
      fr: "Tronquer les données à partir de cette période (AAAAMM, laisser vide pour toutes les données)",
      pt: "Truncar os dados a partir deste período (AAAAMM, deixar em branco para todos os dados)",
    },
    input: {
      inputType: "number",
      defaultValue: "",
    },
  },
  {
    replacementString: "FOURIER_K",
    description: {
      en: "Number of Fourier seasonality terms (K). Raise to 3-4 for campaign-style indicators with sharp seasonal spikes (e.g. Vitamin A)",
      fr: "Nombre de termes de saisonnalité de Fourier (K). Augmenter à 3-4 pour les indicateurs de type campagne avec des pics saisonniers marqués (ex. vitamine A)",
      pt: "Número de termos de sazonalidade de Fourier (K). Aumentar para 3-4 para indicadores de tipo campanha com picos sazonais acentuados (ex. vitamina A)",
    },
    input: {
      inputType: "number",
      defaultValue: "2",
    },
  },
  {
    replacementString: "USE_POSTERIOR_DRAWS",
    description: {
      en: "Use joint posterior draws for proper subnational credible intervals (slower, memory-heavy on large countries)",
      fr: "Utiliser des tirages a posteriori conjoints pour des intervalles de crédibilité sous-nationaux corrects (plus lent, gourmand en mémoire pour les grands pays)",
      pt: "Utilizar amostras conjuntas da posteriori para intervalos de credibilidade subnacionais corretos (mais lento, pesado em memória para países grandes)",
    },
    input: {
      inputType: "select",
      valueType: "number",
      options: [
        {
          value: "TRUE",
          label: "Yes",
        },
        {
          value: "FALSE",
          label: "No",
        },
      ],
      defaultValue: "FALSE",
    },
  },
  {
    replacementString: "ZEROS_REAL",
    description: {
      en: "How to treat explicit zeros: all = real reports of zero activity (canonical); detect = pick all/auto/strict automatically from the data's zero signature; auto = real only when the facility submitted another indicator that month; strict = like auto but the other indicator must be > 0 (for systems that auto-fill zeros across all indicators); none = all zeros treated as missing",
      fr: "Traitement des zéros explicites : all = vrais rapports d'activité nulle (canonique) ; detect = choix automatique all/auto/strict selon la signature des zéros dans les données ; auto = réel seulement si l'établissement a soumis un autre indicateur ce mois-là ; strict = comme auto mais l'autre indicateur doit être > 0 (pour les systèmes qui remplissent des zéros sur tous les indicateurs) ; none = tous les zéros traités comme manquants",
      pt: "Tratamento dos zeros explícitos: all = relatórios reais de atividade nula (canónico); detect = escolha automática all/auto/strict segundo a assinatura dos zeros nos dados; auto = real apenas se a unidade submeteu outro indicador nesse mês; strict = como auto mas o outro indicador deve ser > 0 (para sistemas que preenchem zeros em todos os indicadores); none = todos os zeros tratados como em falta",
    },
    input: {
      inputType: "select",
      valueType: "string",
      options: [
        {
          value: "all",
          label: "All zeros are real (canonical)",
        },
        {
          value: "detect",
          label: "Auto-detect from the data",
        },
        {
          value: "auto",
          label: "Cross-indicator detection (corrected)",
        },
        {
          value: "strict",
          label: "Cross-indicator detection (other indicator > 0)",
        },
        {
          value: "none",
          label: "All zeros treated as missing",
        },
      ],
      defaultValue: "all",
    },
  },
  {
    replacementString: "INLA_THREADS",
    description: {
      en: "INLA threading, e.g. 1:1 to force single-threaded model fits — several-fold lower peak memory at some speed cost. Leave blank for INLA's default",
      fr: "Threads INLA, ex. 1:1 pour forcer des ajustements mono-thread — pic mémoire réduit plusieurs fois au prix de la vitesse. Laisser vide pour la valeur par défaut d'INLA",
      pt: "Threads INLA, ex. 1:1 para forçar ajustes mono-thread — pico de memória várias vezes menor à custa de velocidade. Deixar em branco para o padrão do INLA",
    },
    input: {
      inputType: "text",
      defaultValue: "",
    },
  },
  {
    replacementString: "POSTCLOSURE_GRACE",
    description: {
      en: "Presumed-closure grace window in months: grid months more than N months after a facility's last-ever report are dropped as presumed closure (the first N silent months still count as non-reporting). 0 = off. Only relevant where many facilities go permanently dark (e.g. Liberia)",
      fr: "Fenêtre de grâce de fermeture présumée en mois : les mois de grille au-delà de N mois après le dernier rapport d'un établissement sont exclus comme fermeture présumée (les N premiers mois de silence comptent toujours comme non-rapportage). 0 = désactivé. Pertinent seulement là où beaucoup d'établissements cessent définitivement de rapporter (ex. Libéria)",
      pt: "Janela de tolerância de encerramento presumido em meses: os meses da grelha para além de N meses após o último relatório de uma unidade são excluídos como encerramento presumido (os primeiros N meses de silêncio contam como não-reporte). 0 = desativado. Relevante apenas onde muitas unidades deixam de reportar definitivamente (ex. Libéria)",
    },
    input: {
      inputType: "number",
      defaultValue: "0",
    },
  },
  {
    replacementString: "EXCLUDE_PRELAUNCH",
    description: {
      en: "Exclude grid months before a facility's first-ever report, so facilities that opened partway through the period (and indicators rolled out late) are not counted as non-reporting before they existed",
      fr: "Exclure les mois de grille antérieurs au tout premier rapport d'un établissement, afin que les établissements ouverts en cours de période (et les indicateurs déployés tardivement) ne soient pas comptés comme non-rapporteurs avant d'exister",
      pt: "Excluir os meses da grelha anteriores ao primeiro relatório de uma unidade, para que as unidades abertas a meio do período (e os indicadores introduzidos tardiamente) não contem como não-reportantes antes de existirem",
    },
    input: {
      inputType: "select",
      valueType: "number",
      options: [
        {
          value: "TRUE",
          label: "Yes",
        },
        {
          value: "FALSE",
          label: "No",
        },
      ],
      defaultValue: "TRUE",
    },
  },
];
