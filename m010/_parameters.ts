import type { ModuleParameter } from "../.validation/_module_definition_github.ts";

export const parameters: ModuleParameter[] = [
  {
    replacementString: "USE_SAMPLE_WEIGHTS",
    description: {
      en: "Use sample weights for calculations",
      fr: "Utiliser les pondérations de l'échantillon pour les calculs",
      pt: "Utilizar ponderações da amostra para os cálculos",
    },
    input: {
      inputType: "boolean",
      defaultValue: "FALSE",
    },
  },
  {
    replacementString: "DONT_KNOW_TREATMENT",
    description: {
      en: 'How to treat "Don\'t know" (-99) responses in binary indicators',
      fr: "Comment traiter les réponses « Ne sait pas » (-99) dans les indicateurs binaires",
      pt: 'Como tratar as respostas "Não sabe" (-99) nos indicadores binários',
    },
    input: {
      inputType: "select",
      valueType: "string",
      options: [
        {
          value: "missing",
          label: "Treat as missing (facility excluded from the indicator)",
        },
        {
          value: "no",
          label:
            'Treat as "No" (facility counted as not meeting the indicator)',
        },
      ],
      defaultValue: "missing",
    },
  },
  {
    replacementString: "STOP_IF_INDICATOR_FAILS",
    description: {
      en: "Stop script if any single indicator fails (no R code or errors). If FALSE, skip failing indicators and only fail if no results are generated.",
      fr: "Arrêter le script si un indicateur échoue (pas de code R ou erreurs). Si FALSE, ignorer les indicateurs qui échouent et échouer uniquement si aucun résultat n'est généré.",
      pt: "Parar o script se algum indicador falhar (sem código R ou erros). Se FALSE, ignorar os indicadores que falham e falhar apenas se não forem gerados resultados.",
    },
    input: {
      inputType: "boolean",
      defaultValue: "TRUE",
    },
  },
];
