import type { ModuleParameter } from "../.validation/_module_definition_github.ts";

export const parameters: ModuleParameter[] = [
  {
    description: {
      en: "Count variable to use for modeling",
      fr: "Variable de comptage à utiliser pour la modélisation",
      pt: "Variável de contagem a utilizar para a modelação",
    },
    replacementString: "SELECTEDCOUNT",
    input: {
      inputType: "select",
      options: [
        {
          value: "count_final_none",
          label: "Count (unadjusted)",
        },
        {
          value: "count_final_outliers",
          label: "Count (adjusted for outliers)",
        },
        {
          value: "count_final_completeness",
          label: "Count (adjusted for completeness)",
        },
        {
          value: "count_final_both",
          label: "Count (adjusted for outliers and completeness)",
        },
      ],
      valueType: "string",
      defaultValue: "count_final_outliers",
    },
  },
  {
    description: {
      en: "Count variable to use for visualization",
      fr: "Variable de comptage à utiliser pour la visualisation",
      pt: "Variável de contagem a utilizar para a visualização",
    },
    replacementString: "VISUALIZATIONCOUNT",
    input: {
      inputType: "select",
      options: [
        {
          value: "count_final_none",
          label: "Count (unadjusted)",
        },
        {
          value: "count_final_outliers",
          label: "Count (adjusted for outliers)",
        },
        {
          value: "count_final_completeness",
          label: "Count (adjusted for completeness)",
        },
        {
          value: "count_final_both",
          label: "Count (adjusted for outliers and completeness)",
        },
      ],
      valueType: "string",
      defaultValue: "count_final_outliers",
    },
  },
  {
    description: {
      en: "Run district-level model (admin_area_3)",
      fr: "Exécuter le modèle au niveau du district (admin_area_3)",
      pt: "Executar o modelo ao nível do distrito (admin_area_3)",
    },
    replacementString: "RUN_DISTRICT_MODEL",
    input: {
      inputType: "select",
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
      valueType: "number",
      defaultValue: "FALSE",
    },
  },
  {
    description: {
      en: "Run admin_area_4 analysis",
      fr: "Exécuter l'analyse admin_area_4",
      pt: "Executar a análise admin_area_4",
    },
    replacementString: "RUN_ADMIN_AREA_4_ANALYSIS",
    input: {
      inputType: "select",
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
      valueType: "number",
      defaultValue: "FALSE",
    },
  },
  {
    description: {
      en: "Threshold for MAD-based control limits",
      fr: "Seuil pour les limites de contrôle basées sur MAD",
      pt: "Limiar para os limites de controlo baseados no MAD",
    },
    replacementString: "MADS_THRESHOLD",
    input: {
      inputType: "number",
      defaultValue: "1.5",
    },
  },
  {
    description: {
      en: "Smoothing window (k)",
      fr: "Fenêtre de lissage (k)",
      pt: "Janela de suavização (k)",
    },
    replacementString: "SMOOTH_K",
    input: {
      inputType: "number",
      defaultValue: "7",
    },
  },
  {
    description: {
      en: "Dip threshold (proportion of expected)",
      fr: "Seuil de baisse (proportion de l'attendu)",
      pt: "Limiar de queda (proporção do esperado)",
    },
    replacementString: "DIP_THRESHOLD",
    input: {
      inputType: "number",
      defaultValue: "0.9",
    },
  },
  {
    description: {
      en: "Difference percent threshold for visualization",
      fr: "Seuil de pourcentage de différence pour la visualisation",
      pt: "Limiar de percentagem de diferença para a visualização",
    },
    replacementString: "DIFFPERCENT",
    input: {
      inputType: "number",
      defaultValue: "10",
    },
  },
];
