import type { ModuleParameter } from "../.validation/_module_definition_github.ts";

export const parameters: ModuleParameter[] = [
  {
    description: {
      en: "Count value to use",
      fr: "Valeur de comptage à utiliser",
    },
    replacementString: "SELECTED_COUNT_VARIABLE",
    input: {
      inputType: "select",
      valueType: "string",
      defaultValue: "count_final_none",
      options: [
        { value: "count_final_none", label: "Count (unadjusted)" },
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
    },
  },
  {
    description: {
      en: "Skip indicators with missing data instead of failing",
      fr: "Ignorer les indicateurs avec données manquantes au lieu d'échouer",
    },
    replacementString: "__SKIP_MISSING_INDICATORS_VALUE__",
    input: {
      inputType: "select",
      valueType: "number",
      defaultValue: "TRUE",
      options: [
        { value: "TRUE", label: "Skip missing (recommended)" },
        { value: "FALSE", label: "Fail on missing" },
      ],
    },
  },
];
