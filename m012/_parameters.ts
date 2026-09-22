import type { ModuleParameter } from "../.validation/_module_definition_github.ts";

// One adjustment basis per results package (PLAN_1a §1.5). Comparing
// adjustment scenarios means regenerating with another basis, which is the
// system's normal lifecycle.
export const parameters: ModuleParameter[] = [
  {
    description: {
      en: "Count variable to use",
      fr: "Variable de comptage à utiliser",
      pt: "Variável de contagem a utilizar",
    },
    replacementString: "SELECTEDCOUNT",
    input: {
      inputType: "select",
      options: [
        { value: "count_final_none", label: "Count (unadjusted)" },
        { value: "count_final_outliers", label: "Count (adjusted for outliers)" },
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
];
