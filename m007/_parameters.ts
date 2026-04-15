import type { ModuleParameter } from "../.validation/module_definition_validator.ts";

export const parameters: ModuleParameter[] = [
  {
    description: "Count value to use",
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
    replacementString: "BIRTHS_PCT",
    description:
      "Quarterly births as proportion of total population (annual rate x 0.25)",
    input: {
      inputType: "number",
      defaultValue: "0.01",
    },
  },
  {
    replacementString: "WOMEN_15_49_PCT",
    description:
      "Quarterly women 15-49 as proportion of total population (annual rate x 0.25)",
    input: {
      inputType: "number",
      defaultValue: "0.055",
    },
  },
  {
    replacementString: "INTERPOLATE_POPULATION",
    description: "Use smoothed (interpolated) population or raw annual values",
    input: {
      inputType: "select",
      valueType: "number",
      defaultValue: "FALSE",
      options: [
        { value: "FALSE", label: "Raw (annual step)" },
        { value: "TRUE", label: "Smoothed (interpolated)" },
      ],
    },
  },
];
