import type { ModuleParameter } from "../.validation/_module_definition_github.ts";

export const parameters: ModuleParameter[] = [
  {
    description: "Count value to use",
    replacementString: "SELECTED_COUNT_VARIABLE",
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
    description: "Level to calculate coverage for",
    replacementString: "ANALYSIS_LEVEL",
    input: {
      inputType: "select",
      options: [
        {
          value: "NATIONAL_ONLY",
          label: "National only",
        },
        {
          value: "NATIONAL_PLUS_AA2",
          label: "National and admin area 2",
        },
        {
          value: "NATIONAL_PLUS_AA2_AA3",
          label: "National, admin area 2, and admin area 3",
        },
      ],
      valueType: "string",
      defaultValue: "NATIONAL_PLUS_AA2",
    },
  },
  {
    description: "Pregnancy loss rate",
    replacementString: "PREGNANCY_LOSS_RATE",
    input: {
      inputType: "number",
      defaultValue: "0.03",
    },
  },
  {
    description: "Twin rate",
    replacementString: "TWIN_RATE",
    input: {
      inputType: "number",
      defaultValue: "0.015",
    },
  },
  {
    description: "Stillbirth rate",
    replacementString: "STILLBIRTH_RATE",
    input: {
      inputType: "number",
      defaultValue: "0.02",
    },
  },
  {
    description: "Neonatal mortality rate",
    replacementString: "P1_NMR",
    input: {
      inputType: "number",
      defaultValue: "0.039",
    },
  },
  {
    description: "Postneonatal mortality rate",
    replacementString: "P2_PNMR",
    input: {
      inputType: "number",
      defaultValue: "0.028",
    },
  },
  {
    description: "Infant mortality rate",
    replacementString: "INFANT_MORTALITY_RATE",
    input: {
      inputType: "number",
      defaultValue: "0.067",
    },
  },
  {
    description: "Under 5 mortality rate",
    replacementString: "UNDER5_MORTALITY_RATE",
    input: {
      inputType: "number",
      defaultValue: "0.103",
    },
  },
];
