import type { ModuleParameter } from "../.validation/module_definition_validator.ts";

export const parameters: ModuleParameter[] = [
  {
    description: "Count variable to use for modeling",
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
    description: "Count variable to use for visualization",
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
    description: "Run district-level model (admin_area_3)",
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
    description: "Run admin_area_4 analysis",
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
    description: "Threshold for MAD-based control limits",
    replacementString: "MADS_THRESHOLD",
    input: {
      inputType: "number",
      defaultValue: "1.5",
    },
  },
  {
    description: "Smoothing window (k)",
    replacementString: "SMOOTH_K",
    input: {
      inputType: "number",
      defaultValue: "7",
    },
  },
  {
    description: "Dip threshold (proportion of expected)",
    replacementString: "DIP_THRESHOLD",
    input: {
      inputType: "number",
      defaultValue: "0.9",
    },
  },
  {
    description: "Difference percent threshold for visualization",
    replacementString: "DIFFPERCENT",
    input: {
      inputType: "number",
      defaultValue: "10",
    },
  },
];
