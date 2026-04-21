import type { ModuleParameter } from "../.validation/_module_definition_github.ts";

export const parameters: ModuleParameter[] = [
  {
    description: "Proportion threshold for outlier detection",
    replacementString: "OUTLIER_PROPORTION_THRESHOLD",
    input: {
      inputType: "number",
      defaultValue: "0.8",
    },
  },
  {
    description: "Minimum count threshold for consideration",
    replacementString: "MINIMUM_COUNT_THRESHOLD",
    input: {
      inputType: "number",
      defaultValue: "100",
    },
  },
  {
    description: "Number of MADs",
    replacementString: "MADS",
    input: {
      inputType: "number",
      defaultValue: "10",
    },
  },
  {
    description: "Indicators subjected to DQA",
    replacementString: "DQA_INDICATORS",
    input: {
      inputType: "select",
      options: [
        {
          value: 'c("anc1", "penta1", "opd")',
          label: "anc1, penta1, opd",
        },
        {
          value: 'c("anc1", "penta1")',
          label: "anc1, penta1",
        },
        {
          value: 'c("anc1", "opd")',
          label: "anc1, opd",
        },
        {
          value: 'c("penta1", "opd")',
          label: "penta1, opd",
        },
        {
          value: 'c("anc1")',
          label: "anc1",
        },
        {
          value: 'c("penta1")',
          label: "penta1",
        },
        {
          value: 'c("opd")',
          label: "opd",
        },
      ],
      valueType: "number",
      defaultValue: 'c("anc1", "penta1", "opd")',
    },
  },
  {
    description: "Consistency pairs used",
    replacementString: "CONSISTENCY_PAIRS_USED",
    input: {
      inputType: "select",
      options: [
        {
          value: 'c("penta", "anc", "delivery", "malaria")',
          label: "penta, anc, delivery, malaria",
        },
        {
          value: 'c("penta", "anc", "delivery")',
          label: "penta, anc, delivery",
        },
        {
          value: 'c("penta", "anc")',
          label: "penta, anc",
        },
        {
          value: 'c("penta", "delivery")',
          label: "penta, delivery",
        },
        {
          value: 'c("anc", "delivery")',
          label: "anc, delivery",
        },
        {
          value: 'c("anc")',
          label: "anc",
        },
        {
          value: 'c("delivery")',
          label: "delivery",
        },
        {
          value: 'c("penta")',
          label: "penta",
        },
        {
          value: 'c("malaria")',
          label: "malaria",
        },
      ],
      valueType: "number",
      defaultValue: 'c("penta", "anc", "delivery")',
    },
  },
  {
    description:
      "Admin level used to join facilities to corresponding geo-consistency",
    replacementString: "GEOLEVEL",
    input: {
      inputType: "select",
      options: [
        {
          value: "admin_area_2",
          label: "admin_area_2",
        },
        {
          value: "admin_area_3",
          label: "admin_area_3",
        },
        {
          value: "admin_area_4",
          label: "admin_area_4",
        },
      ],
      valueType: "string",
      defaultValue: "admin_area_3",
    },
  },
];
