import type { MetricDefinitionJSON } from "../.validation/module_definition_validator.ts";

const SCORECARD_LABEL_REPLACEMENTS: Record<string, string> = {
  anc4_anc1_before20_ratio: "ANC4 / ANC1 <20wks",
  anc4_anc1_ratio: "ANC4 / ANC1",
  skilled_birth_attendance: "Skilled Birth Attendant / Reported Deliveries",
  new_fp_acceptors_rate: "New FP Acceptors / Women of Reproductive Age",
  act_malaria_treatment: "ACT for Uncomplicated Malaria",
  penta3_coverage: "Penta 3",
  fully_immunized_coverage: "Fully Immunized",
  htn_new_per_10000: "HTN New per 10,000 person-years",
  diabetes_new_per_10000: " Diabetes New per 10,000 person-years",
  nhmis_data_timeliness_final: "NHMIS reports on time with content",
};

const PLACEHOLDER_AI_DESCRIPTION = {
  summary: { en: "xxx", fr: "xxx" },
  methodology: { en: "xxx", fr: "xxx" },
  interpretation: { en: "xxx", fr: "xxx" },
  typicalRange: { en: "xxx", fr: "xxx" },
  useCases: [
    { en: "xxx", fr: "xxx" },
    { en: "xxx", fr: "xxx" },
  ],
  disaggregationGuidance: { en: "xxx", fr: "xxx" },
};

export const metrics: MetricDefinitionJSON[] = [
  {
    id: "m7-01-01",
    resultsObjectId: "M7_output_scorecard_admin_area_2.csv",
    label: { en: "Scorecard (AA2)", fr: "Scorecard (AA2)" },
    valueProps: ["value"],
    valueFunc: "identity",
    valueLabelReplacements: SCORECARD_LABEL_REPLACEMENTS,
    requiredDisaggregationOptions: [
      "admin_area_2",
      "quarter_id",
      "indicator_common_id",
    ],
    formatAs: "percent",
    periodOptions: ["quarter_id"],
    aiDescription: PLACEHOLDER_AI_DESCRIPTION,
  },
  {
    id: "m7-01-02",
    resultsObjectId: "M7_output_scorecard_admin_area_3.csv",
    label: { en: "Scorecard (AA3)", fr: "Scorecard (AA3)" },
    valueProps: ["value"],
    valueFunc: "identity",
    valueLabelReplacements: SCORECARD_LABEL_REPLACEMENTS,
    requiredDisaggregationOptions: [
      "admin_area_3",
      "quarter_id",
      "indicator_common_id",
    ],
    formatAs: "percent",
    periodOptions: ["quarter_id"],
    aiDescription: PLACEHOLDER_AI_DESCRIPTION,
  },
  {
    id: "m7-01-03",
    resultsObjectId: "M7_output_scorecard_admin_area_4.csv",
    label: { en: "Scorecard (AA4)", fr: "Scorecard (AA4)" },
    valueProps: ["value"],
    valueFunc: "identity",
    valueLabelReplacements: SCORECARD_LABEL_REPLACEMENTS,
    requiredDisaggregationOptions: [
      "admin_area_4",
      "quarter_id",
      "indicator_common_id",
    ],
    formatAs: "percent",
    periodOptions: ["quarter_id"],
    aiDescription: PLACEHOLDER_AI_DESCRIPTION,
  },
];
