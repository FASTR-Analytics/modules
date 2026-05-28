import type { ResultsObjectDefinitionGithub } from "../.validation/_module_definition_github.ts";

export const resultsObjects: ResultsObjectDefinitionGithub[] = [
  {
    id: "M10_hfa_results.csv",
    createTableStatementPossibleColumns: {
      facility_id: "TEXT NOT NULL",
      admin_area_4: "TEXT NOT NULL",
      admin_area_3: "TEXT NOT NULL",
      admin_area_2: "TEXT NOT NULL",
      admin_area_1: "TEXT NOT NULL",
      hfa_indicator: "TEXT NOT NULL",
      hfa_category: "TEXT",
      hfa_sub_category: "TEXT",
      time_point: "TEXT NOT NULL",
      facility_ownership: "TEXT NOT NULL",
      facility_type: "TEXT NOT NULL",
      facility_custom_1: "TEXT NOT NULL",
      facility_custom_2: "TEXT NOT NULL",
      facility_custom_3: "TEXT NOT NULL",
      facility_custom_4: "TEXT NOT NULL",
      facility_custom_5: "TEXT NOT NULL",
      numeric_sum: "NUMERIC",
      numeric_avg: "NUMERIC",
      boolean_sum: "NUMERIC",
      boolean_avg: "NUMERIC",
    },
  },
];
