import type { ResultsObjectDefinitionJSON } from "../.validation/module_definition_validator.ts";

export const resultsObjects: ResultsObjectDefinitionJSON[] = [
  {
    id: "M7_output_scorecard_admin_area_2.csv",
    description: "Scorecard indicators at zone level (long format)",
    createTableStatementPossibleColumns: {
      admin_area_2: "TEXT NOT NULL",
      quarter_id: "INTEGER NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      value: "NUMERIC",
    },
  },
  {
    id: "M7_output_scorecard_admin_area_3.csv",
    description: "Scorecard indicators at state level (long format)",
    createTableStatementPossibleColumns: {
      admin_area_2: "TEXT NOT NULL",
      admin_area_3: "TEXT NOT NULL",
      quarter_id: "INTEGER NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      value: "NUMERIC",
    },
  },
  {
    id: "M7_output_scorecard_admin_area_4.csv",
    description: "Scorecard indicators at LGA level (long format)",
    createTableStatementPossibleColumns: {
      admin_area_2: "TEXT NOT NULL",
      admin_area_3: "TEXT NOT NULL",
      admin_area_4: "TEXT NOT NULL",
      quarter_id: "INTEGER NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      value: "NUMERIC",
    },
  },
];
