import type { ResultsObjectDefinitionGithub } from "../.validation/_module_definition_github.ts";

export const resultsObjects: ResultsObjectDefinitionGithub[] = [
  {
    id: "M7_output_scorecard_admin_area_2.csv",
    createTableStatementPossibleColumns: {
      admin_area_2: "TEXT NOT NULL",
      quarter_id: "INTEGER NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      value: "NUMERIC",
    },
  },
  {
    id: "M7_output_scorecard_admin_area_3.csv",
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
