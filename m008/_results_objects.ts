import type { ResultsObjectDefinitionGithub } from "../.validation/_module_definition_github.ts";

export const resultsObjects: ResultsObjectDefinitionGithub[] = [
  {
    id: "M8_output_scorecard.csv",
    createTableStatementPossibleColumns: {
      admin_area_2: "TEXT NOT NULL",
      admin_area_3: "TEXT NOT NULL",
      admin_area_4: "TEXT NOT NULL",
      period_id: "INTEGER NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      numerator: "NUMERIC",
      denominator: "NUMERIC",
    },
  },
];
