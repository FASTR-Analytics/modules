import type { ResultsObjectDefinitionGithub } from "../.validation/_module_definition_github.ts";

export const resultsObjects: ResultsObjectDefinitionGithub[] = [
  {
    id: "M9_iceh_data.csv",
    createTableStatementPossibleColumns: {
      iceh_indicator: "TEXT NOT NULL",
      year: "INTEGER NOT NULL",
      source: "TEXT NOT NULL",
      strat: "TEXT NOT NULL",
      level: "TEXT NOT NULL",
      estimate: "NUMERIC",
      standard_error: "NUMERIC",
      sample_size: "NUMERIC",
    },
  },
  {
    id: "M9_iceh_inequality.csv",
    createTableStatementPossibleColumns: {
      iceh_indicator: "TEXT NOT NULL",
      year: "INTEGER NOT NULL",
      source: "TEXT NOT NULL",
      strat: "TEXT NOT NULL",
      ratio: "NUMERIC",
      difference: "NUMERIC",
      cix: "NUMERIC",
      sii: "NUMERIC",
    },
  },
];
