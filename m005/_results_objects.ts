import type { ResultsObjectDefinitionJSON } from "../.validation/_module_definition_github.ts";

export const resultsObjects: ResultsObjectDefinitionJSON[] = [
  {
    id: "M5_denominators_national.csv",
    description: "Selected denominators (National)",
    createTableStatementPossibleColumns: {
      admin_area_1: "TEXT NOT NULL",
      year: "INTEGER NOT NULL",
      denominator: "TEXT NOT NULL",
      value: "NUMERIC NOT NULL",
      source_indicator: "TEXT",
      target_population: "TEXT",
    },
  },
  {
    id: "M5_denominators_admin2.csv",
    description: "Selected denominators (Admin area 2)",
    createTableStatementPossibleColumns: {
      admin_area_1: "TEXT NOT NULL",
      admin_area_2: "TEXT NOT NULL",
      year: "INTEGER NOT NULL",
      denominator: "TEXT NOT NULL",
      value: "NUMERIC NOT NULL",
      source_indicator: "TEXT",
      target_population: "TEXT",
    },
  },
  {
    id: "M5_denominators_admin3.csv",
    description: "Selected denominators (Admin area 3)",
    createTableStatementPossibleColumns: {
      admin_area_1: "TEXT NOT NULL",
      admin_area_3: "TEXT NOT NULL",
      year: "INTEGER NOT NULL",
      denominator: "TEXT NOT NULL",
      value: "NUMERIC NOT NULL",
      source_indicator: "TEXT",
      target_population: "TEXT",
    },
  },
  {
    id: "M5_combined_results_national.csv",
    description: "Selected denominators (National)",
    createTableStatementPossibleColumns: {
      admin_area_1: "TEXT NOT NULL",
      year: "INTEGER NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      denominator_best_or_survey: "TEXT NOT NULL",
      value: "NUMERIC NOT NULL",
    },
  },
  {
    id: "M5_combined_results_admin2.csv",
    description: "Selected denominators (National)",
    createTableStatementPossibleColumns: {
      admin_area_1: "TEXT NOT NULL",
      admin_area_2: "TEXT NOT NULL",
      year: "INTEGER NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      denominator_best_or_survey: "TEXT NOT NULL",
      value: "NUMERIC NOT NULL",
    },
  },
  {
    id: "M5_combined_results_admin3.csv",
    description: "Selected denominators (National)",
    createTableStatementPossibleColumns: {
      admin_area_1: "TEXT NOT NULL",
      admin_area_3: "TEXT NOT NULL",
      year: "INTEGER NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      denominator_best_or_survey: "TEXT NOT NULL",
      value: "NUMERIC NOT NULL",
    },
  },
  {
    id: "M5_selected_denominator_per_indicator.csv",
    description: "Selected denominators",
    createTableStatementPossibleColumns: {
      indicator_common_id: "TEXT NOT NULL",
      denominator_national: "TEXT NOT NULL",
      denominator_admin2: "TEXT NOT NULL",
      denominator_admin3: "TEXT NOT NULL",
    },
  },
];
