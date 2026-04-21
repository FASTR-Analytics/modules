import type { ResultsObjectDefinitionGithub } from "../.validation/_module_definition_github.ts";

export const resultsObjects: ResultsObjectDefinitionGithub[] = [
  {
    id: "M4_coverage_estimation.csv",
    description: "Coverage estimates (National)",
    createTableStatementPossibleColumns: {
      indicator_common_id: "TEXT NOT NULL",
      year: "INTEGER NOT NULL",
      coverage_original_estimate: "NUMERIC",
      coverage_avgsurveyprojection: "NUMERIC",
      coverage_cov: "NUMERIC",
    },
  },
  {
    id: "M4_coverage_estimation_admin_area_2.csv",
    description: "Coverage results (sub-national level)",
    createTableStatementPossibleColumns: {
      admin_area_2: "TEXT NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      year: "INTEGER NOT NULL",
      coverage_cov: "NUMERIC",
    },
  },
  {
    id: "M4_coverage_estimation_admin_area_3.csv",
    description: "Coverage results (sub-national level)",
    createTableStatementPossibleColumns: {
      admin_area_3: "TEXT NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      year: "INTEGER NOT NULL",
      coverage_cov: "NUMERIC",
    },
  },
  {
    id: "M4_selected_denominator_per_indicator.csv",
    description: "Selected denominators",
    createTableStatementPossibleColumns: {
      indicator_common_id: "TEXT NOT NULL",
      denominator_national: "TEXT NOT NULL",
      denominator_admin2: "TEXT NOT NULL",
      denominator_admin3: "TEXT NOT NULL",
    },
  },
];
