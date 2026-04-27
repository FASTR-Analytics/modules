import type { ResultsObjectDefinitionGithub } from "../.validation/_module_definition_github.ts";

export const resultsObjects: ResultsObjectDefinitionGithub[] = [
  {
    id: "M6_coverage_estimation_national.csv",
    createTableStatementPossibleColumns: {
      admin_area_1: "TEXT NOT NULL",
      year: "INTEGER NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      denominator: "TEXT NOT NULL",
      coverage_original_estimate: "NUMERIC",
      coverage_avgsurveyprojection: "NUMERIC",
      coverage_cov: "NUMERIC",
      survey_raw_source: "TEXT",
      survey_raw_source_detail: "TEXT",
    },
  },
  {
    id: "M6_coverage_estimation_admin2.csv",
    createTableStatementPossibleColumns: {
      admin_area_1: "TEXT NOT NULL",
      admin_area_2: "TEXT NOT NULL",
      year: "INTEGER NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      denominator: "TEXT NOT NULL",
      coverage_original_estimate: "NUMERIC",
      coverage_avgsurveyprojection: "NUMERIC",
      coverage_cov: "NUMERIC",
      survey_raw_source: "TEXT",
      survey_raw_source_detail: "TEXT",
    },
  },
  {
    id: "M6_coverage_estimation_admin3.csv",
    createTableStatementPossibleColumns: {
      admin_area_1: "TEXT NOT NULL",
      admin_area_3: "TEXT NOT NULL",
      year: "INTEGER NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      denominator: "TEXT NOT NULL",
      coverage_original_estimate: "NUMERIC",
      coverage_avgsurveyprojection: "NUMERIC",
      coverage_cov: "NUMERIC",
      survey_raw_source: "TEXT",
      survey_raw_source_detail: "TEXT",
    },
  },
];
