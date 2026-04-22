import type { ResultsObjectDefinitionGithub } from "../.validation/_module_definition_github.ts";

export const resultsObjects: ResultsObjectDefinitionGithub[] = [
  {
    id: "M3_service_utilization.csv",
    description:
      "Service utilization data with adjusted volumes for all adjustment scenarios",
    createTableStatementPossibleColumns: {
      facility_id: "TEXT NOT NULL",
      admin_area_4: "TEXT NOT NULL",
      admin_area_3: "TEXT NOT NULL",
      admin_area_2: "TEXT NOT NULL",
      period_id: "INTEGER NOT NULL",
      quarter_id: "INTEGER NOT NULL",
      year: "INTEGER NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      count_final_none: "NUMERIC",
      count_final_outliers: "NUMERIC",
      count_final_completeness: "NUMERIC",
      count_final_both: "NUMERIC",
    },
  },
  {
    id: "M3_disruptions_analysis_admin_area_1.csv",
    description: "National-level disruption analysis results",
    createTableStatementPossibleColumns: {
      admin_area_1: "TEXT NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      period_id: "INTEGER NOT NULL",
      quarter_id: "INTEGER NOT NULL",
      year: "INTEGER NOT NULL",
      count_sum: "NUMERIC",
      count_expect_sum: "NUMERIC",
      count_expected_if_above_diff_threshold: "NUMERIC",
    },
  },
  {
    id: "M3_disruptions_analysis_admin_area_2.csv",
    description: "Admin area 2 level disruption analysis results",
    createTableStatementPossibleColumns: {
      admin_area_2: "TEXT NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      period_id: "INTEGER NOT NULL",
      quarter_id: "INTEGER NOT NULL",
      year: "INTEGER NOT NULL",
      count_sum: "NUMERIC",
      count_expect_sum: "NUMERIC",
      count_expected_if_above_diff_threshold: "NUMERIC",
    },
  },
  {
    id: "M3_disruptions_analysis_admin_area_3.csv",
    description: "Admin area 3 level disruption analysis results",
    createTableStatementPossibleColumns: {
      admin_area_2: "TEXT NOT NULL",
      admin_area_3: "TEXT NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      period_id: "INTEGER NOT NULL",
      quarter_id: "INTEGER NOT NULL",
      year: "INTEGER NOT NULL",
      count_sum: "NUMERIC",
      count_expect_sum: "NUMERIC",
      count_expected_if_above_diff_threshold: "NUMERIC",
    },
  },
  {
    id: "M3_disruptions_analysis_admin_area_4.csv",
    description: "Admin area 4 level disruption analysis results",
    createTableStatementPossibleColumns: {
      admin_area_2: "TEXT NOT NULL",
      admin_area_3: "TEXT NOT NULL",
      admin_area_4: "TEXT NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      period_id: "INTEGER NOT NULL",
      quarter_id: "INTEGER NOT NULL",
      year: "INTEGER NOT NULL",
      count_sum: "NUMERIC",
      count_expect_sum: "NUMERIC",
      count_expected_if_above_diff_threshold: "NUMERIC",
    },
  },
  {
    id: "M3_chartout.csv",
    description: "Control chart analysis results with tagged anomalies",
    createTableStatementPossibleColumns: false,
  },
  {
    id: "M3_all_indicators_shortfalls_admin_area_1.csv",
    description:
      "Shortfall and surplus calculations for all indicators (Admin area 1)",
    createTableStatementPossibleColumns: false,
  },
  {
    id: "M3_all_indicators_shortfalls_admin_area_2.csv",
    description:
      "Shortfall and surplus calculations for all indicators (Admin area 2)",
    createTableStatementPossibleColumns: false,
  },
  {
    id: "M3_all_indicators_shortfalls_admin_area_3.csv",
    description:
      "Shortfall and surplus calculations for all indicators (Admin area 3)",
    createTableStatementPossibleColumns: false,
  },
  {
    id: "M3_all_indicators_shortfalls_admin_area_4.csv",
    description:
      "Shortfall and surplus calculations for all indicators (Admin area 4)",
    createTableStatementPossibleColumns: false,
  },
];
