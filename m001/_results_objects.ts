import type { ResultsObjectDefinitionGithub } from "../.validation/_module_definition_github.ts";

export const resultsObjects: ResultsObjectDefinitionGithub[] = [
  {
    id: "M1_output_outliers.csv",
    description:
      "Detailed facility-level data with identified outliers and adjusted volume..",
    createTableStatementPossibleColumns: {
      facility_id: "TEXT NOT NULL",
      admin_area_4: "TEXT NOT NULL",
      admin_area_3: "TEXT NOT NULL",
      admin_area_2: "TEXT NOT NULL",
      period_id: "INTEGER NOT NULL",
      quarter_id: "INTEGER NOT NULL",
      year: "INTEGER NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      outlier_flag: "INTEGER NOT NULL",
    },
  },
  {
    id: "M1_output_completeness.csv",
    description:
      "Facility-level completeness data in a detailed long format, including reported and expected months",
    createTableStatementPossibleColumns: {
      facility_id: "TEXT NOT NULL",
      admin_area_4: "TEXT NOT NULL",
      admin_area_3: "TEXT NOT NULL",
      admin_area_2: "TEXT NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      period_id: "INTEGER NOT NULL",
      quarter_id: "INTEGER NOT NULL",
      year: "INTEGER NOT NULL",
      completeness_flag: "INTEGER NOT NULL",
    },
  },
  {
    id: "M1_output_consistency_geo.csv",
    description: "District-level consistency results",
    createTableStatementPossibleColumns: {
      admin_area_4: "TEXT NOT NULL",
      admin_area_3: "TEXT NOT NULL",
      admin_area_2: "TEXT NOT NULL",
      period_id: "INTEGER NOT NULL",
      quarter_id: "INTEGER NOT NULL",
      year: "INTEGER NOT NULL",
      ratio_type: "TEXT NOT NULL",
      sconsistency: "INTEGER",
    },
  },
  {
    id: "M1_output_dqa.csv",
    description: "Facility-level results from DQA analysis",
    createTableStatementPossibleColumns: {
      facility_id: "TEXT NOT NULL",
      admin_area_4: "TEXT NOT NULL",
      admin_area_3: "TEXT NOT NULL",
      admin_area_2: "TEXT NOT NULL",
      period_id: "INTEGER NOT NULL",
      quarter_id: "INTEGER NOT NULL",
      year: "INTEGER NOT NULL",
      dqa_mean: "NUMERIC NOT NULL",
      dqa_score: "NUMERIC NOT NULL",
    },
  },
  {
    id: "M1_output_outlier_list.csv",
    description: "Outlier list",
    createTableStatementPossibleColumns: {
      facility_id: "TEXT NOT NULL",
      admin_area_4: "TEXT NOT NULL",
      admin_area_3: "TEXT NOT NULL",
      admin_area_2: "TEXT NOT NULL",
      admin_area_1: "TEXT NOT NULL",
      period_id: "INTEGER NOT NULL",
      indicator_common_id: "TEXT NOT NULL",
      count: "NUMERIC NOT NULL",
    },
  },
];
