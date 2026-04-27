import type { ResultsObjectDefinitionGithub } from "../.validation/_module_definition_github.ts";

export const resultsObjects: ResultsObjectDefinitionGithub[] = [
  {
    id: "M2_adjusted_data.csv",
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
    id: "M2_adjusted_data_admin_area.csv",
    createTableStatementPossibleColumns: false,
  },
  {
    id: "M2_adjusted_data_national.csv",
    createTableStatementPossibleColumns: false,
  },
  {
    id: "M2_low_volume_exclusions.csv",
    createTableStatementPossibleColumns: {
      indicator_common_id: "TEXT NOT NULL",
      low_volume_exclude: "TEXT NOT NULL",
    },
  },
];
