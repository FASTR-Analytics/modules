import type { ResultsObjectDefinitionJSON } from "../.validation/module_definition_validator.ts";

export const resultsObjects: ResultsObjectDefinitionJSON[] = [
  {
    "id": "M6_coverage_estimation_national.csv",
    "description": "Coverage estimates (National)",
    "createTableStatementPossibleColumns": {
      "admin_area_1": "TEXT NOT NULL",
      "year": "INTEGER NOT NULL",
      "indicator_common_id": "TEXT NOT NULL",
      "denominator": "TEXT NOT NULL",
      "coverage_original_estimate": "NUMERIC",
      "coverage_avgsurveyprojection": "NUMERIC",
      "coverage_cov": "NUMERIC",
      "survey_raw_source": "TEXT",
      "survey_raw_source_detail": "TEXT"
    }
  },
  {
    "id": "M6_coverage_estimation_admin2.csv",
    "description": "Coverage results (Admin area 2)",
    "createTableStatementPossibleColumns": {
      "admin_area_1": "TEXT NOT NULL",
      "admin_area_2": "TEXT NOT NULL",
      "year": "INTEGER NOT NULL",
      "indicator_common_id": "TEXT NOT NULL",
      "denominator": "TEXT NOT NULL",
      "coverage_original_estimate": "NUMERIC",
      "coverage_avgsurveyprojection": "NUMERIC",
      "coverage_cov": "NUMERIC",
      "survey_raw_source": "TEXT",
      "survey_raw_source_detail": "TEXT"
    }
  },
  {
    "id": "M6_coverage_estimation_admin3.csv",
    "description": "Coverage results (Admin area 3)",
    "createTableStatementPossibleColumns": {
      "admin_area_1": "TEXT NOT NULL",
      "admin_area_3": "TEXT NOT NULL",
      "year": "INTEGER NOT NULL",
      "indicator_common_id": "TEXT NOT NULL",
      "denominator": "TEXT NOT NULL",
      "coverage_original_estimate": "NUMERIC",
      "coverage_avgsurveyprojection": "NUMERIC",
      "coverage_cov": "NUMERIC",
      "survey_raw_source": "TEXT",
      "survey_raw_source_detail": "TEXT"
    }
  }
];
