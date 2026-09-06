import type { ResultsObjectDefinitionGithub } from "../.validation/_module_definition_github.ts";

// One row per indicator × month × admin area at the instance's population
// level (the level of the stored population data, else the HMIS depth),
// carrying that indicator's additive ingredients in numbered slots. AREA ×
// MONTH grain uniformly: rows are summed across facilities and any finer
// admin level, so there is no facility_id and no facility column (facility
// analysis stays on m3-01-01 and the quality metrics). The admin columns
// present in a package are those of its person-years file.
//
// The slots are the generalisation of m008's numerator/denominator pair from
// two columns to eight, which is what lets an indicator be an arbitrary
// expression rather than a ratio.
export const resultsObjects: ResultsObjectDefinitionGithub[] = [
  {
    id: "M12_indicator_values.csv",
    createTableStatementPossibleColumns: {
      indicator_common_id: "TEXT NOT NULL",
      period_id: "INTEGER NOT NULL",
      admin_area_2: "TEXT NOT NULL",
      admin_area_3: "TEXT NOT NULL",
      admin_area_4: "TEXT NOT NULL",
      ing1: "NUMERIC",
      ing2: "NUMERIC",
      ing3: "NUMERIC",
      ing4: "NUMERIC",
      ing5: "NUMERIC",
      ing6: "NUMERIC",
      ing7: "NUMERIC",
      ing8: "NUMERIC",
    },
  },
];
