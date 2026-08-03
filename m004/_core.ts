import type { ModuleDefinitionCore } from "../.validation/_module_definition_github.ts";

export const core: ModuleDefinitionCore = {
  label: {
    en: "M4. Coverage estimates",
    fr: "M4. Estimations de couverture",
    pt: "M4. Estimativas de cobertura",
  },
  prerequisites: ["m002"],
  scriptGenerationType: "template",
  dataSources: [
    {
      replacementString: "M2_adjusted_data_national.csv",
      sourceType: "results_object",
      resultsObjectId: "M2_adjusted_data_national.csv",
      moduleId: "m002",
    },
    {
      replacementString: "M2_adjusted_data_admin_area.csv",
      sourceType: "results_object",
      resultsObjectId: "M2_adjusted_data_admin_area.csv",
      moduleId: "m002",
    },
  ],
  assetsToImport: [
    {
      name: "survey_data_unified.csv",
      repoPath: "survey_data_unified.csv",
    },
    {
      name: "population_estimates_only.csv",
      repoPath: "population_estimates_only.csv",
    },
  ],
};
