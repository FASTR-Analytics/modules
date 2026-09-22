import type { ModuleDefinitionCore } from "../.validation/_module_definition_github.ts";

export const core: ModuleDefinitionCore = {
  label: {
    en: "Coverage estimates",
    fr: "Estimations de couverture",
    pt: "Estimativas de cobertura",
  },
  family: "hmis",
  tier: "secondary",
  sortOrder: 5,
  prerequisites: ["m005"],
  scriptGenerationType: "template",
  dataSources: [
    {
      replacementString: "M5_combined_results_national.csv",
      sourceType: "results_object",
      resultsObjectId: "M5_combined_results_national.csv",
      moduleId: "m005",
    },
    {
      replacementString: "M5_combined_results_admin2.csv",
      sourceType: "results_object",
      resultsObjectId: "M5_combined_results_admin2.csv",
      moduleId: "m005",
    },
    {
      replacementString: "M5_combined_results_admin3.csv",
      sourceType: "results_object",
      resultsObjectId: "M5_combined_results_admin3.csv",
      moduleId: "m005",
    },
  ],
  assetsToImport: [],
};
