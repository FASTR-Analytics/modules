import type { ModuleDefinitionCore } from "../.validation/module_definition_validator.ts";

export const core: ModuleDefinitionCore = {
  label: {
    en: "M7. National Health Sector Scorecard (NHSS)",
    fr: "M7. Tableau de bord du secteur national de la sante (NHSS)",
  },
  prerequisites: ["m002"],
  scriptGenerationType: "template",
  assetsToImport: ["total_population_NGA.csv"],
  dataSources: [
    {
      sourceType: "results_object",
      replacementString: "M2_adjusted_data.csv",
      resultsObjectId: "M2_adjusted_data.csv",
      moduleId: "m002",
    },
  ],
};
