import type { ModuleDefinitionCore } from "../.validation/_module_definition_github.ts";

export const core: ModuleDefinitionCore = {
  label: {
    en: "M8. National Health Sector Scorecard (Catalog-Driven)",
    fr: "M8. Tableau de bord du secteur national de la santé (piloté par catalogue)",
  },
  prerequisites: ["m002"],
  scriptGenerationType: "calculated_indicators",
  assetsToImport: ["population.csv"],
  dataSources: [
    {
      sourceType: "results_object",
      replacementString: "M2_adjusted_data.csv",
      resultsObjectId: "M2_adjusted_data.csv",
      moduleId: "m002",
    },
  ],
};
