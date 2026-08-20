import type { ModuleDefinitionCore } from "../.validation/_module_definition_github.ts";

export const core: ModuleDefinitionCore = {
  label: {
    en: "M11. Bayesian disruption detection (LI model)",
    fr: "M11. Détection bayésienne des perturbations (modèle LI)",
  },
  prerequisites: ["m002"],
  scriptGenerationType: "template",
  dataSources: [
    {
      sourceType: "dataset",
      replacementString: "PROJECT_DATA_HMIS",
      datasetType: "hmis",
    },
    {
      replacementString: "M2_adjusted_data.csv",
      sourceType: "results_object",
      resultsObjectId: "M2_adjusted_data.csv",
      moduleId: "m002",
    },
  ],
  assetsToImport: [],
};
