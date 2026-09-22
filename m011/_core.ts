import type { ModuleDefinitionCore } from "../.validation/_module_definition_github.ts";

export const core: ModuleDefinitionCore = {
  label: {
    en: "Disruption detection",
    fr: "Détection des perturbations",
    pt: "Detecção de perturbações",
  },
  family: "hmis",
  tier: "secondary",
  sortOrder: 3,
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
