import type { ModuleDefinitionCore } from "../.validation/_module_definition_github.ts";

export const core: ModuleDefinitionCore = {
  label: {
    en: "M1. Data quality assessment x",
    fr: "M1. Évaluation de la qualité des données",
  },
  prerequisites: [],
  scriptGenerationType: "template",
  dataSources: [
    {
      sourceType: "dataset",
      replacementString: "PROJECT_DATA_HMIS",
      datasetType: "hmis",
    },
  ],
  assetsToImport: [],
};
