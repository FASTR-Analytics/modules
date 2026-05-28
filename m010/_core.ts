import type { ModuleDefinitionCore } from "../.validation/_module_definition_github.ts";

export const core: ModuleDefinitionCore = {
  label: {
    en: "M10. Health facility assessment",
    fr: "M10. Évaluation des établissements de santé",
  },
  prerequisites: [],
  scriptGenerationType: "hfa",
  dataSources: [
    {
      sourceType: "dataset",
      replacementString: "PROJECT_DATA_HFA",
      datasetType: "hfa",
    },
  ],
  assetsToImport: [],
};
