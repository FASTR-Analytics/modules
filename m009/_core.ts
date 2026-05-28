import type { ModuleDefinitionCore } from "../.validation/_module_definition_github.ts";

export const core: ModuleDefinitionCore = {
  label: {
    en: "M9. ICEH Survey Data Analysis",
    fr: "M9. Analyse des données d'enquête ICEH",
  },
  prerequisites: [],
  scriptGenerationType: "template",
  dataSources: [
    {
      sourceType: "dataset",
      replacementString: "PROJECT_DATA_ICEH",
      datasetType: "iceh",
    },
  ],
  assetsToImport: [],
};
