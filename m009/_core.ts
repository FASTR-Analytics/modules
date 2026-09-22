import type { ModuleDefinitionCore } from "../.validation/_module_definition_github.ts";

export const core: ModuleDefinitionCore = {
  label: {
    en: "ICEH survey analysis",
    fr: "Analyse de l'enquête ICEH",
    pt: "Análise do inquérito ICEH",
  },
  family: "iceh",
  tier: "primary",
  sortOrder: 1,
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
