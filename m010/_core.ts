import type { ModuleDefinitionCore } from "../.validation/_module_definition_github.ts";

export const core: ModuleDefinitionCore = {
  label: {
    en: "Health facility assessment",
    fr: "Évaluation des établissements de santé",
    pt: "Avaliação de unidades sanitárias",
  },
  family: "hfa",
  tier: "primary",
  sortOrder: 1,
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
