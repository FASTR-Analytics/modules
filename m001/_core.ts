import type { ModuleDefinitionCore } from "../.validation/_module_definition_github.ts";

export const core: ModuleDefinitionCore = {
  label: {
    en: "M1. Data quality assessment",
    fr: "M1. Évaluation de la qualité des données",
    pt: "M1. Avaliação da qualidade dos dados",
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
