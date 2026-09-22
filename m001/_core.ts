import type { ModuleDefinitionCore } from "../.validation/_module_definition_github.ts";

export const core: ModuleDefinitionCore = {
  label: {
    en: "Data quality assessment",
    fr: "Évaluation de la qualité des données",
    pt: "Avaliação da qualidade dos dados",
  },
  family: "hmis",
  tier: "secondary",
  sortOrder: 1,
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
