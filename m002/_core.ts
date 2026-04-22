import type { ModuleDefinitionCore } from "../.validation/_module_definition_github.ts";

export const core: ModuleDefinitionCore = {
  label: {
    en: "M2. Data quality adjustments x",
    fr: "M2. Ajustements de la qualité des données",
  },
  prerequisites: ["m001"],
  scriptGenerationType: "template",
  dataSources: [
    {
      sourceType: "dataset",
      replacementString: "PROJECT_DATA_HMIS",
      datasetType: "hmis",
    },
    {
      replacementString: "M1_output_outliers.csv",
      sourceType: "results_object",
      resultsObjectId: "M1_output_outliers.csv",
      moduleId: "m001",
    },
    {
      replacementString: "M1_output_completeness.csv",
      sourceType: "results_object",
      resultsObjectId: "M1_output_completeness.csv",
      moduleId: "m001",
    },
  ],
  assetsToImport: [],
};
