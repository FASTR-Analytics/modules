import type { ModuleDefinitionCore } from "../.validation/module_definition_validator.ts";

export const core: ModuleDefinitionCore = {
  label: {
    en: "M3. Service utilization",
    fr: "M3. Utilisation des services",
  },
  prerequisites: ["m001", "m002"],
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
      replacementString: "M2_adjusted_data_admin_area.csv",
      sourceType: "results_object",
      resultsObjectId: "M2_adjusted_data_admin_area.csv",
      moduleId: "m002",
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
