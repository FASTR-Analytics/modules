import type { ModuleDefinitionCore } from "../.validation/module_definition_validator.ts";

export const core: ModuleDefinitionCore = {
  "label": {
    "en": "M5. Coverage estimates ~ new, part 1",
    "fr": "M5. Estimations de couverture ~ nouveau, partie 1"
  },
  "prerequisites": [
    "m002"
  ],
  "scriptGenerationType": "template",
  "dataSources": [
    {
      "replacementString": "M2_adjusted_data_national.csv",
      "sourceType": "results_object",
      "resultsObjectId": "M2_adjusted_data_national.csv",
      "moduleId": "m002"
    },
    {
      "replacementString": "M2_adjusted_data_admin_area.csv",
      "sourceType": "results_object",
      "resultsObjectId": "M2_adjusted_data_admin_area.csv",
      "moduleId": "m002"
    }
  ],
  "assetsToImport": [
    "survey_data_unified.csv",
    "population_estimates_only.csv",
    "ng_province_denominators_corrected.csv",
    "ng_national_denominators_corrected.csv",
    "chmis_national_for_module4.csv",
    "chmis_admin_area_for_module4.csv"
  ]
};
