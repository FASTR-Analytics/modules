import type { ModuleDefinitionCore } from "../.validation/_module_definition_github.ts";

export const core: ModuleDefinitionCore = {
  label: {
    en: "M5. Coverage estimates ~ new, part 1",
    fr: "M5. Estimations de couverture ~ nouveau, partie 1",
    pt: "M5. Estimativas de cobertura ~ novo, parte 1",
  },
  prerequisites: ["m002"],
  scriptGenerationType: "template",
  dataSources: [
    {
      replacementString: "M2_adjusted_data_national.csv",
      sourceType: "results_object",
      resultsObjectId: "M2_adjusted_data_national.csv",
      moduleId: "m002",
    },
    {
      replacementString: "M2_adjusted_data_admin_area.csv",
      sourceType: "results_object",
      resultsObjectId: "M2_adjusted_data_admin_area.csv",
      moduleId: "m002",
    },
  ],
  assetsToImport: [
    {
      name: "survey_data_unified.csv",
      repoPath: "survey_data_unified.csv",
      commit: "baa11785b67cbf8eed858e61c1a995ca9bb029bd",
    },
    {
      name: "population_estimates_only.csv",
      repoPath: "population_estimates_only.csv",
      commit: "baa11785b67cbf8eed858e61c1a995ca9bb029bd",
    },
  ],
};
