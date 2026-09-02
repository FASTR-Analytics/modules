import type { ModuleDefinitionCore } from "../.validation/_module_definition_github.ts";

// M12 materialises the additive INGREDIENTS of every common indicator, so the
// app can apply each indicator's own formula after aggregation and get an
// exact answer at any grouping (PLAN_1a §0).
//
// One upstream module, one upstream file, and the run's person-years file.
// The ingredient table is NOT a data source: the app substitutes it into
// script.R as a tribble literal in place of the INDICATOR_INGREDIENTS token
// (PLAN_1a §1.5), the same channel as COUNTRY_ISO3 and every module
// parameter. It therefore rides in the script text the memoization key
// already hashes, so an indicator edit re-runs this module by construction.
//
// The population source (PLAN_1b) is the app-written inputs/population.csv:
// monthly person-years per population type at the finest admin level,
// expanded from the instance population store at capture. It IS a data
// source — its content hash enters this module's memoization key, so a
// population edit re-runs the module and an unchanged store does not.
export const core: ModuleDefinitionCore = {
  label: {
    en: "M12. Indicator values",
    fr: "M12. Valeurs des indicateurs",
    pt: "M12. Valores dos indicadores",
  },
  prerequisites: ["m002"],
  scriptGenerationType: "template",
  assetsToImport: [],
  dataSources: [
    {
      sourceType: "results_object",
      replacementString: "M2_adjusted_data.csv",
      resultsObjectId: "M2_adjusted_data.csv",
      moduleId: "m002",
    },
    {
      sourceType: "population",
      replacementString: "POPULATION_PERSON_YEARS",
    },
  ],
};
