import type { ModuleDefinitionCore } from "../.validation/_module_definition_github.ts";

// M12 materialises the additive INGREDIENTS of every common indicator, so the
// app can apply each indicator's own formula after aggregation and get an
// exact answer at any grouping (PLAN_1a §0).
//
// One upstream module and one upstream file. The ingredient table is NOT a
// data source: the app substitutes it into script.R as a tribble literal in
// place of the INDICATOR_INGREDIENTS token (PLAN_1a §1.5), the same channel as
// COUNTRY_ISO3 and every module parameter. It therefore rides in the script
// text the memoization key already hashes, so an indicator edit re-runs this
// module by construction.
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
  ],
};
