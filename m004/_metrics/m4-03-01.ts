import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [];

///////////////////////////////////////////////////////////////
//  __       __              __                __            //
// /  \     /  |            /  |              /  |           //
// $$  \   /$$ |  ______   _$$ |_     ______  $$/   _______  //
// $$$  \ /$$$ | /      \ / $$   |   /      \ /  | /       | //
// $$$$  /$$$$ |/$$$$$$  |$$$$$$/   /$$$$$$  |$$ |/$$$$$$$/  //
// $$ $$ $$/$$ |$$    $$ |  $$ | __ $$ |  $$/ $$ |$$ |       //
// $$ |$$$/ $$ |$$$$$$$$/   $$ |/  |$$ |      $$ |$$ \_____  //
// $$ | $/  $$ |$$       |  $$  $$/ $$ |      $$ |$$       | //
// $$/      $$/  $$$$$$$/    $$$$/  $$/       $$/  $$$$$$$/  //
//                                                           //
///////////////////////////////////////////////////////////////

export const metric: MetricDefinitionGithub = {
  id: "m4-03-01",
  resultsObjectId: "M4_coverage_estimation_admin_area_3.csv",
  valueProps: ["coverage_cov"],
  valueFunc: "AVG",
  valueLabelReplacements: {
    coverage_cov: "Coverage calculated from HMIS data",
  },
  label: {
    en: "Coverage calculated from HMIS data",
    fr: "Couverture calculée à partir des données HMIS",
  },
  variantLabel: {
    en: "Admin Area 3",
    fr: "Zone administrative 3",
  },
  requiredDisaggregationOptions: [
    "indicator_common_id",
    "admin_area_3",
    "year",
  ],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "HMIS-derived health service coverage at admin area 3 (district) level.",
      fr: "Couverture des services de santé dérivée du HMIS au niveau de la zone administrative 3 (district).",
    },
    methodology: {
      en: "AVG of district-level coverage calculated as HMIS service volumes divided by district population denominators. Finest geographic resolution for coverage monitoring.",
      fr: "Moyenne de la couverture au niveau du district calculée comme volumes de services HMIS divisés par dénominateurs de population du district.",
    },
    interpretation: {
      en: "Enables district-level targeting and operational planning. Useful for identifying micro-level coverage gaps. Interpret with caution if sample sizes small.",
      fr: "Permet le ciblage au niveau du district et la planification opérationnelle. Utile pour identifier les lacunes de couverture au micro-niveau.",
    },
    typicalRange: {
      en: "0-100%. Greater variation expected than higher geographic levels due to smaller denominators and sample sizes.",
      fr: "0-100%. Plus grande variation attendue que les niveaux géographiques supérieurs en raison de plus petits dénominateurs.",
    },
    caveats: {
      en: "District-level denominators may have substantial uncertainty. Population mobility and small sample sizes increase volatility. Only available when ANALYSIS_LEVEL includes admin_area_3.",
      fr: "Les dénominateurs au niveau du district peuvent avoir une incertitude substantielle. Disponible uniquement lorsque ANALYSIS_LEVEL inclut admin_area_3.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id, admin_area_3, and year (all required). Consider aggregating to admin_area_2 if district-level estimates appear unstable.",
      fr: "Toujours désagréger par indicator_common_id, admin_area_3 et year (tous requis). Considérer l'agrégation à admin_area_2 si les estimations au niveau district semblent instables.",
    },
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
