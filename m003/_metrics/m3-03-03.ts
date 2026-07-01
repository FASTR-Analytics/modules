import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [];

export const metric: MetricDefinitionGithub = {
  id: "m3-03-03",
  resultsObjectId: "M3_disruptions_analysis_admin_area_2.csv",
  label: {
    en: "Actual vs expected service volume",
    fr: "Volume de services réel vs attendu",
  },
  variantLabel: {
    en: "Admin area 2",
    fr: "Zone administrative 2",
  },
  valueProps: ["count_sum", "count_expect_sum"],
  valueFunc: "SUM",
  valueLabelReplacements: {
    count_sum: "Actual service volume",
    count_expect_sum: "Expected service volume",
  },
  requiredDisaggregationOptions: ["indicator_common_id", "admin_area_2"],
  formatAs: "number",
  aiDescription: {
    summary: {
      en: "Comparison of actual reported service volumes against model-predicted expected volumes at admin area 2 (province/state) level.",
      fr: "Comparaison des volumes de services réels contre volumes attendus au niveau de la zone administrative 2 (province/état).",
    },
    methodology: {
      en: "SUM of actual vs expected counts from area-specific robust regression models. Expected volumes account for local time trends and seasonal patterns. The full expected series is reported for every period (no difference threshold is applied).",
      fr: "Somme des comptes réels vs attendus des modèles de régression robustes spécifiques à la zone. Les volumes attendus tiennent compte des tendances temporelles locales et des modèles saisonniers. La série attendue complète est rapportée pour chaque période (aucun seuil de différence n'est appliqué).",
    },
    interpretation: {
      en: "Enables subnational identification of service disruptions. Compare across admin areas to identify geographic hotspots of service delivery problems. Areas with persistent negative deviations need targeted support.",
      fr: "Permet l'identification sous-nationale des perturbations de services. Comparer entre les zones administratives pour identifier les points chauds géographiques de problèmes de prestation de services. Les zones avec des écarts négatifs persistants nécessitent un soutien ciblé.",
    },
    typicalRange: {
      en: "Expected to closely match actual in stable regions. Deviations indicate local disruptions, data issues, or demand changes.",
      fr: "Attendu pour correspondre étroitement au réel dans les régions stables. Les écarts indiquent des perturbations locales, des problèmes de données ou des changements de la demande.",
    },
    caveats: {
      en: "Smaller geographic areas may have more volatile patterns, making model predictions less reliable. Consider aggregating to higher levels if area-specific models perform poorly.",
      fr: "Les zones géographiques plus petites peuvent avoir des modèles plus volatils, rendant les prédictions du modèle moins fiables. Considérer l'agrégation à des niveaux supérieurs si les modèles spécifiques aux zones fonctionnent mal.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id and admin_area_2 (both required). Time series reveals when and where disruptions occurred. Map visualization effectively shows geographic distribution of service gaps.",
      fr: "Toujours désagréger par indicator_common_id et admin_area_2 (tous deux requis). Les séries temporelles révèlent quand et où les perturbations se sont produites. La visualisation cartographique montre efficacement la distribution géographique des écarts de services.",
    },
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
