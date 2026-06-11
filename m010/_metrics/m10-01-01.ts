import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { CF_80_70 } from "../../.validation/cf_presets.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "hfa-category-indicator-table",
    label: {
      en: "Indicators by time point (per service category)",
      fr: "Indicateurs par point temporel (par catégorie de service)",
    },
    description: {
      en: "Table of HFA indicator values for a selected service category, with indicators as rows and time points as columns. Switch the category using the replicant selector.",
      fr: "Tableau des valeurs des indicateurs HFA pour une catégorie de service sélectionnée, avec les indicateurs en lignes et les points temporels en colonnes. Changez de catégorie à l'aide du sélecteur de réplicant.",
    },
    allowedFilters: [
      "hfa_category",
      "hfa_sub_category",
      "hfa_indicator",
      "time_point",
    ],
    createDefaultVisualizationOnInstall: "82617c89-803b-4139-b000-ebddf2469dfd",
    config: {
      d: {
        type: "table",
        valuesDisDisplayOpt: "col",
        disaggregateBy: [
          { disOpt: "hfa_category", disDisplayOpt: "replicant" },
          { disOpt: "hfa_indicator", disDisplayOpt: "row" },
          { disOpt: "time_point", disDisplayOpt: "col" },
        ],
        selectedReplicantValue: "financing",
        filterBy: [],
      },
      s: {
        ...CF_80_70,
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "HFA indicators",
          fr: "Indicateurs HFA",
        },
        captionRelFontSize: null,
        subCaption: null,
        subCaptionRelFontSize: null,
        footnote: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
];

export const metric: MetricDefinitionGithub = {
  id: "m10-01-01",
  hide: false,
  resultsObjectId: "M10_hfa_results.csv",
  valueProps: ["value"],
  valueFunc: "identity",
  valueLabelReplacements: {},
  label: {
    en: "HFA indicators",
    fr: "Indicateurs HFA",
  },
  requiredDisaggregationOptions: ["hfa_indicator", "time_point"],
  formatAs: "number",
  postAggregationExpression: {
    ingredientValues: [
      { prop: "sum_val", func: "SUM" },
      { prop: "avg_num", func: "SUM" },
      { prop: "avg_weight", func: "SUM" },
    ],
    // Bare division only: the PO query evaluator auto-wraps "/column" in a
    // NULLIF guard; writing NULLIF here would be double-wrapped into bad SQL
    expression: "value = COALESCE(sum_val, avg_num / avg_weight)",
  },
  aiDescription: {
    summary: {
      en: "Health Facility Assessment indicator values aggregated according to each indicator's defined aggregation method (sum or average) and type (binary or numeric).",
      fr: "Valeurs des indicateurs d'évaluation des établissements de santé agrégées selon la méthode d'agrégation définie pour chaque indicateur (somme ou moyenne) et le type (binaire ou numérique).",
    },
    methodology: {
      en: "Each indicator is automatically aggregated using its configured method: binary indicators with avg aggregation yield percentages (0-1), binary with sum yield facility counts, numeric with avg yield mean values, numeric with sum yield totals.",
      fr: "Chaque indicateur est automatiquement agrégé selon sa méthode configurée: les indicateurs binaires avec agrégation moyenne donnent des pourcentages (0-1), binaires avec somme donnent des comptes d'établissements, numériques avec moyenne donnent des valeurs moyennes, numériques avec somme donnent des totaux.",
    },
    interpretation: {
      en: "Interpretation depends on the indicator type and aggregation method. Check the indicator definition to understand whether values represent percentages, counts, averages, or totals.",
      fr: "L'interprétation dépend du type d'indicateur et de la méthode d'agrégation. Vérifiez la définition de l'indicateur pour comprendre si les valeurs représentent des pourcentages, des comptes, des moyennes ou des totaux.",
    },
    typicalRange: {
      en: "Varies by indicator type: binary avg indicators range 0-1 (percentage), others vary based on what is being measured.",
      fr: "Varie selon le type d'indicateur: les indicateurs binaires moyens varient de 0 à 1 (pourcentage), les autres varient selon ce qui est mesuré.",
    },
    caveats: {
      en: "HFA data is typically cross-sectional and may not reflect temporal trends. Survey timing and facility sampling affect comparability across assessments. When sampling weights are enabled, sums are weighted population-total estimates (non-integer) and may be unreliable at fine admin-area disaggregation if the weights were designed for national or stratum-level estimates.",
      fr: "Les données HFA sont généralement transversales et peuvent ne pas refléter les tendances temporelles. Lorsque les pondérations d'échantillonnage sont activées, les sommes sont des estimations pondérées de totaux de population (non entières) et peuvent être peu fiables à une désagrégation administrative fine si les pondérations ont été conçues pour des estimations nationales ou par strate.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by hfa_indicator and time_point (both required). Use admin_area to compare regional results. Disaggregate by facility_type or facility_ownership to identify disparities.",
      fr: "Toujours désagréger par hfa_indicator et time_point (tous deux requis). Utiliser admin_area pour comparer les résultats régionaux.",
    },
  },
  variantLabel: null,
  importantNotes: null,
  vizPresets,
};
