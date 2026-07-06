/**
 * METRICS TEMPLATE
 *
 * Copy this file to your module directory as `_metrics.ts` and fill in the values.
 * See DOC_METRICS_REFERENCE.md for detailed explanations of each field.
 *
 * Delete this comment block when done.
 */

import type { MetricDefinitionGithub } from "../.validation/_module_definition_github.ts";
// import { CF_01_03, CF_80_70, CF_90_80 } from "../.validation/cf_presets.ts";

export const metrics: MetricDefinitionGithub[] = [
  //
  //
  //
  //
  //
  //
  //
  //
  //
  // ███╗   ███╗███████╗████████╗██████╗ ██╗ ██████╗
  // ████╗ ████║██╔════╝╚══██╔══╝██╔══██╗██║██╔════╝
  // ██╔████╔██║█████╗     ██║   ██████╔╝██║██║
  // ██║╚██╔╝██║██╔══╝     ██║   ██╔══██╗██║██║
  // ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║██║╚██████╗
  // ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝ ╚═════╝
  //
  // mX-01-01: [Metric name]
  //
  // Source: MX_output_filename.csv
  //
  //
  //
  //
  //
  //
  //
  //
  //
  {
    id: "mX-01-01",
    label: {
      en: "Metric label",
      fr: "Étiquette de la métrique",
    },
    variantLabel: null,

    resultsObjectId: "MX_output_filename.csv",
    valueProps: ["column_name"],
    valueFunc: "AVG",
    valueLabelReplacements: {
      // column_name: "Display Name",
    },

    requiredDisaggregationOptions: [],
    postAggregationExpression: null,
    // Example postAggregationExpression:
    // {
    //   ingredientValues: [
    //     { prop: "numerator", func: "SUM" },
    //     { prop: "denominator", func: "SUM" },
    //   ],
    //   // Must be "name = arithmetic" with exactly one "=". Division is
    //   // NULLIF-guarded automatically — don't write NULLIF yourself.
    //   expression: "value = numerator / denominator",
    // },

    formatAs: "percent",
    hide: false,
    importantNotes: null,

    aiDescription: {
      summary: {
        en: "Brief one-sentence description of what this metric measures.",
        fr: "Description brève en une phrase de ce que cette métrique mesure.",
      },
      methodology: {
        en: "How the metric is calculated. Reference the valueFunc and any parameters.",
        fr: "Comment la métrique est calculée.",
      },
      interpretation: {
        en: "What values mean. What is good/bad? What thresholds matter?",
        fr: "Ce que les valeurs signifient.",
      },
      typicalRange: {
        en: "Expected range of values (e.g., '0-5% is good, >10% is concerning').",
        fr: "Plage de valeurs attendue.",
      },
      caveats: null,
      disaggregationGuidance: {
        en: "How to best disaggregate this metric for insights.",
        fr: "Comment désagréger cette métrique pour obtenir des informations.",
      },
    },

    //
    //
    //
    //
    //
    //
    //
    //
    // ██╗   ██╗██╗███████╗    ██████╗ ██████╗ ███████╗███████╗███████╗████████╗███████╗
    // ██║   ██║██║╚══███╔╝    ██╔══██╗██╔══██╗██╔════╝██╔════╝██╔════╝╚══██╔══╝██╔════╝
    // ██║   ██║██║  ███╔╝     ██████╔╝██████╔╝█████╗  ███████╗█████╗     ██║   ███████╗
    // ╚██╗ ██╔╝██║ ███╔╝      ██╔═══╝ ██╔══██╗██╔══╝  ╚════██║██╔══╝     ██║   ╚════██║
    //  ╚████╔╝ ██║███████╗    ██║     ██║  ██║███████╗███████║███████╗   ██║   ███████║
    //   ╚═══╝  ╚═╝╚══════╝    ╚═╝     ╚═╝  ╚═╝╚══════╝╚══════╝╚══════╝   ╚═╝   ╚══════╝
    //
    //
    //
    //
    //
    //
    //
    //

    vizPresets: [
      //
      //
      //
      //
      //
      //
      //
      //
      // ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
      // ┃                                                                         ┃
      // ┃   PRESET: table-by-region                                               ┃
      // ┃   Table by region                                                       ┃
      // ┃                                                                         ┃
      // ┃   Type: table                                                           ┃
      // ┃   Rows: admin_area_2                                                    ┃
      // ┃   Cols: indicator_common_id                                             ┃
      // ┃                                                                         ┃
      // ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
      //
      //
      //
      //
      //
      //
      //
      //
      {
        id: "table-by-region",
        label: {
          en: "Table by region",
          fr: "Tableau par région",
        },
        description: {
          en: "Table showing metric values by indicator and region",
          fr: "Tableau montrant les valeurs par indicateur et région",
        },
        importantNotes: null,
        allowedFilters: ["indicator_common_id", "admin_area_2"],
        createDefaultVisualizationOnInstall: null,
        config: {
          d: {
            type: "table",
            timeseriesGrouping: "period_id",
            valuesDisDisplayOpt: "col",
            disaggregateBy: [
              { disOpt: "indicator_common_id", disDisplayOpt: "col" },
              { disOpt: "admin_area_2", disDisplayOpt: "row" },
            ],
            filterBy: [],
            periodFilter: { filterType: "last_n_months", nMonths: 12 },
          },
          s: {
            content: "lines",
            decimalPlaces: 1,
            // Conditional formatting - uncomment and import CF_* preset:
            // ...CF_90_80,
          },
          t: {
            caption: {
              en: "Table Title",
              fr: "Titre du tableau",
            },
            captionRelFontSize: null,
            subCaption: {
              en: "Subtitle with DATE_RANGE placeholder",
              fr: "Sous-titre avec DATE_RANGE",
            },
            subCaptionRelFontSize: null,
            footnote: {
              en: "Methodology footnote explaining how to interpret the data.",
              fr: "Note de bas de page méthodologique.",
            },
            footnoteRelFontSize: null,
          },
        },
      },

      //
      //
      //
      //
      //
      //
      //
      //
      // ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
      // ┃                                                                         ┃
      // ┃   PRESET: timeseries-by-indicator                                       ┃
      // ┃   Trend over time                                                       ┃
      // ┃                                                                         ┃
      // ┃   Type: timeseries                                                      ┃
      // ┃   Series: indicator_common_id                                           ┃
      // ┃                                                                         ┃
      // ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
      //
      //
      //
      //
      //
      //
      //
      //
      {
        id: "timeseries-by-indicator",
        label: {
          en: "Trend over time",
          fr: "Tendance dans le temps",
        },
        description: {
          en: "Line chart showing trends over time by indicator",
          fr: "Graphique linéaire montrant les tendances par indicateur",
        },
        importantNotes: null,
        allowedFilters: ["indicator_common_id"],
        createDefaultVisualizationOnInstall: null,
        config: {
          d: {
            type: "timeseries",
            timeseriesGrouping: "period_id",
            valuesDisDisplayOpt: "series",
            disaggregateBy: [
              { disOpt: "indicator_common_id", disDisplayOpt: "series" },
            ],
            filterBy: [],
          },
          s: {
            content: "lines",
            decimalPlaces: 1,
          },
          t: {
            caption: {
              en: "Chart Title",
              fr: "Titre du graphique",
            },
            captionRelFontSize: null,
            subCaption: {
              en: "Subtitle DATE_RANGE",
              fr: "Sous-titre DATE_RANGE",
            },
            subCaptionRelFontSize: null,
            footnote: null,
            footnoteRelFontSize: null,
          },
        },
      },
    ],
  },
];
