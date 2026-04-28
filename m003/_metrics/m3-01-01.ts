import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "volume-monthly",
    label: {
      en: "Service volume over time (monthly)",
      fr: "Volume de services dans le temps (mensuel)",
    },
    description: {
      en: "Line chart showing monthly service volume by indicator",
      fr: "Graphique linéaire montrant le volume de services mensuel par indicateur",
    },
    createDefaultVisualizationOnInstall: "45f2bcd8-879d-4423-a4b0-a84127e168bf",
    allowedFilters: ["indicator_common_id"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "cell",
          },
        ],
        filterBy: [],
        valuesFilter: ["count_final_both"],
      },
      s: {
        content: "lines",
      },
      t: {
        caption: {
          en: "Service utilization over time",
          fr: "Utilisation des services dans le temps",
        },
        subCaption: {
          en: "DATE_RANGE",
          fr: "DATE_RANGE",
        },
        footnote: {
          en: "Yearly volume is adjusted for both outliers and completeness.",
          fr: "Le volume annuel est ajusté à la fois pour les valeurs aberrantes et la complétude.",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
  {
    id: "volume-quarterly",
    label: {
      en: "Volume quarterly change",
      fr: "Variation trimestrielle du volume",
    },
    description: {
      en: "Bar chart showing quarterly volume with quarter-on-quarter change",
      fr: "Diagramme à barres montrant le volume trimestriel avec variation trimestre par trimestre",
    },
    createDefaultVisualizationOnInstall: null,
    allowedFilters: ["indicator_common_id"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "quarter_id",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "row",
          },
        ],
        filterBy: [],
        periodFilter: {
          filterType: "last_n_months",
          nMonths: 12,
        },
        valuesFilter: ["count_final_outliers"],
      },
      s: {
        specialBarChart: true,
        specialBarChartDataLabels: "all-values",
      },
      t: {
        caption: {
          en: "Service volume by quarter & quarter-on-quarter change",
          fr: "Volume de services par trimestre et variation trimestre par trimestre",
        },
        subCaption: {
          en: "DATE_RANGE",
          fr: "DATE_RANGE",
        },
        footnote: {
          en: "Service volume is adjusted for outliers.",
          fr: "Le volume de services est ajusté pour les valeurs aberrantes.",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
  {
    id: "volume-annual",
    label: {
      en: "Volume annual change",
      fr: "Variation annuelle du volume",
    },
    description: {
      en: "Bar chart showing annual volume with year-on-year change",
      fr: "Diagramme à barres montrant le volume annuel avec variation d'une année sur l'autre",
    },
    createDefaultVisualizationOnInstall: null,
    allowedFilters: ["indicator_common_id"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "row",
          },
        ],
        filterBy: [],
        valuesFilter: ["count_final_outliers"],
      },
      s: {
        specialBarChart: true,
        specialBarChartDataLabels: "all-values",
      },
      t: {
        caption: {
          en: "Service volume by year & year-on-year change",
          fr: "Volume de services par année et variation d'une année sur l'autre",
        },
        subCaption: {
          en: "DATE_RANGE",
          fr: "DATE_RANGE",
        },
        footnote: {
          en: "Service volume is adjusted for outliers.",
          fr: "Le volume de services est ajusté pour les valeurs aberrantes.",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
  {
    id: "volume-subnational",
    label: {
      en: "Volume annual change by region",
      fr: "Variation annuelle du volume par région",
    },
    description: {
      en: "Bar chart showing annual volume change by indicator and admin area",
      fr: "Diagramme à barres montrant la variation annuelle du volume par indicateur et zone administrative",
    },
    createDefaultVisualizationOnInstall: "20658bc8-2b24-4adc-8090-407c6e34f22a",
    allowedFilters: ["indicator_common_id", "admin_area_2"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "row",
          },
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "col",
          },
        ],
        filterBy: [],
        valuesFilter: ["count_final_outliers"],
      },
      s: {
        scale: 1.7,
        specialBarChart: true,
        specialBarChartDataLabels: "all-values",
      },
      t: {
        caption: {
          en: "Service volume by year & year-on-year change",
          fr: "Volume de services par année et variation d'une année sur l'autre",
        },
        subCaption: {
          en: "DATE_RANGE",
          fr: "DATE_RANGE",
        },
        footnote: {
          en: "Yearly volume is adjusted for outliers.",
          fr: "Le volume annuel est ajusté pour les valeurs aberrantes.",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
  {
    id: "dq-comparison",
    label: {
      en: "Data quality adjustment comparison",
      fr: "Comparaison des ajustements de qualité des données",
    },
    description: {
      en: "Line chart comparing volume under different adjustment scenarios",
      fr: "Graphique linéaire comparant le volume selon différents scénarios d'ajustement",
    },
    createDefaultVisualizationOnInstall: "508f17cc-fbfd-4585-a2e8-8242234898c3",
    allowedFilters: ["indicator_common_id"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "col",
          },
        ],
        filterBy: [],
      },
      s: {
        scale: 1.8,
        colorScale: "custom",
        customSeriesStyles: [
          {
            color: "#00897b",
            lineStyle: "solid",
            strokeWidth: 5,
          },
          {
            color: "#757575",
            lineStyle: "solid",
            strokeWidth: 5,
          },
          {
            color: "#8e24aa",
            lineStyle: "solid",
            strokeWidth: 5,
          },
          {
            color: "#7cb342",
            lineStyle: "solid",
            strokeWidth: 5,
          },
        ],
      },
      t: {
        caption: {
          en: "Change in volume due to data quality adjustments",
          fr: "Variation du volume due aux ajustements de qualité des données",
        },
        subCaption: {
          en: "DATE_RANGE",
          fr: "DATE_RANGE",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnote: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
];

export const metric: MetricDefinitionGithub = {
  id: "m3-01-01",
  resultsObjectId: "M3_service_utilization.csv",
  valueProps: [
    "count_final_none",
    "count_final_outliers",
    "count_final_completeness",
    "count_final_both",
  ],
  valueFunc: "SUM",
  valueLabelReplacements: {
    count_final_none: "Number of services reported",
    count_final_outliers: "Number of services after outlier adjustment",
    count_final_completeness:
      "Number of services after completeness adjustment",
    count_final_both:
      "Number of services after both outlier and completeness adjustment",
  },
  label: {
    en: "Number of services reported, by adjustment type",
    fr: "Nombre de services déclarés, par type d'ajustement",
  },
  requiredDisaggregationOptions: ["indicator_common_id"],
  formatAs: "number",
  aiDescription: {
    summary: {
      en: "Total service volumes across four data quality adjustment scenarios: unadjusted, outlier-adjusted, completeness-adjusted, and both adjustments combined.",
      fr: "Volumes totaux de services à travers quatre scénarios d'ajustement de qualité: non ajusté, ajusté pour aberrants, ajusté pour complétude, et les deux ajustements combinés.",
    },
    methodology: {
      en: "SUM of service counts under each adjustment type. Four values are presented: (1) raw reported volumes, (2) outlier-adjusted volumes, (3) completeness-adjusted volumes, (4) fully-adjusted volumes with both corrections applied.",
      fr: "Somme des comptes de services sous chaque type d'ajustement. Quatre valeurs sont présentées: (1) volumes bruts déclarés, (2) volumes ajustés pour aberrants, (3) volumes ajustés pour complétude, (4) volumes totalement ajustés avec les deux corrections appliquées.",
    },
    interpretation: {
      en: "Comparing the four adjustment types reveals the impact of data quality corrections on service totals. Large differences between unadjusted and adjusted values indicate significant data quality issues. Users can select which adjustment scenario to use for subsequent analysis based on their data quality tolerance.",
      fr: "La comparaison des quatre types d'ajustement révèle l'impact des corrections de qualité des données sur les totaux de services. De grandes différences entre les valeurs non ajustées et ajustées indiquent des problèmes de qualité des données importants. Les utilisateurs peuvent sélectionner le scénario d'ajustement à utiliser pour l'analyse ultérieure en fonction de leur tolérance à la qualité des données.",
    },
    typicalRange: {
      en: "Varies by service type and dataset quality. Completeness adjustment typically increases totals; outlier adjustment may increase or decrease totals.",
      fr: "Varie selon le type de service et la qualité du jeu de données. L'ajustement de complétude augmente généralement les totaux; l'ajustement des valeurs aberrantes peut augmenter ou diminuer les totaux.",
    },
    caveats: {
      en: "The appropriate adjustment type depends on analytical goals. Conservative analyses may prefer minimal adjustment; comprehensive coverage estimates may require full adjustment. Maternal/neonatal/under-5 deaths are excluded from all adjustments.",
      fr: "Le type d'ajustement approprié dépend des objectifs analytiques. Les analyses conservatrices peuvent préférer un ajustement minimal; les estimations de couverture complète peuvent nécessiter un ajustement complet. Les décès maternels/néonatals/moins de 5 ans sont exclus de tous les ajustements.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id (required) to see adjustment impact per service. Time series reveals if data quality improves over time. Regional disaggregation shows geographic variation in adjustment needs.",
      fr: "Toujours désagréger par indicator_common_id (requis) pour voir l'impact de l'ajustement par service. Les séries temporelles révèlent si la qualité des données s'améliore au fil du temps. La désagrégation régionale montre la variation géographique des besoins d'ajustement.",
    },
  },
  importantNotes: {
    en: "The intention of the first four visualization presets (volume-monthly, volume-quarterly, volume-annual, volume-subnational) is that only one of the values is selected: count_final_none, count_final_outliers, count_final_completeness, or count_final_both. You should use valuesFilter to select it (i.e. filter so it is only one value). To determine which to use, run get_module_settings for module id 'm003' and see what value is set for 'Count variable to use for visualization', which should be the default value that you use (although the user may ask for a different one, which you can use instead). If the user wants to compare different adjustments, you should use preset dq-comparison.",
    fr: "L'intention des quatre premiers préréglages de visualisation (volume-monthly, volume-quarterly, volume-annual, volume-subnational) est qu'une seule des valeurs soit sélectionnée : count_final_none, count_final_outliers, count_final_completeness, ou count_final_both. Vous devez utiliser valuesFilter pour la sélectionner (c'est-à-dire filtrer pour n'avoir qu'une seule valeur). Pour déterminer laquelle utiliser, exécutez get_module_settings pour le module 'm003' et consultez la valeur définie pour 'Count variable to use for visualization', qui devrait être la valeur par défaut que vous utilisez (bien que l'utilisateur puisse en demander une différente, que vous pouvez utiliser à la place). Si l'utilisateur souhaite comparer différents ajustements, vous devez utiliser le préréglage dq-comparison.",
  },
  variantLabel: null,
  postAggregationExpression: null,
  hide: false,
  vizPresets,
};
