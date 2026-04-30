import type { MetricDefinitionGithub, VizPreset } from "../../.validation/_module_definition_github.ts";
import { CF_01_03, CF_05_10, CF_10_20, CF_80_70, CF_90_80, CF_NEG10_POS10 } from "../../.validation/cf_presets.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "adjustment-table",
    label: {
      en: "Combined adjustment impact table",
      fr: "Tableau d'impact de l'ajustement combiné",
    },
    description: {
      en: "Table showing percent change due to combined outlier and completeness adjustment by indicator and Admin Area 2",
      fr: "Tableau montrant le changement en pourcentage dû à l'ajustement combiné des valeurs aberrantes et de la complétude par indicateur et région",
    },
    createDefaultVisualizationOnInstall: "5337d614-02b8-4de8-abcb-f390d2b7a714",
    allowedFilters: ["indicator_common_id", "admin_area_2"],
    config: {
      d: {
        type: "table",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "col",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "col",
          },
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "row",
          },
        ],
        filterBy: [],
        periodFilter: {
          filterType: "last_n_months",
          nMonths: 12,
        },
      },
      s: {
        content: "lines",
        ...CF_01_03,
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Deviance Due to Incompleteness and Outliers",
          fr: "Déviance due à l'incomplétude et aux valeurs aberrantes",
        },
        subCaption: {
          en: "Percent change in volume due to both outlier and completeness adjustment, DATE_RANGE",
          fr: "Changement en pourcentage du volume dû à l'ajustement combiné des valeurs aberrantes et de la complétude, PLAGE_DE_DATES",
        },
        footnote: {
          en: "Combined adjustment applies both outlier removal and completeness imputation. Outliers are suspiciously high values identified using median absolute deviation. Completeness adjustment imputes missing data from non-reporting facilities. The deviance shows the total change in volume after both adjustments. High levels of deviance indicate substantial data quality issues.",
          fr: "L'ajustement combiné applique à la fois la suppression des valeurs aberrantes et l'imputation de la complétude. Les valeurs aberrantes sont des valeurs suspectement élevées identifiées à l'aide de l'écart absolu médian. L'ajustement de complétude impute les données manquantes des établissements non déclarants. La déviance montre le changement total de volume après les deux ajustements. Des niveaux élevés de déviance indiquent des problèmes substantiels de qualité des données.",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
];

export const metric: MetricDefinitionGithub = {
  id: "m2-01-03",
  resultsObjectId: "M2_adjusted_data.csv",
  label: {
    en: "Percent change in volume due to both outlier and completeness adjustment",
    fr: "Changement en pourcentage du volume dû à l'ajustement combiné des valeurs aberrantes et de la complétude",
  },
  valueProps: ["pct_change"],
  valueFunc: "identity",
  valueLabelReplacements: {
    pct_change: "Percent change",
  },
  postAggregationExpression: {
    ingredientValues: [
      {
        prop: "count_final_none",
        func: "SUM",
      },
      {
        prop: "count_final_both",
        func: "SUM",
      },
    ],
    expression: "pct_change = ABS(count_final_none-count_final_both)/count_final_none",
  },
  requiredDisaggregationOptions: [],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Combined magnitude of change in reported service volumes after both outlier removal and missing data imputation.",
      fr: "Ampleur combinée du changement dans les volumes après suppression des aberrants et imputation des données manquantes.",
    },
    methodology: {
      en: "Calculated as ABS(unadjusted - both_adjusted) / unadjusted. Applies both outlier correction (replacing extreme values) and completeness adjustment (filling missing records) sequentially to produce fully-adjusted estimates.",
      fr: "Calculé comme ABS(non ajusté - ajusté pour les deux) / non ajusté. Applique à la fois la correction des aberrants et l'ajustement de complétude séquentiellement.",
    },
    interpretation: {
      en: "Represents the total impact of data quality corrections on service volumes. Higher percentages indicate that raw data required substantial adjustment. Compare with individual adjustment metrics (m2-01-01, m2-01-02) to understand whether outliers or completeness drove the change.",
      fr: "Représente l'impact total des corrections de qualité des données. Des pourcentages plus élevés indiquent que les données brutes nécessitaient un ajustement substantiel.",
    },
    typicalRange: {
      en: "0-10% indicates minor corrections; 10-25% moderate adjustment; >25% indicates major data quality issues requiring substantial correction.",
      fr: "0-10% indique des corrections mineures; 10-25% ajustement modéré; >25% indique des problèmes majeurs nécessitant correction substantielle.",
    },
    caveats: {
      en: "Combined adjustment effect is not simply additive - outlier and completeness adjustments interact. Large combined changes may indicate the need to verify that adjustments are appropriate rather than over-correcting.",
      fr: "L'effet d'ajustement combiné n'est pas simplement additif - les ajustements interagissent. Les grands changements combinés peuvent indiquer la nécessité de vérifier les ajustements.",
    },
    disaggregationGuidance: {
      en: "Disaggregate by indicator_common_id to identify indicators requiring most adjustment. Time series shows if data quality improves over time (decreasing adjustment percentages). Regional disaggregation reveals geographic patterns in data quality.",
      fr: "Désagréger par indicator_common_id pour identifier les indicateurs nécessitant le plus d'ajustement. Les séries temporelles montrent si la qualité s'améliore.",
    },
  },
  variantLabel: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
