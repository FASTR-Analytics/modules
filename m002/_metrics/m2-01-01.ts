import type { MetricDefinitionGithub, VizPreset } from "../../.validation/_module_definition_github.ts";
import { CF_01_03, CF_05_10, CF_10_20, CF_80_70, CF_90_80, CF_NEG10_POS10 } from "../../.validation/cf_presets.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "adjustment-table",
    label: {
      en: "Outlier adjustment impact table",
      fr: "Tableau d'impact de l'ajustement des valeurs aberrantes",
    },
    description: {
      en: "Table showing percent change due to outlier adjustment by indicator and region",
      fr: "Tableau montrant le changement en pourcentage dû à l'ajustement des valeurs aberrantes par indicateur et région",
    },
    createDefaultVisualizationOnInstall: "e5edce68-369c-498e-a4b0-03ba73d31d6c",
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
          en: "Deviance Due to Outliers",
          fr: "Déviance due aux valeurs aberrantes",
        },
        subCaption: {
          en: "Percent change in volume due to outlier adjustment, DATE_RANGE",
          fr: "Changement en pourcentage du volume dû à l'ajustement des valeurs aberrantes, DATE_RANGE",
        },
        footnote: {
          en: "Outliers are reports which are suspiciously high compared to the usual volume reported by the facility in other months. Outliers are identified by assessing the within-facility variation in monthly reporting for each indicator. Outliers are defined observations which are greater than 10 times the median absolute deviation (MAD) from the monthly median value for the indicator in each time period, OR a value for which the proportional contribution in volume for a facility, indicator, and time period is greater than 80%. Outliers are only identified for indicators where the volume is greater than or equal to the median, the volume is not missing, and the average volume is greater than 100. The deviance is the difference in volume after removing the outlier. High levels of deviance can affect the plausiability of the data.",
          fr: "Les valeurs aberrantes sont des rapports anormalement élevés par rapport au volume habituel déclaré par l'établissement au cours des autres mois. Elles sont identifiées en évaluant la variation intra-établissement des déclarations mensuelles pour chaque indicateur. Les valeurs aberrantes sont définies comme des observations supérieures à 10 fois l'écart absolu médian (MAD) par rapport à la valeur médiane mensuelle de l'indicateur pour chaque période, OU une valeur dont la contribution proportionnelle au volume pour un établissement, indicateur et période est supérieure à 80%. Les valeurs aberrantes ne sont identifiées que pour les indicateurs dont le volume est supérieur ou égal à la médiane, le volume n'est pas manquant, et le volume moyen est supérieur à 100. La déviance est la différence de volume après suppression de la valeur aberrante. Des niveaux élevés de déviance peuvent affecter la plausibilité des données.",
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
  id: "m2-01-01",
  resultsObjectId: "M2_adjusted_data.csv",
  label: {
    en: "Percent change in volume due to outlier adjustment",
    fr: "Changement en pourcentage du volume dû à l'ajustement des valeurs aberrantes",
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
        prop: "count_final_outliers",
        func: "SUM",
      },
    ],
    expression: "pct_change = ABS(count_final_none-count_final_outliers)/count_final_none",
  },
  requiredDisaggregationOptions: [],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Magnitude of change in reported service volumes after removing or adjusting outlier values.",
      fr: "Ampleur du changement dans les volumes de services déclarés après suppression ou ajustement des valeurs aberrantes.",
    },
    methodology: {
      en: "Calculated as ABS(unadjusted - outlier_adjusted) / unadjusted. Outlier adjustment replaces flagged extreme values using rolling mean imputation (6-month centered, forward, or backward windows) or facility-level means.",
      fr: "Calculé comme ABS(non ajusté - ajusté pour aberrants) / non ajusté. L'ajustement des aberrants remplace les valeurs extrêmes par imputation par moyenne mobile.",
    },
    interpretation: {
      en: "Higher percentages indicate that outliers had significant impact on reported volumes. Values above 10% suggest substantial data quality issues that could bias analysis if left unadjusted. Compare across indicators to identify those most affected.",
      fr: "Des pourcentages plus élevés indiquent que les aberrants ont eu un impact significatif. Les valeurs supérieures à 10% suggèrent des problèmes de qualité importants.",
    },
    typicalRange: {
      en: "0-5% indicates minimal outlier impact; 5-15% moderate; >15% suggests major quality issues.",
      fr: "0-5% indique un impact aberrant minimal; 5-15% modéré; >15% suggère des problèmes majeurs.",
    },
    caveats: {
      en: "Percentage change reflects magnitude but not direction. Large changes may be appropriate if outliers were genuine data errors, but could also indicate over-correction if flagged values were valid.",
      fr: "Le changement en pourcentage reflète l'ampleur mais pas la direction. Les grands changements peuvent être appropriés si les aberrants étaient des erreurs.",
    },
    disaggregationGuidance: {
      en: "Disaggregate by indicator_common_id to identify which services are most affected by outlier adjustment. Use admin_area to find regions where outliers have larger impact. Time series analysis shows if adjustment impact changes over time.",
      fr: "Désagréger par indicator_common_id pour identifier quels services sont les plus affectés. Utiliser admin_area pour trouver les régions où l'impact est plus important.",
    },
  },
  variantLabel: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
