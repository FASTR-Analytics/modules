import type { MetricDefinitionGithub, VizPreset } from "../../.validation/_module_definition_github.ts";
import { CF_01_03, CF_05_10, CF_10_20, CF_80_70, CF_90_80, CF_NEG10_POS10 } from "../../.validation/cf_presets.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "adjustment-table",
    label: {
      en: "Completeness adjustment impact table",
      fr: "Tableau d'impact de l'ajustement de complétude",
    },
    description: {
      en: "Table showing percent change due to completeness adjustment by indicator and region",
      fr: "Tableau montrant le changement en pourcentage dû à l'ajustement de complétude par indicateur et région",
    },
    createDefaultVisualizationOnInstall: "b4750223-9ffd-43f6-958b-0ba9c0412df4",
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
          en: "Deviance Due to Incompleteness",
          fr: "Déviance due à l'incomplétude",
        },
        subCaption: {
          en: "Percent change in volume due to completeness adjustment, DATE_RANGE",
          fr: "Changement en pourcentage du volume dû à l'ajustement de complétude, DATE_RANGE",
        },
        footnote: {
          en: "Completeness is defined as the percentage of reporting facilities each month out of the total number of facilities expected to report. A facility is expected to report if it has reported any volume for each indicator anytime within a year. The deviance is the difference in volume after imputing incomplete data. High levels of deviance can affect the plausiability of the data.",
          fr: "La complétude est définie comme le pourcentage d'établissements déclarants chaque mois par rapport au nombre total d'établissements censés déclarer. Un établissement est censé déclarer s'il a déclaré un volume pour chaque indicateur à tout moment au cours de l'année. La déviance est la différence de volume après imputation des données incomplètes. Des niveaux élevés de déviance peuvent affecter la plausibilité des données.",
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
  id: "m2-01-02",
  resultsObjectId: "M2_adjusted_data.csv",
  label: {
    en: "Percent change in volume due to completeness adjustment",
    fr: "Changement en pourcentage du volume dû à l'ajustement de complétude",
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
        prop: "count_final_completeness",
        func: "SUM",
      },
    ],
    expression: "pct_change = ABS(count_final_none-count_final_completeness)/count_final_none",
  },
  requiredDisaggregationOptions: [],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Magnitude of change in reported service volumes after imputing missing facility reports.",
      fr: "Ampleur du changement dans les volumes de services déclarés après imputation des rapports d'établissement manquants.",
    },
    methodology: {
      en: "Calculated as ABS(unadjusted - completeness_adjusted) / unadjusted. Completeness adjustment fills missing facility-period records using rolling mean imputation to account for non-reporting facilities.",
      fr: "Calculé comme ABS(non ajusté - ajusté pour complétude) / non ajusté. L'ajustement de complétude remplit les enregistrements manquants par imputation.",
    },
    interpretation: {
      en: "Higher percentages indicate significant missing data that affects totals. Values above 20% suggest incomplete reporting that could substantially underestimate service coverage. This adjustment increases volumes by filling gaps.",
      fr: "Des pourcentages plus élevés indiquent des données manquantes significatives. Les valeurs supérieures à 20% suggèrent une déclaration incomplète.",
    },
    typicalRange: {
      en: "0-10% indicates good reporting completeness; 10-30% moderate gaps; >30% indicates major reporting issues.",
      fr: "0-10% indique une bonne complétude; 10-30% lacunes modérées; >30% indique des problèmes de déclaration majeurs.",
    },
    caveats: {
      en: "Completeness adjustment assumes missing facilities have similar service volumes to reporting facilities. If non-reporting facilities systematically differ (e.g., closed or low-functioning), imputation may over- or under-estimate totals.",
      fr: "L'ajustement de complétude suppose que les établissements manquants ont des volumes similaires. Si les établissements non déclarants diffèrent systématiquement, l'imputation peut surestimer.",
    },
    disaggregationGuidance: {
      en: "Disaggregate by indicator_common_id to see which services have most missing data. Use time periods to identify when reporting completeness deteriorated. Admin_area disaggregation reveals geographic reporting patterns.",
      fr: "Désagréger par indicator_common_id pour voir quels services ont le plus de données manquantes. Utiliser les périodes pour identifier quand la complétude s'est détériorée.",
    },
  },
  variantLabel: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
