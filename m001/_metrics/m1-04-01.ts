import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [
  ///////////////////////////////////////////////
  //  ________         __        __            //
  // /        |       /  |      /  |           //
  // $$$$$$$$/______  $$ |____  $$ |  ______   //
  //    $$ | /      \ $$      \ $$ | /      \  //
  //    $$ | $$$$$$  |$$$$$$$  |$$ |/$$$$$$  | //
  //    $$ | /    $$ |$$ |  $$ |$$ |$$    $$ | //
  //    $$ |/$$$$$$$ |$$ |__$$ |$$ |$$$$$$$$/  //
  //    $$ |$$    $$ |$$    $$/ $$ |$$       | //
  //    $$/  $$$$$$$/ $$$$$$$/  $$/  $$$$$$$/  //
  //                                           //
  ///////////////////////////////////////////////
  {
    id: "dqa-score-table",
    label: {
      en: "Overall DQA score table",
      fr: "Tableau du score EQD global",
    },
    description: {
      en: "Table showing DQA scores by region and year",
      fr: "Tableau montrant les scores EQD par région et année",
    },
    createDefaultVisualizationOnInstall: "d46e1957-09dd-41c3-b7dc-b4409da23bbe",
    allowedFilters: ["admin_area_2"],
    config: {
      d: {
        type: "table",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "col",
        disaggregateBy: [
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "row",
          },
          {
            disOpt: "year",
            disDisplayOpt: "col",
          },
        ],
        filterBy: [],
      },
      s: {
        content: "lines",
        cfMode: "thresholds",
        cfThresholdCutoffs: [0.7, 0.8],
        cfThresholdBuckets: [
          {
            color: "#F18989",
          },
          {
            color: "#F6D982",
          },
          {
            color: "#68C690",
          },
        ],
        cfThresholdDirection: "higher-is-better",
        cfThresholdNoDataColor: "#ffffff",
        cfScalePaletteKind: "preset",
        cfScalePalettePreset: "",
        cfScaleCustomFrom: "",
        cfScaleCustomMid: "",
        cfScaleCustomTo: "",
        cfScaleReverse: false,
        cfScaleSteps: 0,
        cfScaleDomainKind: "auto",
        cfScaleDomainMin: 0,
        cfScaleDomainMax: 1,
        cfScaleNoDataColor: "",
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Overall DQA score",
          fr: "Score EQD global",
        },
        subCaption: {
          en: "Percentage of facility-months with adequate data quality over time",
          fr: "Pourcentage de mois-établissements avec une qualité des données adéquate dans le temps",
        },
        footnote: {
          en: "Adequate data quality is defined as: 1) No missing data or outliers for OPD, Penta1, and ANC1, where available 2) Consistent reporting between Penta1/Penta3 and ANC1/ANC4.",
          fr: "La qualité adéquate des données est définie comme : 1) Pas de données manquantes ou de valeurs aberrantes pour OPD, Penta1 et ANC1, lorsque disponibles 2) Déclaration cohérente entre Penta1/Penta3 et ANC1/ANC4.",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
];

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
  id: "m1-04-01",
  resultsObjectId: "M1_output_dqa.csv",
  valueProps: ["dqa_score"],
  valueFunc: "AVG",
  valueLabelReplacements: {
    dqa_score: "Binary variable indicating adequate data quality",
  },
  label: {
    en: "Proportion of facilities with adequate data quality",
    fr: "Proportion d'établissements avec une qualité des données adéquate",
  },
  requiredDisaggregationOptions: [],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Proportion of facilities meeting the composite data quality assessment threshold.",
      fr: "Proportion d'établissements atteignant le seuil d'évaluation composite de la qualité des données.",
    },
    methodology: {
      en: "AVG of binary dqa_score based on composite assessment of completeness, outliers, and consistency.",
      fr: "Moyenne du score DQA binaire basé sur une évaluation composite de la complétude, des valeurs aberrantes et de la cohérence.",
    },
    interpretation: {
      en: "Higher values indicate more facilities with trustworthy data. Facilities below threshold may need additional support or data verification.",
      fr: "Des valeurs plus élevées indiquent plus d'établissements avec des données fiables.",
    },
    typicalRange: {
      en: "70-100% is acceptable; <70% indicates widespread quality issues.",
      fr: "70-100% est acceptable; <70% indique des problèmes de qualité généralisés.",
    },
    caveats: {
      en: "Composite score weights may need adjustment based on local context. Consider individual components for detailed diagnosis.",
      fr: "Les poids du score composite peuvent nécessiter un ajustement selon le contexte local.",
    },
    disaggregationGuidance: {
      en: "Disaggregate by admin_area to identify regions needing quality improvement support. Use facility_type to see if certain facility levels have more challenges.",
      fr: "Désagréger par zone administrative pour identifier les régions nécessitant un soutien.",
    },
  },
  variantLabel: null,
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
