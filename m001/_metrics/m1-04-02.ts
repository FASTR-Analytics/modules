import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { CF_80_70 } from "../../.validation/cf_presets.ts";

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
    id: "mean-dqa-table",
    label: {
      en: "Mean DQA score table",
      fr: "Tableau du score EQD moyen",
    },
    description: {
      en: "Table showing mean DQA scores by region and year",
      fr: "Tableau montrant les scores EQD moyens par région et année",
    },
    createDefaultVisualizationOnInstall: "4dc02c21-29da-4a01-9812-469deedaaac8",
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
        ...CF_80_70,
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Mean DQA score",
          fr: "Score EQD moyen",
        },
        subCaption: {
          en: "Average data quality score across facility-months",
          fr: "Score moyen de qualité des données à travers les mois-établissements",
        },
        footnote: {
          en: "Items included in the DQA score include: No missing data for 1) OPD, 2) Penta1, and 3) ANC1, where available; No outliers for 4) OPD, 5) Penta1, and 6) ANC1, where available; Consistent reporting between 7) Penta1/Penta3, 8) ANC1/ANC4, 9)BCG/Delivery, where available.",
          fr: "Les éléments inclus dans le score EQD comprennent : Pas de données manquantes pour 1) OPD, 2) Penta1 et 3) ANC1, lorsque disponibles ; Pas de valeurs aberrantes pour 4) OPD, 5) Penta1 et 6) ANC1, lorsque disponibles ; Déclaration cohérente entre 7) Penta1/Penta3, 8) ANC1/ANC4, 9) BCG/Accouchement, lorsque disponibles.",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
  ///////////////////////////////////////
  //  __       __                      //
  // /  \     /  |                     //
  // $$  \   /$$ |  ______    ______   //
  // $$$  \ /$$$ | /      \  /      \  //
  // $$$$  /$$$$ | $$$$$$  |/$$$$$$  | //
  // $$ $$ $$/$$ | /    $$ |$$ |  $$ | //
  // $$ |$$$/ $$ |/$$$$$$$ |$$ |__$$ | //
  // $$ | $/  $$ |$$    $$ |$$    $$/  //
  // $$/      $$/  $$$$$$$/ $$$$$$$/   //
  //                        $$ |       //
  //                        $$ |       //
  //                        $$/        //
  //                                   //
  ///////////////////////////////////////
  {
    id: "mean-dqa-map",
    label: {
      en: "Mean DQA score map",
      fr: "Carte du score EQD moyen",
    },
    description: {
      en: "Map showing mean DQA scores by region",
      fr: "Carte montrant les scores EQD moyens par région",
    },
    createDefaultVisualizationOnInstall: "0730bd53-f951-4286-83ef-5b299395912c",
    allowedFilters: [],
    config: {
      d: {
        type: "map",
        valuesDisDisplayOpt: "mapArea",
        disaggregateBy: [
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "mapArea",
          },
        ],
        filterBy: [],
        periodFilter: {
          filterType: "last_n_months",
          nMonths: 12,
        },
      },
      s: {
        ...CF_80_70,
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Mean DQA score",
          fr: "Score EQD moyen",
        },
        subCaption: {
          en: "Average data quality score across facility-months, DATE_RANGE",
          fr: "Score moyen de qualité des données à travers les mois-établissements, DATE_RANGE",
        },
        footnote: {
          en: "Items included in the DQA score include: No missing data for 1) OPD, 2) Penta1, and 3) ANC1, where available; No outliers for 4) OPD, 5) Penta1, and 6) ANC1, where available; Consistent reporting between 7) Penta1/Penta3, 8) ANC1/ANC4, 9)BCG/Delivery, where available.",
          fr: "Les éléments inclus dans le score EQD comprennent : Pas de données manquantes pour 1) OPD, 2) Penta1 et 3) ANC1, lorsque disponibles ; Pas de valeurs aberrantes pour 4) OPD, 5) Penta1 et 6) ANC1, lorsque disponibles ; Déclaration cohérente entre 7) Penta1/Penta3, 8) ANC1/ANC4, 9) BCG/Accouchement, lorsque disponibles.",
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
  id: "m1-04-02",
  resultsObjectId: "M1_output_dqa.csv",
  valueProps: ["dqa_mean"],
  valueFunc: "AVG",
  valueLabelReplacements: {
    dqa_mean: "Data quality score across facilities",
  },
  label: {
    en: "Average data quality score across facilities",
    fr: "Score moyen de qualité des données à travers les établissements",
  },
  requiredDisaggregationOptions: [],
  formatAs: "percent",
  aiDescription: {
    summary: {
      en: "Average composite data quality score across all facilities.",
      fr: "Score moyen composite de qualité des données à travers tous les établissements.",
    },
    methodology: {
      en: "AVG of continuous dqa_mean score. Combines completeness, outlier, and consistency assessments.",
      fr: "Moyenne du score dqa_mean continu. Combine les évaluations de complétude, valeurs aberrantes et cohérence.",
    },
    interpretation: {
      en: "Higher values indicate better overall data quality. Use alongside m1-04-01 to understand both average performance and threshold compliance.",
      fr: "Des valeurs plus élevées indiquent une meilleure qualité globale des données.",
    },
    typicalRange: {
      en: "0.7-1.0 is good; 0.5-0.7 moderate; <0.5 indicates significant issues.",
      fr: "0.7-1.0 est bon; 0.5-0.7 modéré; <0.5 indique des problèmes significatifs.",
    },
    disaggregationGuidance: {
      en: "Disaggregate by admin_area for regional comparison. Use time series to track improvement over time.",
      fr: "Désagréger par zone administrative pour comparaison régionale.",
    },
    caveats: null,
  },
  variantLabel: null,
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
