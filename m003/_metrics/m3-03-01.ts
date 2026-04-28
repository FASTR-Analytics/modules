import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //  ________  __                                                        __                      //
  // /        |/  |                                                      /  |                     //
  // $$$$$$$$/ $$/  _____  ____    ______    _______   ______    ______  $$/   ______    _______  //
  //    $$ |   /  |/     \/    \  /      \  /       | /      \  /      \ /  | /      \  /       | //
  //    $$ |   $$ |$$$$$$ $$$$  |/$$$$$$  |/$$$$$$$/ /$$$$$$  |/$$$$$$  |$$ |/$$$$$$  |/$$$$$$$/  //
  //    $$ |   $$ |$$ | $$ | $$ |$$    $$ |$$      \ $$    $$ |$$ |  $$/ $$ |$$    $$ |$$      \  //
  //    $$ |   $$ |$$ | $$ | $$ |$$$$$$$$/  $$$$$$  |$$$$$$$$/ $$ |      $$ |$$$$$$$$/  $$$$$$  | //
  //    $$ |   $$ |$$ | $$ | $$ |$$       |/     $$/ $$       |$$ |      $$ |$$       |/     $$/  //
  //    $$/    $$/ $$/  $$/  $$/  $$$$$$$/ $$$$$$$/   $$$$$$$/ $$/       $$/  $$$$$$$/ $$$$$$$/   //
  //                                                                                              //
  //////////////////////////////////////////////////////////////////////////////////////////////////
  {
    id: "disruption-chart",
    label: {
      en: "Disruptions and surpluses (multiple Admin Area 2, multiple indicators)",
      fr: "Perturbations et excédents (plusieurs Zones administratives 2, plusieurs indicateurs)",
    },
    description: {
      en: "Area chart showing actual vs expected service volume, with multiple Admin Areas 2 and multiple indicators",
      fr: "Graphique en aires montrant le volume de services réel vs attendu, avec plusieurs Zones administratives 2 et plusieurs indicateurs",
    },
    createDefaultVisualizationOnInstall: "e1916b10-433a-4b19-b376-491a66b81f11",
    allowedFilters: ["indicator_common_id", "admin_area_2"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "series",
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
      },
      s: {
        specialDisruptionsChart: true,
      },
      t: {
        caption: {
          en: "Disruptions and surpluses in service volume, sub-nationally",
          fr: "Perturbations et excédents du volume de services, au niveau sous-national",
        },
        subCaption: {
          en: "DATE_RANGE",
          fr: "DATE_RANGE",
        },
        footnote: {
          en: "This graph quantifies changes in service volume compared to historical trends and accounting for seasonality. These signals should be triangulated to other data and contextual knowledge to determine if the results are an artifact of data quality. Unexpected volume changes are estimated by comparing the observed volume to the expected volume based on historical trends and seasonality. Previous large unexpected changes in the historical data are removed. This analysis is an interrupted time series regression with facility-level fixed effects.",
          fr: "Ce graphique quantifie les changements du volume de services par rapport aux tendances historiques et en tenant compte de la saisonnalité. Ces signaux doivent être triangulés avec d'autres données et connaissances contextuelles pour déterminer si les résultats sont un artefact de la qualité des données. Les changements de volume inattendus sont estimés en comparant le volume observé au volume attendu basé sur les tendances historiques et la saisonnalité. Les grands changements inattendus précédents dans les données historiques sont supprimés. Cette analyse est une régression de séries temporelles interrompues avec effets fixes au niveau de l'établissement.",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
  //////////////////////////////////////////////////////////////////////////////////////////
  //   ______   __                      __                   ______    ______    ______   //
  //  /      \ /  |                    /  |                 /      \  /      \  /      \  //
  // /$$$$$$  |$$/  _______    ______  $$ |  ______        /$$$$$$  |/$$$$$$  |/$$$$$$  | //
  // $$ \__$$/ /  |/       \  /      \ $$ | /      \       $$ |__$$ |$$ |__$$ |$$____$$ | //
  // $$      \ $$ |$$$$$$$  |/$$$$$$  |$$ |/$$$$$$  |      $$    $$ |$$    $$ | /    $$/  //
  //  $$$$$$  |$$ |$$ |  $$ |$$ |  $$ |$$ |$$    $$ |      $$$$$$$$ |$$$$$$$$ |/$$$$$$/   //
  // /  \__$$ |$$ |$$ |  $$ |$$ \__$$ |$$ |$$$$$$$$/       $$ |  $$ |$$ |  $$ |$$ |_____  //
  // $$    $$/ $$ |$$ |  $$ |$$    $$ |$$ |$$       |      $$ |  $$ |$$ |  $$ |$$       | //
  //  $$$$$$/  $$/ $$/   $$/  $$$$$$$ |$$/  $$$$$$$/       $$/   $$/ $$/   $$/ $$$$$$$$/  //
  //                         /  \__$$ |                                                   //
  //                         $$    $$/                                                    //
  //                          $$$$$$/                                                     //
  //                                                                                      //
  //////////////////////////////////////////////////////////////////////////////////////////
  {
    id: "disruption-chart-single-admin-area-2",
    label: {
      en: "Disruptions and surpluses (single Admin Area 2, multiple indicators)",
      fr: "Perturbations et excédents (unique Zone administrative 2, plusieurs indicateurs)",
    },
    description: {
      en: "Area chart showing actual vs expected service volume, for a single Admin Area 2 and multiple indicators",
      fr: "Graphique en aires montrant le volume de services réel vs attendu, pour une unique Zone administrative 2 et plusieurs indicateurs",
    },
    allowedFilters: ["indicator_common_id"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "series",
        selectedReplicantValue: "",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "cell",
          },
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "replicant",
          },
        ],
        filterBy: [],
      },
      s: {
        specialDisruptionsChart: true,
      },
      t: {
        caption: null,
        captionRelFontSize: null,
        subCaption: null,
        subCaptionRelFontSize: null,
        footnote: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
    createDefaultVisualizationOnInstall: "dee25984-afe7-4643-a823-fd0e26b3cbc0",
  },
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //   ______   __                      __                  ______                  __  __                        __                          //
  //  /      \ /  |                    /  |                /      |                /  |/  |                      /  |                         //
  // /$$$$$$  |$$/  _______    ______  $$ |  ______        $$$$$$/  _______    ____$$ |$$/   _______   ______   _$$ |_     ______    ______   //
  // $$ \__$$/ /  |/       \  /      \ $$ | /      \         $$ |  /       \  /    $$ |/  | /       | /      \ / $$   |   /      \  /      \  //
  // $$      \ $$ |$$$$$$$  |/$$$$$$  |$$ |/$$$$$$  |        $$ |  $$$$$$$  |/$$$$$$$ |$$ |/$$$$$$$/  $$$$$$  |$$$$$$/   /$$$$$$  |/$$$$$$  | //
  //  $$$$$$  |$$ |$$ |  $$ |$$ |  $$ |$$ |$$    $$ |        $$ |  $$ |  $$ |$$ |  $$ |$$ |$$ |       /    $$ |  $$ | __ $$ |  $$ |$$ |  $$/  //
  // /  \__$$ |$$ |$$ |  $$ |$$ \__$$ |$$ |$$$$$$$$/        _$$ |_ $$ |  $$ |$$ \__$$ |$$ |$$ \_____ /$$$$$$$ |  $$ |/  |$$ \__$$ |$$ |       //
  // $$    $$/ $$ |$$ |  $$ |$$    $$ |$$ |$$       |      / $$   |$$ |  $$ |$$    $$ |$$ |$$       |$$    $$ |  $$  $$/ $$    $$/ $$ |       //
  //  $$$$$$/  $$/ $$/   $$/  $$$$$$$ |$$/  $$$$$$$/       $$$$$$/ $$/   $$/  $$$$$$$/ $$/  $$$$$$$/  $$$$$$$/    $$$$/   $$$$$$/  $$/        //
  //                         /  \__$$ |                                                                                                       //
  //                         $$    $$/                                                                                                        //
  //                          $$$$$$/                                                                                                         //
  //                                                                                                                                          //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  {
    id: "disruption-chart-single-indicator",
    label: {
      en: "Disruptions and surpluses (single indicator, multiple Admin Area 2)",
      fr: "Perturbations et excédents (unique indicateur, plusieurs Zones administratives 2)",
    },
    description: {
      en: "Area chart showing actual vs expected service volume, for a single indicator and multiple Admin Area 2",
      fr: "Graphique en aires montrant le volume de services réel vs attendu, pour un unique indicateur et plusieurs Zones administratives 2",
    },
    allowedFilters: ["admin_area_2"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "series",
        selectedReplicantValue: "anc1",
        disaggregateBy: [
          {
            disOpt: "indicator_common_id",
            disDisplayOpt: "replicant",
          },
          {
            disOpt: "admin_area_2",
            disDisplayOpt: "cell",
          },
        ],
        filterBy: [],
      },
      s: {
        specialDisruptionsChart: true,
      },
      t: {
        caption: null,
        captionRelFontSize: null,
        subCaption: null,
        subCaptionRelFontSize: null,
        footnote: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
    createDefaultVisualizationOnInstall: "d95a382e-557c-44da-92e9-63a7abd1c7df",
  },
];

export const metric: MetricDefinitionGithub = {
  id: "m3-03-01",
  resultsObjectId: "M3_disruptions_analysis_admin_area_2.csv",
  label: {
    en: "Actual vs expected service volume",
    fr: "Volume de services réel vs attendu",
  },
  variantLabel: {
    en: "Admin area 2",
    fr: "Zone administrative 2",
  },
  valueProps: ["count_sum", "count_expected_if_above_diff_threshold"],
  valueFunc: "SUM",
  valueLabelReplacements: {
    count_sum: "Actual service volume",
    count_expected_if_above_diff_threshold: "Expected service volume",
  },
  requiredDisaggregationOptions: ["indicator_common_id", "admin_area_2"],
  formatAs: "number",
  aiDescription: {
    summary: {
      en: "Comparison of actual reported service volumes against model-predicted expected volumes at admin area 2 (province/state) level.",
      fr: "Comparaison des volumes de services réels contre volumes attendus au niveau de la zone administrative 2 (province/état).",
    },
    methodology: {
      en: "SUM of actual vs expected counts from area-specific robust regression models. Expected volumes account for local time trends and seasonal patterns. Only periods exceeding the difference threshold are shown.",
      fr: "Somme des comptes réels vs attendus des modèles de régression robustes spécifiques à la zone. Les volumes attendus tiennent compte des tendances temporelles locales et des modèles saisonniers. Seules les périodes dépassant le seuil de différence sont affichées.",
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
