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
      pt: "Perturbações e excedentes (várias Áreas administrativas 2, vários indicadores)",
    },
    description: {
      en: "Area chart showing actual vs expected service volume, with multiple Admin Areas 2 and multiple indicators",
      fr: "Graphique en aires montrant le volume de services réel vs attendu, avec plusieurs Zones administratives 2 et plusieurs indicateurs",
      pt: "Gráfico de áreas que mostra o volume de serviços real vs esperado, com várias Áreas administrativas 2 e vários indicadores",
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
          pt: "Perturbações e excedentes no volume de serviços, a nível subnacional",
        },
        subCaption: {
          en: "DATE_RANGE",
          fr: "PLAGE_DE_DATES",
          pt: "INTERVALO_DE_DATAS",
        },
        footnote: {
          en: "This graph quantifies changes in service volume compared to historical trends and accounting for seasonality. These signals should be triangulated to other data and contextual knowledge to determine if the results are an artifact of data quality. Unexpected volume changes are estimated by comparing the observed volume to the expected volume based on historical trends and seasonality. Previous large unexpected changes in the historical data are removed. This analysis is an interrupted time series regression with facility-level fixed effects.",
          fr: "Ce graphique quantifie les changements du volume de services par rapport aux tendances historiques et en tenant compte de la saisonnalité. Ces signaux doivent être triangulés avec d'autres données et connaissances contextuelles pour déterminer si les résultats sont un artefact de la qualité des données. Les changements de volume inattendus sont estimés en comparant le volume observé au volume attendu basé sur les tendances historiques et la saisonnalité. Les grands changements inattendus précédents dans les données historiques sont supprimés. Cette analyse est une régression de séries temporelles interrompues avec effets fixes au niveau de l'établissement.",
          pt: "Este gráfico quantifica as alterações no volume de serviços em comparação com as tendências históricas e tendo em conta a sazonalidade. Estes sinais devem ser triangulados com outros dados e conhecimento contextual para determinar se os resultados são um artefacto da qualidade dos dados. As alterações de volume inesperadas são estimadas comparando o volume observado com o volume esperado com base nas tendências históricas e na sazonalidade. As grandes alterações inesperadas anteriores nos dados históricos são removidas. Esta análise é uma regressão de séries temporais interrompidas com efeitos fixos ao nível do estabelecimento.",
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
      pt: "Perturbações e excedentes (uma única Área administrativa 2, vários indicadores)",
    },
    description: {
      en: "Area chart showing actual vs expected service volume, for a single Admin Area 2 and multiple indicators",
      fr: "Graphique en aires montrant le volume de services réel vs attendu, pour une unique Zone administrative 2 et plusieurs indicateurs",
      pt: "Gráfico de áreas que mostra o volume de serviços real vs esperado, para uma única Área administrativa 2 e vários indicadores",
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
        subCaption: {
          en: "DATE_RANGE",
          fr: "PLAGE_DE_DATES",
          pt: "INTERVALO_DE_DATAS",
        },
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
      pt: "Perturbações e excedentes (um único indicador, várias Áreas administrativas 2)",
    },
    description: {
      en: "Area chart showing actual vs expected service volume, for a single indicator and multiple Admin Area 2",
      fr: "Graphique en aires montrant le volume de services réel vs attendu, pour un unique indicateur et plusieurs Zones administratives 2",
      pt: "Gráfico de áreas que mostra o volume de serviços real vs esperado, para um único indicador e várias Áreas administrativas 2",
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
        subCaption: {
          en: "DATE_RANGE",
          fr: "PLAGE_DE_DATES",
          pt: "INTERVALO_DE_DATAS",
        },
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
    en: "Disruptions and surpluses",
    fr: "Perturbations et excédents",
    pt: "Perturbações e excedentes",
  },
  variantLabel: {
    en: "Admin area 2",
    fr: "Zone administrative 2",
    pt: "Área administrativa 2",
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
      pt: "Comparação dos volumes de serviços reais notificados com os volumes esperados previstos pelo modelo ao nível da área administrativa 2 (província/estado).",
    },
    methodology: {
      en: "SUM of actual vs expected counts from area-specific robust regression models. Expected volumes account for local time trends and seasonal patterns. Only periods exceeding the difference threshold are shown.",
      fr: "Somme des comptes réels vs attendus des modèles de régression robustes spécifiques à la zone. Les volumes attendus tiennent compte des tendances temporelles locales et des modèles saisonniers. Seules les périodes dépassant le seuil de différence sont affichées.",
      pt: "SOMA das contagens reais vs esperadas a partir de modelos de regressão robustos específicos da área. Os volumes esperados têm em conta as tendências temporais locais e os padrões sazonais. São apresentados apenas os períodos que excedem o limiar de diferença.",
    },
    interpretation: {
      en: "Enables subnational identification of service disruptions. Compare across admin areas to identify geographic hotspots of service delivery problems. Areas with persistent negative deviations need targeted support.",
      fr: "Permet l'identification sous-nationale des perturbations de services. Comparer entre les zones administratives pour identifier les points chauds géographiques de problèmes de prestation de services. Les zones avec des écarts négatifs persistants nécessitent un soutien ciblé.",
      pt: "Permite a identificação subnacional das perturbações dos serviços. Comparar entre áreas administrativas para identificar os focos geográficos de problemas na prestação de serviços. As áreas com desvios negativos persistentes necessitam de apoio direcionado.",
    },
    typicalRange: {
      en: "Expected to closely match actual in stable regions. Deviations indicate local disruptions, data issues, or demand changes.",
      fr: "Attendu pour correspondre étroitement au réel dans les régions stables. Les écarts indiquent des perturbations locales, des problèmes de données ou des changements de la demande.",
      pt: "Espera-se que corresponda de perto ao real nas regiões estáveis. Os desvios indicam perturbações locais, problemas de dados ou alterações da procura.",
    },
    caveats: {
      en: "Smaller geographic areas may have more volatile patterns, making model predictions less reliable. Consider aggregating to higher levels if area-specific models perform poorly.",
      fr: "Les zones géographiques plus petites peuvent avoir des modèles plus volatils, rendant les prédictions du modèle moins fiables. Considérer l'agrégation à des niveaux supérieurs si les modèles spécifiques aux zones fonctionnent mal.",
      pt: "As áreas geográficas mais pequenas podem ter padrões mais voláteis, tornando as previsões do modelo menos fiáveis. Considerar a agregação a níveis superiores se os modelos específicos das áreas tiverem um desempenho fraco.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id and admin_area_2 (both required). Time series reveals when and where disruptions occurred. Map visualization effectively shows geographic distribution of service gaps.",
      fr: "Toujours désagréger par indicator_common_id et admin_area_2 (tous deux requis). Les séries temporelles révèlent quand et où les perturbations se sont produites. La visualisation cartographique montre efficacement la distribution géographique des écarts de services.",
      pt: "Desagregar sempre por indicator_common_id e admin_area_2 (ambos obrigatórios). As séries temporais revelam quando e onde ocorreram as perturbações. A visualização cartográfica mostra eficazmente a distribuição geográfica das lacunas de serviços.",
    },
  },
  postAggregationExpression: null,
  importantNotes: null,
  hide: false,
  vizPresets,
};
