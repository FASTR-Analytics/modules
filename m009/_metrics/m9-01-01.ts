import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "iceh-equiplot",
    label: {
      en: "ICEH equiplot",
      fr: "Équigraphique ICEH",
    },
    description: {
      en: "Coverage by subgroup as an equiplot: one row per indicator, one coloured dot per subgroup. Switch the stratifier with the replicant selector.",
      fr: "Couverture par sous-groupe sous forme d'équigraphique : une ligne par indicateur, un point coloré par sous-groupe. Changez le stratificateur avec le sélecteur de réplicant.",
    },
    createDefaultVisualizationOnInstall: "8cf07c58-5de6-4b48-b361-55b33b51de37",
    allowedFilters: ["iceh_indicator", "strat", "level", "year"],
    config: {
      d: {
        type: "chart",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "col",
        disaggregateBy: [
          {
            disOpt: "iceh_indicator",
            disDisplayOpt: "indicator",
          },
          {
            disOpt: "strat",
            disDisplayOpt: "replicant",
          },
          {
            disOpt: "level",
            disDisplayOpt: "series",
          },
          {
            disOpt: "year",
            disDisplayOpt: "cell",
          },
        ],
        filterBy: [],
        selectedReplicantValue: "wealth_quintiles",
      },
      s: {
        horizontal: true,
        content: "points",
      },
      t: {
        caption: {
          en: "ICEH survey estimates",
          fr: "Estimations d'enquête ICEH",
        },
        subCaption: {
          en: "Estimates by equity stratifier, DATE_RANGE",
          fr: "Estimations par stratificateur d'équité, PLAGE_DE_DATES",
        },
        footnote: {
          en: "",
          fr: "",
        },
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
  {
    id: "iceh-coverage-table",
    label: {
      en: "ICEH coverage table",
      fr: "Tableau de couverture ICEH",
    },
    description: {
      en: "Table of ICEH coverage estimates with indicators as rows and subgroup levels as columns. Switch the stratifier (wealth, area, education, region, ...) with the replicant selector.",
      fr: "Tableau des estimations de couverture ICEH avec les indicateurs en lignes et les niveaux de sous-groupe en colonnes. Changez le stratificateur (richesse, zone, éducation, région, ...) avec le sélecteur de réplicant.",
    },
    allowedFilters: ["iceh_indicator", "strat", "level", "year"],
    createDefaultVisualizationOnInstall: "e7b3c1d4-5a82-4f9e-b6c3-2d8f1a7e90b5",
    config: {
      d: {
        type: "table",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "col",
        disaggregateBy: [
          { disOpt: "iceh_indicator", disDisplayOpt: "row" },
          { disOpt: "strat", disDisplayOpt: "replicant" },
          { disOpt: "level", disDisplayOpt: "col" },
          { disOpt: "year", disDisplayOpt: "colGroup" },
        ],
        filterBy: [],
        selectedReplicantValue: "wealth_quintiles",
      },
      s: {
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "ICEH coverage",
          fr: "Couverture ICEH",
        },
        captionRelFontSize: null,
        subCaption: {
          en: "Coverage by equity stratifier, by indicator, DATE_RANGE",
          fr: "Couverture par stratificateur d'équité, par indicateur, PLAGE_DE_DATES",
        },
        subCaptionRelFontSize: null,
        footnote: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
  {
    id: "iceh-coverage-trend",
    label: {
      en: "ICEH coverage trend",
      fr: "Tendance de couverture ICEH",
    },
    description: {
      en: "Coverage of a selected indicator over survey years, one line per wealth quintile (watch the rich–poor gap change over time). Switch the indicator with the replicant selector.",
      fr: "Couverture d'un indicateur sélectionné au fil des années d'enquête, une ligne par quintile de richesse. Changez l'indicateur avec le sélecteur de réplicant.",
    },
    allowedFilters: ["iceh_indicator", "strat", "level", "year"],
    createDefaultVisualizationOnInstall: "a3f8e1d7-6c24-4b9e-8a51-2f7d9c3b6e08",
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "series",
        disaggregateBy: [
          { disOpt: "iceh_indicator", disDisplayOpt: "replicant" },
          { disOpt: "level", disDisplayOpt: "series" },
        ],
        filterBy: [{ disOpt: "strat", values: ["wealth_quintiles"] }],
        selectedReplicantValue: "vdpt",
      },
      s: {
        content: "lines",
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Coverage trend",
          fr: "Tendance de couverture",
        },
        captionRelFontSize: null,
        subCaption: {
          en: "By wealth quintile, DATE_RANGE",
          fr: "Par quintile de richesse, PLAGE_DE_DATES",
        },
        subCaptionRelFontSize: null,
        footnote: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
];

export const metric: MetricDefinitionGithub = {
  id: "m9-01-01",
  hide: false,
  resultsObjectId: "M9_iceh_data.csv",
  valueProps: ["estimate"],
  valueFunc: "identity",
  valueLabelReplacements: {},
  label: {
    en: "ICEH survey estimate",
    fr: "Estimation d'enquête ICEH",
  },
  requiredDisaggregationOptions: ["iceh_indicator", "level", "year"],
  formatAs: "percent",
  postAggregationExpression: null,
  aiDescription: {
    summary: {
      en: "ICEH Retriever data: standardized estimates from DHS, MICS, and other nationally representative household surveys. Covers 100+ RMNCH+N indicators across multiple years, with equity stratifications (wealth, education, residence, etc.). Use for population-level coverage and prevalence—complements HMIS facility data.",
      fr: "Données ICEH Retriever : estimations standardisées provenant des enquêtes DHS, MICS et autres enquêtes ménages représentatives. Couvre plus de 100 indicateurs SRMNEA sur plusieurs années, avec stratifications d'équité (richesse, éducation, résidence, etc.). Utilisez pour la couverture et prévalence au niveau population—complémente les données HMIS.",
    },
    methodology: {
      en: "Data from the ICEH Retriever (International Center for Equity in Health, Federal University of Pelotas). Harmonized estimates from Demographic and Health Surveys (DHS), Multiple Indicator Cluster Surveys (MICS), and comparable nationally representative surveys. All indicators standardized for cross-country and temporal comparability. Values are survey-weighted estimates with associated standard errors and sample sizes.",
      fr: "Données du ICEH Retriever (Centre International pour l'Équité en Santé, Université Fédérale de Pelotas). Estimations harmonisées des Enquêtes Démographiques et de Santé (EDS), Enquêtes par Grappes à Indicateurs Multiples (MICS) et enquêtes nationales comparables. Tous les indicateurs sont standardisés pour la comparabilité entre pays et dans le temps. Les valeurs sont des estimations pondérées avec erreurs-types et tailles d'échantillon.",
    },
    interpretation: {
      en: "Estimates represent population-level coverage or prevalence from household surveys—different from HMIS facility-reported data. Survey data captures actual population health behaviors (e.g., children fully vaccinated) rather than services delivered. Compare across stratifiers to identify equity gaps. Standard error indicates precision; smaller samples have wider uncertainty.",
      fr: "Les estimations représentent la couverture ou prévalence au niveau population issues d'enquêtes ménages—différent des données HMIS déclarées par les établissements. Les enquêtes capturent les comportements de santé réels (ex: enfants complètement vaccinés) plutôt que les services fournis. Comparez entre stratificateurs pour identifier les écarts d'équité. L'erreur-type indique la précision.",
    },
    typicalRange: {
      en: "Most indicators 0-100% (coverage/prevalence). Mortality rates per 1,000 live births. Fertility rates per 1,000 women.",
      fr: "La plupart des indicateurs 0-100% (couverture/prévalence). Taux de mortalité pour 1 000 naissances vivantes. Taux de fécondité pour 1 000 femmes.",
    },
    caveats: {
      en: "NA = estimate suppressed (sample size <25, or <125 for fertility, <250 for mortality). SS = small sample warning. Survey years may not align with HMIS periods. Recall periods vary by indicator (e.g., births in last 2-5 years). Subnational estimates may have wide confidence intervals.",
      fr: "NA = estimation supprimée (taille d'échantillon <25, ou <125 pour fécondité, <250 pour mortalité). SS = avertissement petit échantillon. Les années d'enquête peuvent ne pas correspondre aux périodes HMIS. Les périodes de rappel varient selon l'indicateur. Les estimations infranationales peuvent avoir de larges intervalles de confiance.",
    },
    disaggregationGuidance: {
      en: "Data structure: each row is one indicator × stratifier × level × year. Stratifiers (strat): national, area (urban/rural), wealth_quintiles, wealth_deciles, womans_education, sex, subnational_unit, etc. Levels: values within each stratifier (e.g., Q1-Q5 for wealth quintiles, urban/rural for area). For equity analysis, compare Q1 vs Q5 (wealth), rural vs urban (area), or none vs secondary+ (education). Filter by strat to focus analysis.",
      fr: "Structure des données : chaque ligne est un indicateur × stratificateur × niveau × année. Stratificateurs (strat) : national, zone (urbain/rural), quintiles_richesse, déciles_richesse, éducation_femme, sexe, unité_infranationale, etc. Niveaux : valeurs au sein de chaque stratificateur (ex: Q1-Q5 pour quintiles de richesse). Pour l'analyse d'équité, comparez Q1 vs Q5 (richesse), rural vs urbain (zone), ou aucune vs secondaire+ (éducation).",
    },
  },
  variantLabel: null,
  importantNotes: null,
  vizPresets,
};
