import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

const vizPresets: VizPreset[] = [
  {
    id: "observed-expected-ci",
    label: {
      en: "Disruptions chart with 95% credible band",
      fr: "Graphique des perturbations avec bande de crédibilité à 95%",
    },
    description: {
      en: "Observed (solid) and expected (dashed) lines per indicator and admin area, by month, over a grey 95% credible band — green shading where observed rises above the band, red where it falls below.",
      fr: "Lignes observée (pleine) et attendue (pointillée) par indicateur et zone administrative, par mois, sur une bande de crédibilité grise à 95% — ombrage vert où l'observé dépasse la bande, rouge où il tombe en dessous.",
    },
    createDefaultVisualizationOnInstall: null,
    allowedFilters: ["indicator_common_id", "admin_area_2"],
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "period_id",
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
      },
      s: {
        content: "lines",
        specialDisruptionsChartV2: true,
      },
      t: {
        caption: {
          en: "Observed vs expected services with 95% credible interval",
          fr: "Services observés vs attendus avec intervalle de crédibilité à 95%",
        },
        subCaption: {
          en: "DATE_RANGE",
          fr: "PLAGE_DE_DATES",
        },
        footnote: null,
        captionRelFontSize: null,
        subCaptionRelFontSize: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
];

export const metric: MetricDefinitionGithub = {
  id: "m11-01-02",
  resultsObjectId: "M11_disruptions_analysis_admin_area_2.csv",
  variantLabel: {
    en: "Admin area 2",
    fr: "Zone administrative 2",
  },
  postAggregationExpression: null,
  hide: false,
  valueProps: ["observed", "expected", "ppi_lwr", "ppi_upr"],
  valueFunc: "SUM",
  valueLabelReplacements: {
    observed: "Observed services",
    expected: "Expected services (Bayesian)",
    ppi_lwr: "Lower 95% credible bound",
    ppi_upr: "Upper 95% credible bound",
  },
  label: {
    en: "Bayesian disruption detection — observed vs expected with 95% CI",
    fr: "Détection bayésienne — observé vs attendu avec IC à 95%",
  },
  requiredDisaggregationOptions: ["indicator_common_id", "admin_area_2"],
  formatAs: "number",
  aiDescription: {
    summary: {
      en: "Observed service volumes vs Bayesian expected baseline with 95% credible interval, isolating real service-rate disruptions from reporting-coverage gaps.",
      fr: "Volumes de services observés vs base de référence bayésienne avec intervalle de crédibilité à 95%, isolant les vraies perturbations de service des écarts de complétude.",
    },
    methodology: {
      en: "Two-part Bayesian model (author: Mustapha Wasseja) fit with INLA. Part 1 (Reporting): a Bernoulli model on the full facility-by-month grid predicts the probability that each facility reports. Part 2 (Service): a NegBin model on reporters only estimates latent intensity mu(f,t). The expected baseline is the sum of mu over the same facilities that reported each month (apples-to-apples). Deviations outside the 95% credible interval are flagged.",
      fr: "Modèle bayésien en deux parties (auteur: Mustapha Wasseja) ajusté avec INLA. Partie 1 (Complétude): un modèle Bernoulli sur la grille complète établissement-mois prédit la probabilité de rapport. Partie 2 (Service): un modèle NegBin sur les rapporteurs estime l'intensité latente mu(f,t). La base de référence est la somme de mu sur les mêmes établissements rapporteurs chaque mois.",
    },
    interpretation: {
      en: "Observed below ppi_lwr indicates a real service decline among reporting facilities (deficit). Observed above ppi_upr indicates a service surge (surplus). Values inside the band are within expected variation given trend, seasonality, and facility heterogeneity — not signal.",
      fr: "Observé sous ppi_lwr indique une vraie baisse de service parmi les rapporteurs (déficit). Observé au-dessus de ppi_upr indique une augmentation (surplus). Les valeurs dans la bande sont dans la variation attendue.",
    },
    typicalRange: {
      en: "Varies by service and country. Mean reporting (mean_p) typically 60-95%.",
      fr: "Varie selon le service et le pays. Complétude moyenne (mean_p) typiquement 60-95%.",
    },
    caveats: {
      en: "The expected baseline is conditional on the same facilities that reported each month, so it isolates service-rate change from reporting-coverage change. The marginal-sum CI is approximate at fine admin levels — set USE_POSTERIOR_DRAWS=TRUE for properly calibrated bands at admin_area_3.",
      fr: "La base de référence attendue est conditionnelle aux mêmes établissements rapportant chaque mois, isolant le changement de taux de service du changement de complétude. L'IC par somme marginale est approximatif aux niveaux administratifs fins; mettre USE_POSTERIOR_DRAWS=TRUE pour des bandes correctement calibrées au niveau admin_area_3.",
    },
    disaggregationGuidance: {
      en: "Always disaggregate by indicator_common_id (required). Subnational comparison reveals geographic concentration of disruptions; time series reveals onset and recovery.",
      fr: "Toujours désagréger par indicator_common_id (requis). La comparaison sous-nationale révèle la concentration géographique; les séries temporelles révèlent début et reprise.",
    },
  },
  importantNotes: null,
  vizPresets,
};
