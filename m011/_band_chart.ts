import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../.validation/_module_definition_github.ts";

// One band-chart metric per admin level, because the credible bounds are the
// level's own: with USE_POSTERIOR_DRAWS they are quantiles of that level's
// summed draws and do not add up from a finer level. Everything the four
// metrics share is built here so they differ only where their grain does.

export type AdminLevel = 2 | 3 | 4;
export type Level = 1 | AdminLevel;

type DisOpt = VizPreset["config"]["d"]["disaggregateBy"][number]["disOpt"];
type DisaggregateBy = VizPreset["config"]["d"]["disaggregateBy"];

const adminArea = (level: AdminLevel): DisOpt => `admin_area_${level}`;

const areaLabel = (level: Level) => ({
  en: `Admin Area ${level}`,
  fr: `Zone administrative ${level}`,
});

const CONFIG_S: VizPreset["config"]["s"] = {
  content: "lines",
  specialDisruptionsChartV2: true,
};

const CONFIG_T: VizPreset["config"]["t"] = {
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
};

const DESCRIPTION_TAIL = {
  en: "green where observed rises above the band, red where it falls below.",
  fr: "vert où l'observé dépasse la bande, rouge où il tombe en dessous.",
};

function preset(p: {
  id: string;
  label: { en: string; fr: string };
  layout: { en: string; fr: string };
  allowedFilters: DisOpt[];
  disaggregateBy: DisaggregateBy;
  selectedReplicantValue?: string;
}): VizPreset {
  return {
    id: p.id,
    label: p.label,
    description: {
      en: `Observed (solid) and expected (dashed) lines with the 95% credible band, ${p.layout.en}; ${DESCRIPTION_TAIL.en}`,
      fr: `Lignes observée (pleine) et attendue (pointillée) avec la bande de crédibilité à 95%, ${p.layout.fr} ; ${DESCRIPTION_TAIL.fr}`,
    },
    importantNotes: null,
    allowedFilters: p.allowedFilters,
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "period_id",
        valuesDisDisplayOpt: "series",
        disaggregateBy: p.disaggregateBy,
        filterBy: [],
        ...(p.selectedReplicantValue === undefined
          ? {}
          : { selectedReplicantValue: p.selectedReplicantValue }),
      },
      s: CONFIG_S,
      t: CONFIG_T,
    },
  };
}

// One panel per indicator, the whole area the metric covers.
export function byIndicatorPreset(): VizPreset {
  return preset({
    id: "observed-expected-ci",
    label: {
      en: "Disruptions chart by indicator",
      fr: "Graphique des perturbations par indicateur",
    },
    layout: { en: "one panel per indicator", fr: "un panneau par indicateur" },
    allowedFilters: ["indicator_common_id"],
    disaggregateBy: [{ disOpt: "indicator_common_id", disDisplayOpt: "cell" }],
  });
}

// Indicators down, the level's areas across.
export function gridPreset(level: 2 | 3 | 4): VizPreset {
  const area = areaLabel(level);
  return preset({
    id: "observed-expected-ci",
    label: {
      en: `Disruptions chart by indicator and ${area.en}`,
      fr: `Graphique des perturbations par indicateur et ${area.fr}`,
    },
    layout: {
      en: `one panel per indicator and ${area.en}`,
      fr: `un panneau par indicateur et ${area.fr}`,
    },
    allowedFilters: ["indicator_common_id", adminArea(level)],
    disaggregateBy: [
      { disOpt: "indicator_common_id", disDisplayOpt: "row" },
      { disOpt: adminArea(level), disDisplayOpt: "col" },
    ],
  });
}

// The same grid for the children of one parent area, chosen as the replicant.
export function childrenOfParentPreset(level: 3 | 4): VizPreset {
  const parent = (level - 1) as 2 | 3;
  const area = areaLabel(level);
  const parentArea = areaLabel(parent);
  return preset({
    id:
      `observed-expected-ci-single-admin-area-${parent}-multiple-admin-area-${level}`,
    label: {
      en: `Disruptions chart by indicator and ${area.en} for one ${parentArea.en}`,
      fr: `Graphique des perturbations par indicateur et ${area.fr} pour une ${parentArea.fr}`,
    },
    layout: {
      en: `one panel per indicator and ${area.en}, for a single ${parentArea.en}`,
      fr: `un panneau par indicateur et ${area.fr}, pour une seule ${parentArea.fr}`,
    },
    allowedFilters: ["indicator_common_id", adminArea(level)],
    disaggregateBy: [
      { disOpt: "indicator_common_id", disDisplayOpt: "row" },
      { disOpt: adminArea(parent), disDisplayOpt: "replicant" },
      { disOpt: adminArea(level), disDisplayOpt: "col" },
    ],
    selectedReplicantValue: "",
  });
}

// One area, chosen as the replicant, one panel per indicator.
export function singleAreaPreset(level: 2 | 3 | 4): VizPreset {
  const area = areaLabel(level);
  return preset({
    id: `observed-expected-ci-single-admin-area-${level}`,
    label: {
      en: `Disruptions chart for one ${area.en}`,
      fr: `Graphique des perturbations pour une ${area.fr}`,
    },
    layout: {
      en: `one panel per indicator, for a single ${area.en}`,
      fr: `un panneau par indicateur, pour une seule ${area.fr}`,
    },
    allowedFilters: ["indicator_common_id"],
    disaggregateBy: [
      { disOpt: "indicator_common_id", disDisplayOpt: "cell" },
      { disOpt: adminArea(level), disDisplayOpt: "replicant" },
    ],
    selectedReplicantValue: "",
  });
}

// One indicator, chosen as the replicant, one panel per area of the level.
export function singleIndicatorPreset(level: 2 | 3 | 4): VizPreset {
  const area = areaLabel(level);
  return preset({
    id: "observed-expected-ci-single-indicator",
    label: {
      en: `Disruptions chart for one indicator by ${area.en}`,
      fr: `Graphique des perturbations pour un indicateur par ${area.fr}`,
    },
    layout: {
      en: `one panel per ${area.en}, for a single indicator`,
      fr: `un panneau par ${area.fr}, pour un seul indicateur`,
    },
    allowedFilters: [adminArea(level)],
    disaggregateBy: [
      { disOpt: "indicator_common_id", disDisplayOpt: "replicant" },
      { disOpt: adminArea(level), disDisplayOpt: "cell" },
    ],
    selectedReplicantValue: "",
  });
}

export function bandChartMetric(p: {
  id: string;
  level: Level;
  vizPresets: VizPreset[];
}): MetricDefinitionGithub {
  const grain = p.level === 1
    ? { en: "indicator_common_id (required)", fr: "indicator_common_id (requis)" }
    : {
      en: `indicator_common_id and admin_area_${p.level} (both required)`,
      fr: `indicator_common_id et admin_area_${p.level} (tous deux requis)`,
    };
  return {
    id: p.id,
    resultsObjectId: `M11_disruptions_analysis_admin_area_${p.level}.csv`,
    variantLabel: p.level === 1
      ? { en: "National", fr: "National" }
      : { en: `Admin area ${p.level}`, fr: `Zone administrative ${p.level}` },
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
    requiredDisaggregationOptions: p.level === 1
      ? ["indicator_common_id"]
      : ["indicator_common_id", adminArea(p.level)],
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
        en: "Observed below ppi_lwr indicates a real service decline among reporting facilities (deficit). Observed above ppi_upr indicates a service surge (surplus). Values inside the band are within expected variation given trend, seasonality, and facility heterogeneity, not signal.",
        fr: "Observé sous ppi_lwr indique une vraie baisse de service parmi les rapporteurs (déficit). Observé au-dessus de ppi_upr indique une augmentation (surplus). Les valeurs dans la bande sont dans la variation attendue.",
      },
      typicalRange: {
        en: "Varies by service and country. Mean reporting (mean_p) typically 60-95%.",
        fr: "Varie selon le service et le pays. Complétude moyenne (mean_p) typiquement 60-95%.",
      },
      caveats: {
        en: "The expected baseline is conditional on the same facilities that reported each month, so it isolates service-rate change from reporting-coverage change. The default marginal-sum credible interval is slightly conservative at fine admin levels; USE_POSTERIOR_DRAWS=TRUE gives calibrated bands at every level.",
        fr: "La base de référence attendue est conditionnelle aux mêmes établissements rapportant chaque mois, isolant le changement de taux de service du changement de complétude. L'intervalle de crédibilité par somme marginale, par défaut, est légèrement conservateur aux niveaux administratifs fins ; USE_POSTERIOR_DRAWS=TRUE donne des bandes calibrées à chaque niveau.",
      },
      disaggregationGuidance: {
        en: `Always disaggregate by ${grain.en}. Subnational comparison reveals geographic concentration of disruptions; time series reveals onset and recovery.`,
        fr: `Toujours désagréger par ${grain.fr}. La comparaison sous-nationale révèle la concentration géographique ; les séries temporelles révèlent début et reprise.`,
      },
    },
    importantNotes: null,
    catalogExpressionEvaluation: null,
    vizPresets: p.vizPresets,
  };
}
