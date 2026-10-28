import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";
import { CF_DIVERGING_5_10_20 } from "../../.validation/cf_presets.ts";

// One metric for every admin level: observed and expected are additive, so
// the level is a disaggregation choice, not a metric (compare _band_chart.ts,
// where the credible bounds force a metric per level).

type AdminLevel = 2 | 3 | 4;
type DisOpt = VizPreset["config"]["d"]["disaggregateBy"][number]["disOpt"];

const adminArea = (level: AdminLevel): DisOpt => `admin_area_${level}`;

const areaLabel = (level: AdminLevel) => ({
  en: `Admin Area ${level}`,
  fr: `Zone administrative ${level}`,
  pt: `Área administrativa ${level}`,
});

const LAST_QUARTER: VizPreset["config"]["d"]["periodFilter"] = {
  filterType: "last_n_calendar_quarters",
  nQuarters: 1,
};

const FOOTNOTE = {
  en:
    "Negative values indicate fewer services than expected. The expected value is the Bayesian model's estimate from trend and seasonality for the facilities that reported in the month.",
  fr:
    "Les valeurs négatives indiquent moins de services qu'attendu. La valeur attendue est l'estimation du modèle bayésien à partir de la tendance et de la saisonnalité pour les établissements ayant rapporté dans le mois.",
  pt:
    "Os valores negativos indicam menos serviços do que o esperado. O valor esperado é a estimativa do modelo bayesiano a partir da tendência e da sazonalidade para os estabelecimentos que reportaram no mês.",
};

const TEXT_SIZES = {
  captionRelFontSize: null,
  subCaptionRelFontSize: null,
  footnoteRelFontSize: null,
};

// Areas of the level down, indicators across; the children of one parent
// area when a parent is given, chosen as the replicant.
function tablePreset(
  level: AdminLevel,
  parent: 2 | 3 | undefined,
): VizPreset {
  const area = areaLabel(level);
  const parentArea = parent === undefined ? undefined : areaLabel(parent);
  return {
    id: `disruption-table-admin-area-${level}`,
    label: parentArea === undefined
      ? {
        en: `Service disruption by ${area.en}`,
        fr: `Perturbation des services par ${area.fr}`,
        pt: `Perturbação dos serviços por ${area.pt}`,
      }
      : {
        en: `Service disruption by ${area.en} for one ${parentArea.en}`,
        fr:
          `Perturbation des services par ${area.fr} pour une ${parentArea.fr}`,
        pt: `Perturbação dos serviços por ${area.pt} para uma ${parentArea.pt}`,
      },
    description: parentArea === undefined
      ? {
        en:
          `Table showing the percentage difference between observed and expected services, ${area.en} down and indicators across`,
        fr:
          `Tableau montrant la différence en pourcentage entre les services observés et attendus, ${area.fr} en lignes et indicateurs en colonnes`,
        pt:
          `Tabela que mostra a diferença em percentagem entre os serviços observados e esperados, ${area.pt} em linhas e indicadores em colunas`,
      }
      : {
        en:
          `Table showing the percentage difference between observed and expected services, ${area.en} down and indicators across, for a single ${parentArea.en}`,
        fr:
          `Tableau montrant la différence en pourcentage entre les services observés et attendus, ${area.fr} en lignes et indicateurs en colonnes, pour une seule ${parentArea.fr}`,
        pt:
          `Tabela que mostra a diferença em percentagem entre os serviços observados e esperados, ${area.pt} em linhas e indicadores em colunas, para uma única ${parentArea.pt}`,
      },
    importantNotes: null,
    allowedFilters: ["indicator_common_id", adminArea(level)],
    config: {
      d: {
        type: "table",
        valuesDisDisplayOpt: "col",
        disaggregateBy: [
          { disOpt: "indicator_common_id", disDisplayOpt: "col" },
          ...(parent === undefined ? [] : [{
            disOpt: adminArea(parent),
            disDisplayOpt: "replicant" as const,
          }]),
          { disOpt: adminArea(level), disDisplayOpt: "row" },
        ],
        filterBy: [],
        periodFilter: LAST_QUARTER,
        ...(parent === undefined ? {} : { selectedReplicantValue: "" }),
      },
      s: { ...CF_DIVERGING_5_10_20 },
      t: {
        caption: {
          en: `Service disruption by ${area.en}`,
          fr: `Perturbation des services par ${area.fr}`,
          pt: `Perturbação dos serviços por ${area.pt}`,
        },
        subCaption: {
          en:
            "Percentage difference between observed and expected services, DATE_RANGE",
          fr:
            "Différence en pourcentage entre les services observés et attendus, PLAGE_DE_DATES",
          pt:
            "Diferença em percentagem entre os serviços observados e esperados, INTERVALO_DE_DATAS",
        },
        footnote: FOOTNOTE,
        ...TEXT_SIZES,
      },
    },
  };
}

function mapPreset(level: AdminLevel): VizPreset {
  const area = areaLabel(level);
  return {
    id: `disruption-map-admin-area-${level}`,
    label: {
      en: `Service disruption map by ${area.en}`,
      fr: `Carte des perturbations par ${area.fr}`,
      pt: `Mapa das perturbações por ${area.pt}`,
    },
    description: {
      en:
        `Map showing the percentage difference between observed and expected services by ${area.en}, for one indicator`,
      fr:
        `Carte montrant la différence en pourcentage entre les services observés et attendus par ${area.fr}, pour un indicateur`,
      pt:
        `Mapa que mostra a diferença em percentagem entre os serviços observados e esperados por ${area.pt}, para um indicador`,
    },
    importantNotes: null,
    allowedFilters: [],
    config: {
      d: {
        type: "map",
        valuesDisDisplayOpt: "cell",
        disaggregateBy: [
          { disOpt: adminArea(level), disDisplayOpt: "mapArea" },
          { disOpt: "indicator_common_id", disDisplayOpt: "replicant" },
        ],
        filterBy: [],
        periodFilter: LAST_QUARTER,
        selectedReplicantValue: "anc1",
      },
      s: { ...CF_DIVERGING_5_10_20 },
      t: {
        caption: {
          en: "Service disruption",
          fr: "Perturbation des services",
          pt: "Perturbação dos serviços",
        },
        subCaption: {
          en:
            `Percentage difference between observed and expected services by ${area.en}, DATE_RANGE`,
          fr:
            `Différence en pourcentage entre les services observés et attendus par ${area.fr}, PLAGE_DE_DATES`,
          pt:
            `Diferença em percentagem entre os serviços observados e esperados por ${area.pt}, INTERVALO_DE_DATAS`,
        },
        footnote: FOOTNOTE,
        ...TEXT_SIZES,
      },
    },
  };
}

export const vizPresets: VizPreset[] = [
  tablePreset(2, undefined),
  mapPreset(2),
  tablePreset(3, 2),
  mapPreset(3),
  tablePreset(4, 3),
];

export const metric: MetricDefinitionGithub = {
  id: "m11-02-01",
  resultsObjectId: "M11_disruptions_analysis.csv",
  valueProps: ["pct_diff"],
  valueFunc: "identity",
  formatAs: "percent",
  postAggregationExpression: {
    ingredientValues: [
      { prop: "observed", func: "SUM" },
      { prop: "expected", func: "SUM" },
    ],
    expression: "pct_diff = (observed - expected) / expected",
  },
  label: {
    en: "Difference between observed and expected services (%)",
    fr: "Différence entre les services observés et attendus (%)",
    pt: "Diferença entre os serviços observados e esperados (%)",
  },
  variantLabel: null,
  requiredDisaggregationOptions: ["indicator_common_id"],
  valueLabelReplacements: {
    pct_diff: {
      en: "Percent difference",
      fr: "Différence en pourcentage",
      pt: "Diferença em percentagem",
    },
  },
  catalogExpressionEvaluation: null,
  aiDescription: {
    summary: {
      en:
        "Percentage difference between observed services and the Bayesian expected baseline, at any admin level and for any period.",
      fr:
        "Différence en pourcentage entre les services observés et la base de référence bayésienne attendue, à tout niveau administratif et pour toute période.",
      pt:
        "Diferença em percentagem entre os serviços observados e a base de referência bayesiana esperada, a qualquer nível administrativo e para qualquer período.",
    },
    methodology: {
      en:
        "Observed and expected services are summed over the grouping shown, then the difference is taken as a share of expected: (observed - expected) / expected. The expected baseline is the Bayesian disruption model's latent intensity summed over the same facilities that reported each month (see the band-chart metrics).",
      fr:
        "Les services observés et attendus sont sommés sur le regroupement affiché, puis la différence est rapportée à l'attendu : (observé - attendu) / attendu. La base de référence attendue est l'intensité latente du modèle bayésien de perturbation, sommée sur les mêmes établissements ayant rapporté chaque mois (voir les métriques du graphique à bande).",
      pt:
        "Os serviços observados e esperados são somados no agrupamento apresentado e a diferença é expressa como proporção do esperado: (observado - esperado) / esperado. A base de referência esperada é a intensidade latente do modelo bayesiano de perturbação, somada sobre os mesmos estabelecimentos que reportaram em cada mês (ver as métricas do gráfico de banda).",
    },
    interpretation: {
      en:
        "Negative values mean fewer services than expected, positive values more. The preset thresholds at 5%, 10% and 20% either side mark increasing departures; whether a departure lies outside the model's credible interval is shown by the band-chart metrics.",
      fr:
        "Les valeurs négatives signifient moins de services qu'attendu, les positives davantage. Les seuils des préréglages à 5%, 10% et 20% de part et d'autre marquent des écarts croissants ; le graphique à bande montre si un écart sort de l'intervalle de crédibilité du modèle.",
      pt:
        "Os valores negativos significam menos serviços do que o esperado, os positivos mais. Os limiares das predefinições a 5%, 10% e 20% de cada lado marcam desvios crescentes; o gráfico de banda mostra se um desvio sai do intervalo de credibilidade do modelo.",
    },
    typicalRange: {
      en:
        "Mostly within 10% either side; sustained values beyond 20% indicate a disruption or surge to triangulate with other data.",
      fr:
        "Le plus souvent à moins de 10% de part et d'autre ; des valeurs durablement au-delà de 20% indiquent une perturbation ou une hausse à trianguler avec d'autres données.",
      pt:
        "Na maioria dos casos dentro de 10% de cada lado; valores persistentemente acima de 20% indicam uma perturbação ou um aumento a triangular com outros dados.",
    },
    caveats: {
      en:
        "A share of expected, so small expected volumes give large percentages; read alongside the observed volume. Only facilities that reported enter both terms, so the value reflects service-rate change, not reporting gaps.",
      fr:
        "Rapportée à l'attendu, donc de petits volumes attendus donnent de grands pourcentages ; à lire avec le volume observé. Seuls les établissements ayant rapporté entrent dans les deux termes, la valeur reflète donc un changement du taux de service, pas des lacunes de rapportage.",
      pt:
        "É uma proporção do esperado, pelo que volumes esperados pequenos dão percentagens grandes; ler em conjunto com o volume observado. Só os estabelecimentos que reportaram entram nos dois termos, pelo que o valor reflete a mudança da taxa de serviço e não lacunas de reporte.",
    },
    disaggregationGuidance: {
      en:
        "Always by indicator_common_id (required). By admin area at any level for tables and maps; by month, quarter or year over time.",
      fr:
        "Toujours par indicator_common_id (requis). Par zone administrative à tout niveau pour les tableaux et les cartes ; par mois, trimestre ou année dans le temps.",
      pt:
        "Sempre por indicator_common_id (obrigatório). Por área administrativa a qualquer nível para tabelas e mapas; por mês, trimestre ou ano ao longo do tempo.",
    },
  },
  importantNotes: null,
  hide: false,
  vizPresets,
};
