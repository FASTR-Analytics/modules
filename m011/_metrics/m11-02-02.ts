import type { MetricDefinitionGithub } from "../../.validation/_module_definition_github.ts";

// The additive pair at any admin level, for custom tables. The band chart is
// on the per-level metrics, which carry the credible bounds.
export const metric: MetricDefinitionGithub = {
  id: "m11-02-02",
  resultsObjectId: "M11_disruptions_analysis.csv",
  valueProps: ["observed", "expected"],
  valueFunc: "SUM",
  formatAs: "number",
  postAggregationExpression: null,
  label: {
    en: "Observed and expected services",
    fr: "Services observés et attendus",
    pt: "Serviços observados e esperados",
  },
  variantLabel: null,
  requiredDisaggregationOptions: ["indicator_common_id"],
  valueLabelReplacements: {
    observed: "Observed services",
    expected: "Expected services (Bayesian)",
  },
  catalogExpressionEvaluation: null,
  aiDescription: {
    summary: {
      en:
        "Observed services and the Bayesian expected baseline, summed to any admin level and period.",
      fr:
        "Services observés et base de référence bayésienne attendue, sommés à tout niveau administratif et pour toute période.",
      pt:
        "Serviços observados e base de referência bayesiana esperada, somados a qualquer nível administrativo e para qualquer período.",
    },
    methodology: {
      en:
        "Both values are sums over facilities of the Bayesian disruption model's inputs and outputs: observed is the reported count, expected is the model's latent intensity for the same facilities that reported each month.",
      fr:
        "Les deux valeurs sont des sommes sur les établissements des entrées et sorties du modèle bayésien de perturbation : l'observé est le nombre rapporté, l'attendu est l'intensité latente du modèle pour les mêmes établissements ayant rapporté chaque mois.",
      pt:
        "Ambos os valores são somas sobre os estabelecimentos das entradas e saídas do modelo bayesiano de perturbação: o observado é a contagem reportada, o esperado é a intensidade latente do modelo para os mesmos estabelecimentos que reportaram em cada mês.",
    },
    interpretation: {
      en:
        "Observed below expected is a shortfall, above is a surplus. The percentage version is the companion metric; whether a gap is outside the credible interval is shown by the band-chart metrics.",
      fr:
        "Un observé inférieur à l'attendu est un déficit, supérieur un excédent. La version en pourcentage est la métrique compagne ; le graphique à bande montre si un écart sort de l'intervalle de crédibilité.",
      pt:
        "Um observado inferior ao esperado é um défice, superior um excedente. A versão em percentagem é a métrica companheira; o gráfico de banda mostra se um desvio sai do intervalo de credibilidade.",
    },
    typicalRange: {
      en:
        "Volumes in the units of the indicator; expected tracks observed closely outside disruptions.",
      fr:
        "Volumes dans les unités de l'indicateur ; l'attendu suit de près l'observé hors perturbation.",
      pt:
        "Volumes nas unidades do indicador; o esperado acompanha de perto o observado fora das perturbações.",
    },
    caveats: {
      en:
        "Only facilities that reported enter both values, so a reporting gap lowers both and does not open a gap between them.",
      fr:
        "Seuls les établissements ayant rapporté entrent dans les deux valeurs ; une lacune de rapportage abaisse les deux sans creuser d'écart entre elles.",
      pt:
        "Só os estabelecimentos que reportaram entram nos dois valores, pelo que uma lacuna de reporte baixa ambos sem abrir um desvio entre eles.",
    },
    disaggregationGuidance: {
      en:
        "Always by indicator_common_id (required). By admin area at any level; by month, quarter or year over time.",
      fr:
        "Toujours par indicator_common_id (requis). Par zone administrative à tout niveau ; par mois, trimestre ou année dans le temps.",
      pt:
        "Sempre por indicator_common_id (obrigatório). Por área administrativa a qualquer nível; por mês, trimestre ou ano ao longo do tempo.",
    },
  },
  importantNotes: null,
  hide: false,
  vizPresets: [],
};
