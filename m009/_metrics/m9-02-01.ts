import type {
  MetricDefinitionGithub,
  VizPreset,
} from "../../.validation/_module_definition_github.ts";

export const vizPresets: VizPreset[] = [
  {
    id: "iceh-inequality-table",
    label: {
      en: "ICEH inequality measures table",
      fr: "Tableau des mesures d'inégalité ICEH",
      pt: "Tabela de medidas de desigualdade ICEH",
    },
    description: {
      en: "Table of wealth-inequality summary measures (Ratio, Difference, CIX, SII) with indicators as rows and the four measures as columns. Filtered to wealth quintiles; switch the survey year with the replicant selector.",
      fr: "Tableau des mesures synthétiques d'inégalité de richesse (Ratio, Différence, CIX, SII) avec les indicateurs en lignes et les quatre mesures en colonnes. Filtré sur les quintiles de richesse ; changez l'année d'enquête avec le sélecteur de réplicant.",
      pt: "Tabela de medidas-resumo de desigualdade de riqueza (Rácio, Diferença, CIX, SII) com os indicadores em linhas e as quatro medidas em colunas. Filtrada para quintis de riqueza; altere o ano de inquérito com o seletor de replicante.",
    },
    allowedFilters: ["iceh_indicator", "strat", "year"],
    createDefaultVisualizationOnInstall: "c4a9f1e2-7b3d-4e6a-8f12-9d05b7c3a6e1",
    config: {
      d: {
        type: "table",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "col",
        periodFilter: {
          filterType: "last_calendar_year",
        },
        disaggregateBy: [
          { disOpt: "iceh_indicator", disDisplayOpt: "row" },
          { disOpt: "strat", disDisplayOpt: "replicant" },
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
          en: "Wealth inequality measures",
          fr: "Mesures d'inégalité de richesse",
          pt: "Medidas de desigualdade de riqueza",
        },
        captionRelFontSize: null,
        subCaption: {
          en: "Ratio (Q5/Q1), Difference (Q5−Q1, pp), CIX and SII by indicator, DATE_RANGE",
          fr: "Ratio (Q5/Q1), Différence (Q5−Q1, pp), CIX et SII par indicateur, PLAGE_DE_DATES",
          pt: "Rácio (Q5/Q1), Diferença (Q5−Q1, pp), CIX e SII por indicador, INTERVALO_DE_DATAS",
        },
        subCaptionRelFontSize: null,
        footnote: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
  {
    id: "iceh-inequality-chart",
    label: {
      en: "ICEH inequality chart",
      fr: "Graphique d'inégalité ICEH",
      pt: "Gráfico de desigualdade ICEH",
    },
    description: {
      en: "Chart ranking indicators by a wealth-inequality measure (concentration index by default; change the measure via the value filter). Stratifier set by the replicant selector (wealth quintiles by default).",
      fr: "Graphique classant les indicateurs selon une mesure d'inégalité de richesse (indice de concentration par défaut ; changez la mesure via le filtre de valeurs). Stratificateur défini par le sélecteur de réplicant (quintiles de richesse par défaut).",
      pt: "Gráfico que ordena os indicadores segundo uma medida de desigualdade de riqueza (índice de concentração por predefinição; altere a medida através do filtro de valores). Estratificador definido pelo seletor de replicante (quintis de riqueza por predefinição).",
    },
    allowedFilters: ["iceh_indicator", "strat", "year"],
    createDefaultVisualizationOnInstall: "f2a8d6c3-9b41-4e7a-8d52-3c1f6b9e04a7",
    config: {
      d: {
        type: "chart",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "indicator",
        valuesFilter: ["cix"],
        periodFilter: {
          filterType: "last_calendar_year",
        },
        disaggregateBy: [
          { disOpt: "iceh_indicator", disDisplayOpt: "indicator" },
          { disOpt: "strat", disDisplayOpt: "replicant" },
          { disOpt: "year", disDisplayOpt: "cell" },
        ],
        filterBy: [],
        selectedReplicantValue: "wealth_quintiles",
      },
      s: {
        content: "points",
        horizontal: true,
        colorScale: "single-grey",
        decimalPlaces: 1,
        sortIndicatorValues: "descending",
      },
      t: {
        caption: {
          en: "Wealth inequality by indicator",
          fr: "Inégalité de richesse par indicateur",
          pt: "Desigualdade de riqueza por indicador",
        },
        captionRelFontSize: null,
        subCaption: {
          en: "Concentration index (CIX) across wealth quintiles, DATE_RANGE",
          fr: "Indice de concentration (CIX) entre quintiles de richesse, PLAGE_DE_DATES",
          pt: "Índice de concentração (CIX) entre quintis de riqueza, INTERVALO_DE_DATAS",
        },
        subCaptionRelFontSize: null,
        footnote: null,
        footnoteRelFontSize: null,
      },
    },
    importantNotes: null,
  },
  {
    id: "iceh-inequality-trend",
    label: {
      en: "ICEH inequality trend",
      fr: "Tendance d'inégalité ICEH",
      pt: "Tendência de desigualdade ICEH",
    },
    description: {
      en: "A wealth-inequality measure (concentration index by default; change via the value filter) over survey years, one line per indicator. Filtered to wealth quintiles via the replicant selector.",
      fr: "Une mesure d'inégalité de richesse (indice de concentration par défaut ; changez via le filtre de valeurs) au fil des années d'enquête, une ligne par indicateur. Filtré sur les quintiles de richesse via le sélecteur de réplicant.",
      pt: "Uma medida de desigualdade de riqueza (índice de concentração por predefinição; altere através do filtro de valores) ao longo dos anos de inquérito, uma linha por indicador. Filtrada para quintis de riqueza através do seletor de replicante.",
    },
    allowedFilters: ["iceh_indicator", "strat", "year"],
    createDefaultVisualizationOnInstall: "b6d2c9a4-8e31-4f7a-9b62-1d4f8c2e07a3",
    config: {
      d: {
        type: "timeseries",
        timeseriesGrouping: "year",
        valuesDisDisplayOpt: "series",
        valuesFilter: ["cix"],
        disaggregateBy: [
          { disOpt: "iceh_indicator", disDisplayOpt: "series" },
          { disOpt: "strat", disDisplayOpt: "replicant" },
        ],
        filterBy: [],
        selectedReplicantValue: "wealth_quintiles",
      },
      s: {
        content: "lines",
        decimalPlaces: 1,
      },
      t: {
        caption: {
          en: "Inequality trend",
          fr: "Tendance d'inégalité",
          pt: "Tendência de desigualdade",
        },
        captionRelFontSize: null,
        subCaption: {
          en: "Concentration index (CIX) over time, DATE_RANGE",
          fr: "Indice de concentration (CIX) au fil du temps, PLAGE_DE_DATES",
          pt: "Índice de concentração (CIX) ao longo do tempo, INTERVALO_DE_DATAS",
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
  id: "m9-02-01",
  hide: false,
  resultsObjectId: "M9_iceh_inequality.csv",
  valueProps: ["ratio", "difference", "cix", "sii"],
  valueFunc: "identity",
  valueLabelReplacements: {
    ratio: "Ratio (Q5/Q1)",
    difference: "Difference (Q5−Q1)",
    cix: "CIX",
    sii: "SII",
  },
  label: {
    en: "ICEH inequality measure",
    fr: "Mesure d'inégalité ICEH",
    pt: "Medida de desigualdade ICEH",
  },
  requiredDisaggregationOptions: ["iceh_indicator", "strat", "year"],
  formatAs: "number",
  postAggregationExpression: null,
  aiDescription: {
    summary: {
      en: "Wealth-inequality summary measures derived from ICEH coverage estimates: Ratio (richest/poorest), Difference (richest minus poorest), the concentration index (CIX), and the slope index of inequality (SII). Computed over the ordered wealth quintiles and deciles for every indicator, including the Composite Coverage Index when present. Quantifies how unequally each intervention is distributed across wealth groups.",
      fr: "Mesures synthétiques d'inégalité de richesse dérivées des estimations de couverture ICEH : Ratio (plus riches/plus pauvres), Différence (plus riches moins plus pauvres), indice de concentration (CIX) et indice de pente d'inégalité (SII). Calculées sur les quintiles et déciles de richesse ordonnés pour chaque indicateur, y compris l'Indice de Couverture Composite lorsqu'il est présent. Quantifie l'inégalité de distribution de chaque intervention entre groupes de richesse.",
      pt: "Medidas-resumo de desigualdade de riqueza derivadas das estimativas de cobertura ICEH: Rácio (mais ricos/mais pobres), Diferença (mais ricos menos mais pobres), índice de concentração (CIX) e índice de inclinação da desigualdade (SII). Calculadas sobre os quintis e decis de riqueza ordenados para cada indicador, incluindo o Índice de Cobertura Composto quando presente. Quantifica o grau de desigualdade na distribuição de cada intervenção entre os grupos de riqueza.",
    },
    methodology: {
      en: "Computed in the m009 R script from the ingested stratum-level estimates (not microdata). For each indicator × year × source × wealth stratifier with a complete ordered ladder: Ratio = top/bottom level and Difference = top − bottom level (both exact); CIX uses the Kakwani convenient form with equal wealth-group weights (1/n) and fractional-rank midpoints; SII is the slope of a weighted linear regression of coverage on rank. CIX and SII are GROUPED approximations of the published profile values, which are estimated upstream from individual-level survey microdata (covariance-based CIX, logistic SII) that FASTR does not have. Values are produced only when the full quintile/decile set is present.",
      fr: "Calculées dans le script R m009 à partir des estimations agrégées par strate (pas de microdonnées). Pour chaque indicateur × année × source × stratificateur de richesse avec une échelle complète : Ratio = niveau supérieur/inférieur et Différence = niveau supérieur − inférieur (tous deux exacts) ; le CIX utilise la forme commode de Kakwani avec des poids égaux par groupe de richesse (1/n) et des rangs fractionnaires médians ; le SII est la pente d'une régression linéaire pondérée de la couverture sur le rang. Le CIX et le SII sont des approximations GROUPÉES des valeurs publiées, estimées en amont à partir des microdonnées d'enquête (CIX par covariance, SII logistique) que FASTR ne possède pas.",
      pt: "Calculadas no script R do m009 a partir das estimativas ao nível dos estratos ingeridas (não microdados). Para cada indicador × ano × fonte × estratificador de riqueza com uma escala ordenada completa: Rácio = nível superior/inferior e Diferença = nível superior − inferior (ambos exatos); o CIX utiliza a forma conveniente de Kakwani com pesos iguais por grupo de riqueza (1/n) e pontos médios de posição fracionária; o SII é o declive de uma regressão linear ponderada da cobertura sobre a posição. O CIX e o SII são aproximações AGRUPADAS dos valores publicados nos perfis, estimados a montante a partir dos microdados de inquérito ao nível individual (CIX por covariância, SII logístico) que o FASTR não possui. Os valores só são produzidos quando o conjunto completo de quintis/decis está presente.",
    },
    interpretation: {
      en: "Ratio > 1 and positive CIX/SII mean coverage is concentrated among the richest (pro-rich inequality); Ratio < 1 and negative values mean it is concentrated among the poorest (e.g. zero-dose children, stunting). Difference and SII are in percentage points; an absolute CIX above 30 indicates high inequality. A value near 0 (Ratio near 1) means coverage is roughly equal across wealth groups.",
      fr: "Un Ratio > 1 et des CIX/SII positifs signifient que la couverture est concentrée chez les plus riches (inégalité pro-riches) ; un Ratio < 1 et des valeurs négatives signifient une concentration chez les plus pauvres (ex : enfants zéro-dose, retard de croissance). La Différence et le SII sont en points de pourcentage ; un CIX absolu supérieur à 30 indique une forte inégalité. Une valeur proche de 0 (Ratio proche de 1) signifie une couverture à peu près égale entre groupes.",
      pt: "Um Rácio > 1 e valores positivos de CIX/SII significam que a cobertura está concentrada entre os mais ricos (desigualdade pró-ricos); um Rácio < 1 e valores negativos significam que está concentrada entre os mais pobres (por exemplo, crianças com dose zero, atraso de crescimento). A Diferença e o SII estão em pontos percentuais; um CIX absoluto superior a 30 indica forte desigualdade. Um valor próximo de 0 (Rácio próximo de 1) significa que a cobertura é aproximadamente igual entre os grupos de riqueza.",
    },
    typicalRange: {
      en: "Ratio is unitless (commonly 0.1–10). Difference and SII are in percentage points (roughly −100 to 100). CIX ranges about −100 to 100, with absolute values above 30 considered high inequality.",
      fr: "Le Ratio est sans unité (généralement 0,1–10). La Différence et le SII sont en points de pourcentage (environ −100 à 100). Le CIX varie d'environ −100 à 100, les valeurs absolues supérieures à 30 étant considérées comme une forte inégalité.",
      pt: "O Rácio é adimensional (geralmente 0,1–10). A Diferença e o SII estão em pontos percentuais (cerca de −100 a 100). O CIX varia aproximadamente entre −100 e 100, sendo os valores absolutos superiores a 30 considerados como forte desigualdade.",
    },
    caveats: {
      en: "CIX and SII are grouped approximations and will NOT exactly reproduce the published ICEH profile figures (which use microdata and a logistic SII); only Ratio and Difference are exact. Measures are emitted only for complete wealth ladders (all 5 quintiles or all 10 deciles) — partial sets are silently omitted. Inequality is computed for wealth stratifiers only. Standard errors are not propagated (stored as NA). Each cell reflects a single survey source; if one indicator-year has both a DHS and a MICS, filter or disaggregate accordingly.",
      fr: "Le CIX et le SII sont des approximations groupées et ne reproduiront PAS exactement les chiffres publiés des profils ICEH (qui utilisent des microdonnées et un SII logistique) ; seuls le Ratio et la Différence sont exacts. Les mesures ne sont produites que pour des échelles de richesse complètes (les 5 quintiles ou les 10 déciles) — les ensembles partiels sont omis. L'inégalité est calculée uniquement pour les stratificateurs de richesse. Les erreurs-types ne sont pas propagées (stockées en NA). Chaque cellule reflète une seule source d'enquête.",
      pt: "O CIX e o SII são aproximações agrupadas e NÃO reproduzem exatamente os valores publicados nos perfis ICEH (que utilizam microdados e um SII logístico); apenas o Rácio e a Diferença são exatos. As medidas só são geradas para escalas de riqueza completas (todos os 5 quintis ou todos os 10 decis) — os conjuntos parciais são omitidos silenciosamente. A desigualdade é calculada apenas para estratificadores de riqueza. Os erros-padrão não são propagados (armazenados como NA). Cada célula reflete uma única fonte de inquérito; se um indicador-ano tiver tanto um DHS como um MICS, filtre ou desagregue em conformidade.",
    },
    disaggregationGuidance: {
      en: "Data structure: one row per iceh_indicator × year × source × strat, with the four measures as value columns (ratio, difference, cix, sii). For the profile's equity table, place iceh_indicator as rows and the four measures as columns (valuesDisDisplayOpt: col), filter strat to wealth_quintiles (or wealth_deciles), and select a year. There is no level dimension — the wealth levels are collapsed into the measures.",
      fr: "Structure des données : une ligne par iceh_indicator × année × source × strat, avec les quatre mesures en colonnes de valeur (ratio, difference, cix, sii). Pour le tableau d'équité du profil, placez iceh_indicator en lignes et les quatre mesures en colonnes (valuesDisDisplayOpt : col), filtrez strat sur wealth_quintiles (ou wealth_deciles) et sélectionnez une année. Il n'y a pas de dimension level — les niveaux de richesse sont condensés dans les mesures.",
      pt: "Estrutura dos dados: uma linha por iceh_indicator × ano × fonte × strat, com as quatro medidas em colunas de valor (ratio, difference, cix, sii). Para a tabela de equidade do perfil, coloque iceh_indicator em linhas e as quatro medidas em colunas (valuesDisDisplayOpt: col), filtre strat para wealth_quintiles (ou wealth_deciles) e selecione um ano. Não existe dimensão level — os níveis de riqueza estão condensados nas medidas.",
    },
  },
  variantLabel: null,
  importantNotes: null,
  vizPresets,
};
