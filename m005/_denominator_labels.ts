// The labels of every value the denominator column can hold, built the way
// script.R builds `denominator_label` (add_denominator_labels): one id per
// source chain and target population, `d<source>_<target>`, plus the
// `source_*` and `target_*` ids the coverage-by-denominator rows carry and
// the best/survey choice. Shared by the six m005 metrics so a value the
// script emits is never shown as its raw id.

type Localized = { en: string; fr: string; pt: string };

// `phrase` carries its own preposition so each language composes without
// elision trouble ("d'accouchements", not "de accouchements").
const SOURCES: Record<string, { label: Localized; phrase: Localized }> = {
  anc1: {
    label: {
      en: "HMIS data on ANC 1st visits",
      fr: "Données HMIS sur les premières consultations prénatales",
      pt: "Dados HMIS sobre as primeiras consultas pré-natais",
    },
    phrase: {
      en: "derived from HMIS data on ANC 1st visits",
      fr: "dérivé des données HMIS sur les premières consultations prénatales",
      pt: "derivado dos dados HMIS sobre as primeiras consultas pré-natais",
    },
  },
  delivery: {
    label: {
      en: "HMIS data on institutional deliveries",
      fr: "Données HMIS sur les accouchements en établissement",
      pt: "Dados HMIS sobre os partos institucionais",
    },
    phrase: {
      en: "derived from HMIS data on institutional deliveries",
      fr: "dérivé des données HMIS sur les accouchements en établissement",
      pt: "derivado dos dados HMIS sobre os partos institucionais",
    },
  },
  sba: {
    label: {
      en: "HMIS data on skilled birth attendance",
      fr: "Données HMIS sur l'assistance qualifiée à l'accouchement",
      pt: "Dados HMIS sobre a assistência qualificada ao parto",
    },
    phrase: {
      en: "derived from HMIS data on skilled birth attendance",
      fr: "dérivé des données HMIS sur l'assistance qualifiée à l'accouchement",
      pt: "derivado dos dados HMIS sobre a assistência qualificada ao parto",
    },
  },
  bcg: {
    label: {
      en: "HMIS data on BCG doses",
      fr: "Données HMIS sur les doses de BCG",
      pt: "Dados HMIS sobre as doses de BCG",
    },
    phrase: {
      en: "derived from HMIS data on BCG doses",
      fr: "dérivé des données HMIS sur les doses de BCG",
      pt: "derivado dos dados HMIS sobre as doses de BCG",
    },
  },
  penta1: {
    label: {
      en: "HMIS data on Penta-1 doses",
      fr: "Données HMIS sur les doses de Penta-1",
      pt: "Dados HMIS sobre as doses de Penta-1",
    },
    phrase: {
      en: "derived from HMIS data on Penta-1 doses",
      fr: "dérivé des données HMIS sur les doses de Penta-1",
      pt: "derivado dos dados HMIS sobre as doses de Penta-1",
    },
  },
  livebirths: {
    label: {
      en: "HMIS data on live births",
      fr: "Données HMIS sur les naissances vivantes",
      pt: "Dados HMIS sobre os nascimentos vivos",
    },
    phrase: {
      en: "derived from HMIS data on live births",
      fr: "dérivé des données HMIS sur les naissances vivantes",
      pt: "derivado dos dados HMIS sobre os nascimentos vivos",
    },
  },
  wpp: {
    label: {
      en: "UN WPP estimates",
      fr: "Estimations WPP de l'ONU",
      pt: "Estimativas WPP da ONU",
    },
    phrase: {
      en: "based on UN WPP estimates",
      fr: "basé sur les estimations WPP de l'ONU",
      pt: "baseado nas estimativas WPP da ONU",
    },
  },
};

const TARGETS: Record<string, { label: Localized; phrase: Localized }> = {
  pregnancy: {
    label: { en: "Pregnancies", fr: "Grossesses", pt: "Gravidezes" },
    phrase: { en: "of pregnancies", fr: "de grossesses", pt: "de gravidezes" },
  },
  delivery: {
    label: { en: "Deliveries", fr: "Accouchements", pt: "Partos" },
    phrase: { en: "of deliveries", fr: "d'accouchements", pt: "de partos" },
  },
  birth: {
    label: {
      en: "Total births (live + stillbirths)",
      fr: "Naissances totales (vivantes + mortinaissances)",
      pt: "Nascimentos totais (vivos + fetais)",
    },
    phrase: {
      en: "of total births (live + stillbirths)",
      fr: "de naissances totales (vivantes + mortinaissances)",
      pt: "de nascimentos totais (vivos + fetais)",
    },
  },
  livebirth: {
    label: {
      en: "Live births",
      fr: "Naissances vivantes",
      pt: "Nascimentos vivos",
    },
    phrase: {
      en: "of live births",
      fr: "de naissances vivantes",
      pt: "de nascimentos vivos",
    },
  },
  dpt: {
    label: {
      en: "Infants eligible for DPT1",
      fr: "Nourrissons éligibles au DPT1",
      pt: "Lactentes elegíveis para DPT1",
    },
    phrase: {
      en: "of infants eligible for DPT1",
      fr: "de nourrissons éligibles au DPT1",
      pt: "de lactentes elegíveis para DPT1",
    },
  },
  measles1: {
    label: {
      en: "Children eligible for measles dose 1 (MCV1)",
      fr: "Enfants éligibles à la 1re dose contre la rougeole (MCV1)",
      pt: "Crianças elegíveis para a 1.ª dose contra o sarampo (MCV1)",
    },
    phrase: {
      en: "of children eligible for measles dose 1 (MCV1)",
      fr: "d'enfants éligibles à la 1re dose contre la rougeole (MCV1)",
      pt: "de crianças elegíveis para a 1.ª dose contra o sarampo (MCV1)",
    },
  },
  measles2: {
    label: {
      en: "Children eligible for measles dose 2 (MCV2)",
      fr: "Enfants éligibles à la 2e dose contre la rougeole (MCV2)",
      pt: "Crianças elegíveis para a 2.ª dose contra o sarampo (MCV2)",
    },
    phrase: {
      en: "of children eligible for measles dose 2 (MCV2)",
      fr: "d'enfants éligibles à la 2e dose contre la rougeole (MCV2)",
      pt: "de crianças elegíveis para a 2.ª dose contra o sarampo (MCV2)",
    },
  },
  vitaminA: {
    label: {
      en: "Children aged 6-59 months eligible for Vitamin A supplementation",
      fr: "Enfants de 6 à 59 mois éligibles à la supplémentation en vitamine A",
      pt: "Crianças de 6 a 59 meses elegíveis para suplementação de vitamina A",
    },
    phrase: {
      en: "of children aged 6-59 months eligible for Vitamin A supplementation",
      fr:
        "d'enfants de 6 à 59 mois éligibles à la supplémentation en vitamine A",
      pt:
        "de crianças de 6 a 59 meses elegíveis para suplementação de vitamina A",
    },
  },
  fully_immunized: {
    label: {
      en: "Children under 1 year eligible for full immunization",
      fr: "Enfants de moins d'un an éligibles à la vaccination complète",
      pt: "Crianças com menos de 1 ano elegíveis para vacinação completa",
    },
    phrase: {
      en: "of children under 1 year eligible for full immunization",
      fr: "d'enfants de moins d'un an éligibles à la vaccination complète",
      pt: "de crianças com menos de 1 ano elegíveis para vacinação completa",
    },
  },
};

const composed = (target: Localized, source: Localized): Localized => ({
  en: `Estimated number ${target.en} ${source.en}.`,
  fr: `Nombre estimé ${target.fr} ${source.fr}.`,
  pt: `Número estimado ${target.pt} ${source.pt}.`,
});

export const DENOMINATOR_VALUE_LABELS: Record<string, Localized> = {
  best: { en: "Best", fr: "Meilleur", pt: "Melhor" },
  survey: { en: "Survey", fr: "Enquête", pt: "Inquérito" },
  ...Object.fromEntries(
    Object.entries(SOURCES).map(([id, s]) => [`source_${id}`, s.label]),
  ),
  ...Object.fromEntries(
    Object.entries(TARGETS).map(([id, t]) => [`target_${id}`, t.label]),
  ),
  ...Object.fromEntries(
    Object.entries(SOURCES).flatMap(([sourceId, s]) =>
      Object.entries(TARGETS).map(([targetId, t]) => [
        `d${sourceId}_${targetId}`,
        composed(t.phrase, s.phrase),
      ])
    ),
  ),
};
