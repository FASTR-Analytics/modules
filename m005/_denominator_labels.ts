// The labels of every value the denominator column can hold, built the way
// script.R builds `denominator_label` (add_denominator_labels): one id per
// source chain and target population, `d<source>_<target>`, plus the
// `source_*` and `target_*` ids the coverage-by-denominator rows carry and
// the best/survey choice. Shared by the six m005 metrics so a value the
// script emits is never shown as its raw id.

const SOURCES: Record<string, { label: string; phrase: string }> = {
  anc1: {
    label: "HMIS data on ANC 1st visits",
    phrase: "derived from HMIS data on ANC 1st visits",
  },
  delivery: {
    label: "HMIS data on institutional deliveries",
    phrase: "derived from HMIS data on institutional deliveries",
  },
  sba: {
    label: "HMIS data on skilled birth attendance",
    phrase: "derived from HMIS data on skilled birth attendance",
  },
  bcg: {
    label: "HMIS data on BCG doses",
    phrase: "derived from HMIS data on BCG doses",
  },
  penta1: {
    label: "HMIS data on Penta-1 doses",
    phrase: "derived from HMIS data on Penta-1 doses",
  },
  livebirths: {
    label: "HMIS data on live births",
    phrase: "derived from HMIS data on live births",
  },
  wpp: { label: "UN WPP estimates", phrase: "based on UN WPP estimates" },
};

const TARGETS: Record<string, { label: string; phrase: string }> = {
  pregnancy: { label: "Pregnancies", phrase: "pregnancies" },
  delivery: { label: "Deliveries", phrase: "deliveries" },
  birth: {
    label: "Total births (live + stillbirths)",
    phrase: "total births (live + stillbirths)",
  },
  livebirth: { label: "Live births", phrase: "live births" },
  dpt: {
    label: "Infants eligible for DPT1",
    phrase: "infants eligible for DPT1",
  },
  measles1: {
    label: "Children eligible for measles dose 1 (MCV1)",
    phrase: "children eligible for measles dose 1 (MCV1)",
  },
  measles2: {
    label: "Children eligible for measles dose 2 (MCV2)",
    phrase: "children eligible for measles dose 2 (MCV2)",
  },
  vitaminA: {
    label: "Children aged 6-59 months eligible for Vitamin A supplementation",
    phrase: "children aged 6-59 months eligible for Vitamin A supplementation",
  },
  fully_immunized: {
    label: "Children under 1 year eligible for full immunization",
    phrase: "children under 1 year eligible for full immunization",
  },
};

export const DENOMINATOR_VALUE_LABELS: Record<string, string> = {
  best: "Best",
  survey: "Survey",
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
        `Estimated number of ${t.phrase} ${s.phrase}.`,
      ])
    ),
  ),
};
