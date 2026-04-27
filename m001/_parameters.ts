import type { ModuleParameter } from "../.validation/_module_definition_github.ts";

export const parameters: ModuleParameter[] = [
  {
    description: {
      en: "Proportion threshold for outlier detection",
      fr: "Seuil de proportion pour la détection des valeurs aberrantes",
    },
    replacementString: "OUTLIER_PROPORTION_THRESHOLD",
    input: {
      inputType: "number",
      defaultValue: "0.8",
    },
  },
  {
    description: {
      en: "Minimum count threshold for consideration",
      fr: "Seuil de comptage minimum pour considération",
    },
    replacementString: "MINIMUM_COUNT_THRESHOLD",
    input: {
      inputType: "number",
      defaultValue: "100",
    },
  },
  {
    description: {
      en: "Number of MADs",
      fr: "Nombre de MADs",
    },
    replacementString: "MADS",
    input: {
      inputType: "number",
      defaultValue: "10",
    },
  },
  {
    description: {
      en: "Indicators subjected to DQA",
      fr: "Indicateurs soumis à l'AQD",
    },
    replacementString: "DQA_INDICATORS",
    input: {
      inputType: "select",
      options: [
        {
          value: 'c("anc1", "penta1", "opd")',
          label: "anc1, penta1, opd",
        },
        {
          value: 'c("anc1", "penta1")',
          label: "anc1, penta1",
        },
        {
          value: 'c("anc1", "opd")',
          label: "anc1, opd",
        },
        {
          value: 'c("penta1", "opd")',
          label: "penta1, opd",
        },
        {
          value: 'c("anc1")',
          label: "anc1",
        },
        {
          value: 'c("penta1")',
          label: "penta1",
        },
        {
          value: 'c("opd")',
          label: "opd",
        },
      ],
      valueType: "number",
      defaultValue: 'c("anc1", "penta1", "opd")',
    },
  },
  {
    description: {
      en: "Consistency pairs used",
      fr: "Paires de cohérence utilisées",
    },
    replacementString: "CONSISTENCY_PAIRS_USED",
    input: {
      inputType: "select",
      options: [
        {
          value: 'c("penta", "anc", "delivery", "malaria")',
          label: "penta, anc, delivery, malaria",
        },
        {
          value: 'c("penta", "anc", "delivery")',
          label: "penta, anc, delivery",
        },
        {
          value: 'c("penta", "anc")',
          label: "penta, anc",
        },
        {
          value: 'c("penta", "delivery")',
          label: "penta, delivery",
        },
        {
          value: 'c("anc", "delivery")',
          label: "anc, delivery",
        },
        {
          value: 'c("anc")',
          label: "anc",
        },
        {
          value: 'c("delivery")',
          label: "delivery",
        },
        {
          value: 'c("penta")',
          label: "penta",
        },
        {
          value: 'c("malaria")',
          label: "malaria",
        },
      ],
      valueType: "number",
      defaultValue: 'c("penta", "anc", "delivery")',
    },
  },
  {
    description: {
      en: "Admin level used to join facilities to corresponding geo-consistency",
      fr: "Niveau administratif utilisé pour joindre les établissements à la cohérence géographique correspondante",
    },
    replacementString: "GEOLEVEL",
    input: {
      inputType: "select",
      options: [
        {
          value: "admin_area_2",
          label: "admin_area_2",
        },
        {
          value: "admin_area_3",
          label: "admin_area_3",
        },
        {
          value: "admin_area_4",
          label: "admin_area_4",
        },
      ],
      valueType: "string",
      defaultValue: "admin_area_3",
    },
  },
];
