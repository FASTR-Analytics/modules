import type { ModuleParameter } from "../.validation/_module_definition_github.ts";

export const parameters: ModuleParameter[] = [
  {
    description: {
      en: "Count value to use",
      fr: "Valeur de comptage à utiliser",
    },
    replacementString: "SELECTED_COUNT_VARIABLE",
    input: {
      inputType: "select",
      options: [
        {
          value: "count_final_none",
          label: "Count (unadjusted)",
        },
        {
          value: "count_final_outliers",
          label: "Count (adjusted for outliers)",
        },
        {
          value: "count_final_completeness",
          label: "Count (adjusted for completeness)",
        },
        {
          value: "count_final_both",
          label: "Count (adjusted for outliers and completeness)",
        },
      ],
      valueType: "string",
      defaultValue: "count_final_outliers",
    },
  },
  {
    description: {
      en: "Level to calculate coverage for",
      fr: "Niveau pour calculer la couverture",
    },
    replacementString: "ANALYSIS_LEVEL",
    input: {
      inputType: "select",
      options: [
        {
          value: "NATIONAL_ONLY",
          label: "National only",
        },
        {
          value: "NATIONAL_PLUS_AA2",
          label: "National and admin area 2",
        },
        {
          value: "NATIONAL_PLUS_AA2_AA3",
          label: "National, admin area 2, and admin area 3",
        },
      ],
      valueType: "string",
      defaultValue: "NATIONAL_PLUS_AA2",
    },
  },
  {
    description: {
      en: "Pregnancy loss rate",
      fr: "Taux de perte de grossesse",
    },
    replacementString: "PREGNANCY_LOSS_RATE",
    input: {
      inputType: "number",
      defaultValue: "0.03",
    },
  },
  {
    description: {
      en: "Twin rate",
      fr: "Taux de jumeaux",
    },
    replacementString: "TWIN_RATE",
    input: {
      inputType: "number",
      defaultValue: "0.015",
    },
  },
  {
    description: {
      en: "Stillbirth rate",
      fr: "Taux de mortinatalité",
    },
    replacementString: "STILLBIRTH_RATE",
    input: {
      inputType: "number",
      defaultValue: "0.02",
    },
  },
  {
    description: {
      en: "Neonatal mortality rate",
      fr: "Taux de mortalité néonatale",
    },
    replacementString: "P1_NMR",
    input: {
      inputType: "number",
      defaultValue: "0.039",
    },
  },
  {
    description: {
      en: "Postneonatal mortality rate",
      fr: "Taux de mortalité postnéonatale",
    },
    replacementString: "P2_PNMR",
    input: {
      inputType: "number",
      defaultValue: "0.028",
    },
  },
  {
    description: {
      en: "Infant mortality rate",
      fr: "Taux de mortalité infantile",
    },
    replacementString: "INFANT_MORTALITY_RATE",
    input: {
      inputType: "number",
      defaultValue: "0.067",
    },
  },
  {
    description: {
      en: "Under 5 mortality rate",
      fr: "Taux de mortalité des moins de 5 ans",
    },
    replacementString: "UNDER5_MORTALITY_RATE",
    input: {
      inputType: "number",
      defaultValue: "0.103",
    },
  },
];
