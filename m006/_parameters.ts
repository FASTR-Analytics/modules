import type { ModuleParameter } from "../.validation/_module_definition_github.ts";

export const parameters: ModuleParameter[] = [
  {
    description: {
      en: "Denominator to use for all coverage calculations",
      fr: "Dénominateur à utiliser pour tous les calculs de couverture",
    },
    replacementString: "DENOMINATOR_CHAIN",
    input: {
      inputType: "select",
      options: [
        {
          value: "auto",
          label: "auto",
        },
        {
          value: "anc1",
          label: "anc1",
        },
        {
          value: "delivery",
          label: "delivery",
        },
        {
          value: "bcg",
          label: "bcg",
        },
        {
          value: "penta1",
          label: "penta1",
        },
      ],
      valueType: "string",
      defaultValue: "auto",
    },
  },
];
