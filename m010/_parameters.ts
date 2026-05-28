import type { ModuleParameter } from "../.validation/_module_definition_github.ts";

export const parameters: ModuleParameter[] = [
  {
    replacementString: "USE_SAMPLE_WEIGHTS",
    description: {
      en: "Use sample weights for calculations",
      fr: "Utiliser les pondérations de l'échantillon pour les calculs",
    },
    input: {
      inputType: "boolean",
      defaultValue: "FALSE",
    },
  },
  {
    replacementString: "STOP_IF_INDICATOR_FAILS",
    description: {
      en: "Stop script if any single indicator fails (no R code or errors). If FALSE, skip failing indicators and only fail if no results are generated.",
      fr: "Arrêter le script si un indicateur échoue (pas de code R ou erreurs). Si FALSE, ignorer les indicateurs qui échouent et échouer uniquement si aucun résultat n'est généré.",
    },
    input: {
      inputType: "boolean",
      defaultValue: "TRUE",
    },
  },
];
