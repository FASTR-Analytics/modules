import { walk } from "jsr:@std/fs@^1/walk";
import { moduleDefinitionGithubSchema } from "./_module_definition_github.ts";

const ROOT = new URL("..", import.meta.url).pathname;

let hasErrors = false;

for await (const entry of walk(ROOT, {
  match: [/\/definition\.json$/],
  skip: [/\/\./, /\/node_modules\//],
})) {
  const relativePath = entry.path.replace(ROOT, "");
  let raw: unknown;

  try {
    raw = JSON.parse(await Deno.readTextFile(entry.path));
  } catch (e) {
    console.error(`FAIL  ${relativePath}`);
    console.error(`      JSON parse error: ${e}`);
    hasErrors = true;
    continue;
  }

  const result = moduleDefinitionGithubSchema.safeParse(raw);
  if (result.success) {
    console.log(`OK    ${relativePath}`);
  } else {
    console.error(`FAIL  ${relativePath}`);
    for (const issue of result.error.issues) {
      const path = issue.path.length > 0 ? issue.path.join(".") + ": " : "";
      console.error(`      ${path}${issue.message}`);
    }
    hasErrors = true;
  }
}

if (hasErrors) {
  console.error("\nValidation failed.");
  Deno.exit(1);
} else {
  console.log("\nAll definition.json files are valid.");
}
