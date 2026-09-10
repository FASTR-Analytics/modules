# Briefing: reviewing module text for AI consumers

Read this before doing any lane work. It states what the review is for, how
the text you are reviewing reaches an AI agent, what "accurate" means here,
and what is already known. It was written from a code read of both repos on
2026-09-07. The orchestrator checks before each run whether the wb-fastr
files named in section 9 still behave as described, and passes a
`briefing_delta.md` to every lane if anything changed. Lanes do not need to
re-check.

Acronyms used here: HMIS (health management information system, the
routine facility reporting data), HFA (health facility assessment survey),
ICEH (International Center for Equity in Health, whose survey-based equity
estimates m009 analyses), MAD (median absolute deviation), DQA (data
quality assessment), ANC1 (first antenatal care visit), NMR (neonatal
mortality rate), SPA (the browser app), MCP (the protocol by which an
external Claude client connects to an instance).

## 1. Purpose and scope

The module definitions in this repo carry human-authored text: metric labels,
`aiDescription` blocks, `importantNotes`, value labels, preset labels and
notes, parameter descriptions, captions and footnotes. Some of that text is
delivered verbatim to AI agents (the in-app FASTR assistant and any Claude
client connected over `/mcp`) as the only description they get of what a
metric is and how to read its numbers. The agent cannot see the R code, the
results object, the aggregation function, or the expression. An error in the
text becomes an error in what the agent tells users.

The review asks, for every piece of text an AI receives: is it accurate
against what the code computes, and is it clear enough that an agent with no
other context would choose the right metric, query it correctly, and
interpret the numbers correctly?

Two tiers:

- **Tier 1 (primary)**: text that reaches an AI. See section 3.
- **Tier 2 (secondary)**: text the AI cannot see but that ends up on figures
  the AI creates for people: preset captions, sub-captions, footnotes, preset
  descriptions. Review these for accuracy only where the lane has already
  established the ground truth for the metric; do not spend effort on style.

Out of scope: R code correctness, the app's own formatters and prompts
(record observations about them separately, do not propose app changes as
findings), and the frozen modules m007 and m008.

## 2. The pipeline in one page

1. Authors edit `mNNN/_core.ts`, `_parameters.ts`, `_results_objects.ts`,
   `_metrics/<id>.ts`, `script.R`. `deno task build` compiles them into
   `mNNN/definition.json`, which is what the app reads. Review the `.ts`
   sources (they are the edit target) but trust `definition.json` as what
   ships.
2. wb-fastr fetches `definition.json` and `script.R` from GitHub at a pinned
   commit (`server/module_loader/load_module.ts`). It validates with a
   strict zod schema and translates three metric fields into the instance
   language at that moment: `label`, `variantLabel`, `importantNotes`, plus
   the module `label` and every parameter `description`. Everything else
   (`aiDescription`, all vizPreset text) stays as `{en, fr, pt}`.
3. A run executes the R scripts and writes an immutable results package. The
   metric rows in its manifest carry the translated strings frozen at
   generation time.
4. Two AI surfaces read the package:
   - the SPA copilot inside the app (42 tools), and
   - the `/mcp` endpoint (6 read-only tools: `get_overview`,
     `get_available_metrics`, `get_metric_data`, two methodology-doc tools,
     `get_info`).
   Both use the same two shared metric formatters. Only the SPA has the module
   tools (`get_available_modules`, `get_module_settings`, script, logs).
5. The AI formatters take `.en` from every `{en, fr, pt}` field that is
   still untranslated at that point. So on a French or Portuguese instance
   the agent reads a French label followed by English methodology. (The
   loader's translation step in item 2 does fall back to `en` when `pt` is
   absent; the formatters never look at `fr` or `pt` at all.)

## 3. Which fields reach an AI, where, and in what language

The two formatters are `lib/ai_tools/format_metrics_list_for_ai.ts`
(`get_available_metrics`, the catalog every conversation starts from) and
`lib/ai_tools/format_metric_data_for_ai.ts` (`get_metric_data`, the header
that precedes the CSV). `render_ai_view.ts` in this directory reproduces
both from `definition.json`; run it and read its output rather than
reconstructing this by hand.

| Field | Catalog (`get_available_metrics`) | Data header (`get_metric_data`) | Language the AI sees |
| --- | --- | --- | --- |
| metric `id` | yes | yes | n/a |
| metric `label` | yes, as `Label [Variant]` | yes | instance language |
| `variantLabel` | yes, in brackets | yes | instance language |
| `formatAs` | yes, in brackets | yes, with an explainer for `indicator` | n/a |
| `aiDescription.summary` | yes, one indented line | **no** | en only |
| `aiDescription.methodology` | no | yes | en only |
| `aiDescription.interpretation` | no | yes | en only |
| `aiDescription.typicalRange` | no | yes | en only |
| `aiDescription.caveats` | no | yes (if not null) | en only |
| `aiDescription.disaggregationGuidance` | no | yes | en only |
| `importantNotes` (metric) | yes, `NOTE:` | yes, `**IMPORTANT:**` | instance language |
| `valueProps` | yes, `Values:` line | yes, `Value properties:` list | n/a |
| `valueLabelReplacements` | yes, as `prop (label)` | yes, as `prop: label` | as authored (single string) |
| `requiredDisaggregationOptions` | yes, `Auto-disaggregated by:` | no explicit line; the required columns are present in the CSV | n/a |
| available optional disaggregations | yes, derived from the results object columns | yes | n/a |
| vizPreset `id` and `label` | yes, one line per preset | no | **`.en` hard-coded** |
| vizPreset `importantNotes` | yes, `NOTE:` under the preset | no | en only |
| vizPreset `allowedFilters` | yes, `filters:` | no | n/a |
| vizPreset replicant requirement (a replicant is the dimension a preset repeats one small figure per value of, such as one map per indicator; the agent must pick one value) | yes, derived from config | no | n/a |
| vizPreset `description` | **never** | **never** | n/a |
| vizPreset `config.t` caption, subCaption, footnote | **never** | **never** | n/a |
| `valueFunc`, `resultsObjectId`, `postAggregationExpression`, `catalogExpressionEvaluation` | never | never | n/a |
| `hide: true` metrics | filtered out entirely | still queryable by id | n/a |
| module `label` | SPA only, `get_available_modules` | n/a | instance language |
| parameter `description` | SPA only, `get_module_settings`, as `description: value` | n/a | instance language |

Consequences that shape the review:

- The catalog is the agent's only basis for choosing among 46 visible
  metrics. It sees `id`, `Label [Variant] [formatAs]`, `summary`, the values
  line, the disaggregation lines, and the preset lines. Nothing else.
- `summary` and the other five `aiDescription` fields never appear together.
  An agent that has only listed metrics has never seen a caveat.
- The agent picks a preset from `id` and `label.en` alone. `description` is
  invisible to it, and so is the caption the preset will render.
- `get_module_settings` prints parameters as `description: value`. The
  `replacementString` (for example `DIFFPERCENT`) is never shown. Text that
  names a parameter must use the description wording, and on a French
  instance the description is French while `aiDescription` is English.
- On `/mcp` there are no module tools at all. Text that says "check the
  module settings" is unactionable there; state the default value instead.

## 4. What a displayed number actually is

The agent needs to know what one value means at any grouping. Establish it
from the definition, not from the prose:

- `valueFunc` is applied by SQL to each `valueProp` over whatever the query
  groups by. `AVG` of a 0/1 flag is a proportion of rows (usually
  facility-months). `SUM` of counts is a total. `identity` means no
  aggregation and only makes sense with `postAggregationExpression` or
  `catalogExpressionEvaluation`, or when the results object is already at
  the displayed grain.
- `postAggregationExpression` aggregates each ingredient with its own
  function, then applies the arithmetic once. A ratio of sums is not an
  average of ratios; the text should say which it is.
- `catalogExpressionEvaluation` (m012 only) sums the ingredient columns and
  applies each indicator's own formula from the instance's indicator
  dictionary.
- `formatAs: "percent"` values travel to the AI as 0-1 fractions with three
  decimals. `"number"` is two decimals. `"indicator"` defers to each
  indicator's own format, listed per indicator in the "Dimension Summary"
  section that get_metric_data prints below its header. The app force-sets `"indicator"` for the ids in wb-fastr
  `lib/indicator_format_metrics.ts`; the authored value is overridden.
- The results object grain sets the meaning of a row. Read the R code that
  writes the CSV to know whether a row is a facility-month, an area-month, an
  area-year, or a national-year, and whether a value in it is already an
  aggregate (in which case re-summing across areas may double count).

Which disaggregations a metric offers is derived from the results object's
declared columns (wb-fastr `server/runs/disaggregation_availability.ts`):

- physical columns present from this list: admin_area_2/3/4,
  indicator_common_id, denominator, denominator_best_or_survey,
  source_indicator, target_population, ratio_type, hfa_*, time_point,
  iceh_indicator, strat, level;
- `facility_type`, `facility_ownership`, `facility_custom_1..5` only if the
  results object has `facility_id` and the instance enables that column;
- `year`, `month`, `quarter_id`, `period_id` if `period_id` is present;
  `quarter_id`, `year` if only `quarter_id`; `year` if only `year`.
- `quarter_id` is hidden from the AI's lists even when available.
- A results object with `createTableStatementPossibleColumns: false` has no
  declared columns; treat any disaggregation guidance for its metrics as
  unverifiable and say so.

`disaggregationGuidance` must name only dimensions that exist for that
metric, using the `disOpt` ids the agent passes (`admin_area_2`, not
"region"). Required options are always included in a query; guidance that
says "always disaggregate by X" where X is required is redundant but
harmless, and where X is not required it is a real instruction.

## 5. Quality criteria, per field

Apply these in order of consequence. A finding must cite the ground truth it
was checked against.

**Catalog line (`label`, `variantLabel`, `formatAs`, `summary`).**
The line must let an agent choose correctly among siblings without opening
the metric. For siblings (same label, different variant), the variant text
on its own must be enough to tell them apart. The summary is one sentence stating what one value is:
the quantity, the unit or format, the grain (per facility-month, per
area-year), and the direction if it is not obvious. It should not repeat the
label. It should not describe the module.

**Methodology.** State the computation as it happens: what the R script
produces per row, and what the app then does (`AVG` of a flag, `SUM` of
counts, ratio of sums, expression). Name parameter defaults where they
change the result, using the number and the parameter's description wording.
Do not describe steps the agent cannot see or that do not affect the
displayed value. Every number in the text must match `_parameters.ts` or a
constant in `script.R`, and the field should say when a parameter can move
it.

**Interpretation.** Direction (is higher good or bad), what pushes values up
or down, the common misreadings for this metric. Must agree with the sign
convention of the expression (for example a negative `pct_diff` is a
shortfall). If the metric can exceed 100% or go negative, say so and why.

**Typical range.** Honest and sourced. A threshold is either a parameter
default, a constant in the script, an established convention (say whose), or
expert judgement (say so). Values reach the agent as 0-1 fractions, and the
tool tells the agent so. A range may therefore be written in percent,
provided it says "percent" or uses the % sign and does not mix the two
forms.

**Caveats.** Real limitations of this metric: sensitivity to a parameter,
data-quality dependence, denominator assumptions, grains at which the
number stops being meaningful. Not a restatement of interpretation. If the
caveat depends on a setting, state the default so an `/mcp` agent can act.

**Disaggregation guidance.** See section 4. Cross-references to other
metric ids must exist and be visible (not `hide: true`); a cross-module
reference should say which module it lives in, because the other module may
not be in the package.

**importantNotes (metric and preset).** AI-only. Use it for instructions the
agent must follow to query correctly (which value to select, which preset to
use when). Anything it names (preset ids, tool names, parameter descriptions)
must be valid on both surfaces, or the note must say which surface. It is
shown in the instance language in both tools, so the fr and pt versions
must carry the same instruction and the same ids.

**valueLabelReplacements.** A short noun phrase, because the same string is
the chart legend, table header and map label for humans, and the `(label)`
annotation for the AI. Not a sentence. Must describe the column's meaning
after aggregation (a flag column labelled "binary variable" is wrong once
the app has averaged it into a proportion; label the quantity the user
sees).

**vizPreset `id` and `label.en`.** Together they are all the agent has to
pick a preset. The label should state the figure type, the grain, and what
varies (series, rows, replicant). Follow DOC_VIZPRESET_STANDARDS.md for
wording. `label` is `.en` even on fr/pt instances, so the English must stand
alone.

**Parameter `description`.** It is the form label for humans and the only
name of the parameter an SPA agent ever sees. It must say what the value
does in the script, in words that let an agent map a user's request onto
it. If the default is important, `aiDescription` should quote it.

**Module `label`.** Shown in a 224px sidebar and as a section header; kept
short. The `Mn.` prefix is rendered verbatim everywhere. Flag only clear
inaccuracy or a label that no longer describes the module.

**Tier 2: captions, sub-captions, footnotes, preset descriptions.** Check
that numbers and method statements match the script and parameters, that
`DATE_RANGE` / `PLAGE_DE_DATES` / `INTERVALO_DE_DATAS` are present in
sub-captions (the build already enforces this), and that a footnote does not
describe a computation this module does not perform. Do not review style.

**Writing conventions (all fields, low severity).** No em-dash (the
character —); use a full stop, comma or colon. Prefer the common word.
Explain a technical term a general reader would not know. No rhetorical
contrast, no build-up sentences. Use "Admin Area 2/3/4", "actual vs
expected", "service volume" as DOC_VIZPRESET_STANDARDS.md prescribes.
Consistent ids: refer to dimensions by their `disOpt` id and to metrics by
their id.

## 6. Severity scale

- **S1 Wrong.** The text states something false about what the number is or
  how it is computed, or names a dimension, id, tool or parameter that does
  not exist. An agent following it would misreport.
- **S2 Misleading or unactionable.** Technically defensible but likely to
  lead an agent to the wrong metric, the wrong query, or a claim it cannot
  support. Includes numbers that drift from defaults, instructions that only
  work on one surface, and missing direction or grain.
- **S3 Unclear.** Correct but hard to use: vague, verbose, inconsistent with
  sibling metrics, or repeating the label instead of adding information.
- **S4 Style.** Writing conventions, typos, punctuation.

## 7. Known items

These were reported by an inventory pass before this briefing was written. They
have not all been verified. The lane that owns the module confirms or
rejects each one and includes it in its report with the usual evidence;
nobody should spend effort re-discovering them.

| Where | Reported item |
| --- | --- |
| m001 `m1-01-01` | `valueLabelReplacements.outlier_flag` is a sentence with a typo ("whether this an outlier") and becomes a chart legend. Verified by direct read. Same pattern likely on `m1-02-02`, `m1-04-01`. |
| m001 `m1-03-01` presets | Footnote states a 30% threshold; no parameter carries it. Check `script.R`. |
| m001 `_parameters.ts` vs `script.R` header | `CONSISTENCY_PAIRS_USED` and `DQA_INDICATORS` defaults differ from the header literals. Header is local-dev only (see section 8), so this is stale, not broken. |
| m002 `adjustment-table` footnote | Copies m001's outlier sentence (10 MADs, 80%, 100) into a module with no parameters. |
| m003 `m3-01-01.importantNotes` | Names four preset ids, tells the agent to call `get_module_settings` (SPA only), and quotes the parameter by its English description inside the fr and pt versions. |
| m003 `m3-02-01.caveats` | Names `DIFFPERCENT` by replacement string; the agent only ever sees the description "Difference percent threshold for visualization". |
| m003 `m3-02-02`, `m3-03-02`, `m3-04-02`, `m3-05-02` `typicalRange` | ±10 / ±10-30 / ±10-40 / ±20-50 escalate by admin level; source unknown. |
| m003 `_parameters.ts` vs header | `SELECTEDCOUNT` and `RUN_DISTRICT_MODEL` defaults differ from header literals. |
| m004, m005 `_parameters.ts` vs header | `SELECTED_COUNT_VARIABLE` and `INFANT_MORTALITY_RATE` (0.067 vs 0.063) differ; a stale inline comment `P1_NMR <- 0.039  #Default = 0.03`. |
| m005 | Metric ids are `m4a-*` inside module `m005` with results objects `M5_*`; the agent sees `m4a` ids beside an "M5." module label. |
| m005 `m4a-01-01.typicalRange` | "2-5% / 2-4% of total population" relates to the seven rate parameters; check. |
| m006 `script.R` header | Says `DENOMINATOR_CHAIN` must match m005's setting; m005 has no such parameter. |
| m009 `m9-01-01.caveats` | Suppression rules "<25, <125 fertility, <250 mortality" are hard-coded; check they match the script. |
| m009 `iceh-equiplot`, m012 `scorecard-table` | Footnote is `""` rather than `null`. |
| m011 | No `pt` on any field, module label included. Falls back to en. Ten fields hard-code "95% CI". |
| m012 `m12-01-01` | `disaggregationGuidance` points to `m3-01-01` (another module) for facility level; `typicalRange` uses an em-dash. |
| all modules | `caveats` is null on 15 of 48 metrics; every other `aiDescription` field is always present. Judge whether a null caveat is a gap for that metric. |
| fr/pt | Translations are markedly shorter than the English on m001 (fr only), m002, m004 `m4-01-01`, m005 `m4a-01-01`, m006 `m6-01-01`. Only relevant to fields the AI reads in the instance language (label, variantLabel, importantNotes) and to Tier 2. |
| authoring docs | `DOC_METRICS_REFERENCE.md` shows a `postAggregationExpression` written as `cases / NULLIF(tests, 0)` while `_metrics_template.ts` requires `name = arithmetic` with automatic NULLIF; its disaggregation table omits most of `disaggregation_options.ts`; `DOC_MODULES.md` describes a single `_metrics.ts` while modules use `_metrics/`. |

## 8. Pitfalls that produce false findings

- **The script header above `#---` is not executed.** wb-fastr strips every
  line above the first `#---` line and substitutes tokens in what remains.
  Values assigned above the marker are local-development stand-ins.
  A mismatch with `_parameters.ts` defaults is stale documentation
  (S4), never a runtime defect. The body below the marker is the truth.
- **`definition.json` is generated.** Cite the `.ts` source for locations.
  The orchestrator has already confirmed the built output matches the
  sources; if you see a mismatch anyway, report it at the top of your file
  and review `definition.json` as shipped.
- **`quarter_id` is deliberately hidden** from the AI's disaggregation
  lists. Its absence is not a finding.
- **Facility options depend on the instance.** `facility_type` and
  `facility_ownership` appear only when the results object has
  `facility_id` and the instance enables the column. Guidance that mentions
  them is fine if the results object carries `facility_id`.
- **The catalog's duplicate-label check is English only.** Two metrics may
  share a French label without a variant and no build error fires. Only
  flag if the fr or pt labels collide.
- **Preset captions never reach the AI.** A caption that disagrees with the
  summary is Tier 2, not Tier 1.
- **`hide: true` metrics** are absent from the catalog but still queryable
  by id, and a preset on a hidden metric with
  `createDefaultVisualizationOnInstall` still produces a default figure.
  Cross-references to a hidden metric are S1.
- **Sample size and roll-up.** The AI's CSV excludes the admin-area roll-up
  row, and HFA values carry `(n=…)` inline. Text need not explain these; the
  tool does.
- **Consumer-side limits are not module findings.** That `description` is
  never shown, that `summary` and `caveats` never co-occur, that `.en` is
  hard-coded: record them under "Observations about the consumer" so the
  owner can decide, and write the module text to work within them.

## 9. Files to read, per lane

Modules repo (this repo):

- `mNNN/_core.ts`, `_parameters.ts`, `_results_objects.ts`,
  `_metrics/*.ts`, `script.R`, `definition.json`
- `_shared/text_presets.ts` (m004 and m006 import footnotes from it)
- `.validation/_module_definition_github.ts` (the schema),
  `.validation/disaggregation_options.ts`
- `DOC_METRICS_REFERENCE.md`, `DOC_VIZPRESET_STANDARDS.md`,
  `DOC_MODULES.md` (authoring conventions)
- the rendered AI view produced by `render_ai_view.ts`

wb-fastr (read-only reference; paths relative to that repo):

- `lib/ai_tools/format_metrics_list_for_ai.ts`
- `lib/ai_tools/format_metric_data_for_ai.ts`
- `lib/ai_tools/build_system_prompt.ts` (the generic interpretation rules
  the agent is given, so module text need not repeat them)
- `server/runs/disaggregation_availability.ts`
- `server/module_loader/load_module.ts` (translation and validation)
- `lib/indicator_format_metrics.ts` (formatAs override list)
- `client/src/components/project_ai/ai_tools/tools/_internal/format_module_settings_for_ai.ts`
- `client/src/components/project_ai/ai_tools/tools/_internal/format_modules_list_for_ai.ts`
- `SYSTEM_13_ai_assistant.md` (architecture and rulings)
