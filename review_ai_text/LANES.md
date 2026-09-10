# Lane prompts

The orchestrator (see RUN.md) launches the module lanes in parallel, then
Lane C and Lane X, each with the common preamble followed by its lane
section, with every placeholder filled in from the table in RUN.md. Every
lane is read-only with respect to the repos.

Finding ids are fixed per lane so the consolidated report can cite them
without collisions: `F-{MODULE}-<n>` for a module lane (for example
`F-m003-4`), `F-C-<n>` for Lane C, `F-X-<n>` for Lane X, `F-T-<n>` for
Lane T.

---

## Common preamble (prepend to every lane)

You are one lane of a review of the text in the wb-fastr-modules repo that
AI agents consume. Your job is to find inaccuracies and unclear wording and
to propose replacement text. You do not edit or create any file in `{REPO}`
or `{APP}` other than your report `{OUT}`.

Paths: modules repo `{REPO}` at commit `{REPO_HEAD}`; wb-fastr `{APP}` at
commit `{APP_HEAD}` (the word `none` means it is unavailable; rely on the
briefing); rendered AI views in `{AI_VIEW}`; all lane reports in
`{REPORTS}`.

Before anything else, read `{REPO}/review_ai_text/BRIEFING.md` in full. It
defines which fields reach an AI, what "accurate" means, the severity
scale, the known items, and the pitfalls that produce false findings.
Follow it. If `{DELTA}` is a path rather than `none`, read that file next;
where it differs from the briefing, it wins.

Rules:

1. Verify against code, not against other prose. A methodology claim is
   checked against `script.R` and the metric's `valueFunc` or expression. A
   number is checked against `_parameters.ts` defaults or a constant in the
   script body (below the `#---` marker). A dimension is checked against the
   results object's declared columns. A metric id is checked to exist and to
   be visible (`hide: false`).
2. Every finding cites its evidence with file and line numbers.
3. Quote the current text verbatim. Propose replacement English text for S1
   to S3 findings; a proposal is not required for S4.
4. Do not report consumer-side behaviour of wb-fastr as a finding. Put it
   under "Observations about the consumer".
5. Do not re-run the inventory. Confirm or reject the known items that fall
   in your scope and include them with evidence.
6. Prefer fewer, well-evidenced findings over many speculative ones. Mark
   confidence honestly.

Report format for `{OUT}` (Markdown):

```text
# Lane <name> report
Commit: {REPO_HEAD} / {APP_HEAD}
Scope: <what you reviewed>
Files read: <list>

## Summary
<counts by severity; two or three sentences on the main problems>

## Findings
### F-<lane>-<n>: <metric, preset or parameter id> / <field>
- Severity: S1 | S2 | S3 | S4
- Confidence: high | medium | low
- Location: <file>:<line> (<field path, e.g. aiDescription.caveats.en>)
- Current: "<verbatim>"
- Problem: <one or two sentences>
- Evidence: <file>:<lines> and what they show
- Proposed: "<replacement English text>" | delete | no proposal
- Translations affected: <fr, pt, none> and whether the existing fr/pt would need redoing
(repeat)

## Verified accurate
<one line per metric/field pair you checked and found correct, so coverage is visible>

## Not verifiable
<anything you could not establish, and why>

## Observations about the consumer
<facts about how wb-fastr renders or omits these fields that limit what the text can achieve>

## One-line truth per metric   (Lane M only)
<for each visible metric: one sentence stating what one displayed value is:
the quantity, the grain it was computed at, the aggregation the app applies,
the unit>
```

---

## Lane M: one module ({MODULE})

Scope: every AI-reaching text field of `{MODULE}` (Tier 1) plus a light
Tier 2 pass over its preset captions and footnotes.

Procedure:

1. Read `{AI_VIEW}/{MODULE}.md`. This is what the AI sees, plus the ground
   truth it does not see, per metric. Read `{AI_VIEW}/modules.md` for the
   module label and the parameter descriptions as the SPA agent sees them.
2. Read `{REPO}/{MODULE}/_results_objects.ts`, `_parameters.ts`,
   `_core.ts`, and every file in `_metrics/`.
3. Read `{REPO}/{MODULE}/script.R`. For each results object the metrics use,
   locate where the script writes it and trace how each `valueProp` column
   and each ingredient column is computed, at what grain (facility-month,
   area-month, area-year, national-year), and under which parameters. Note
   every parameter's actual effect in the body (below `#---`). If the
   module has prerequisites, read only as much of the upstream module's
   `_results_objects.ts` as needed to know what the inputs are.
4. For each visible metric, work through the fields in the order of section
   5 of the briefing: catalog line, methodology, interpretation, typical
   range, caveats, disaggregation guidance, importantNotes, value labels,
   preset id and label, preset importantNotes. For each, decide: accurate
   against step 3, or not; clear for an agent that has nothing else, or not.
   Check every number against `_parameters.ts` defaults or script constants.
   Check every named dimension against section 4 of the briefing. Check
   every metric or preset id named in text exists and is visible.
5. For hidden metrics, check only that nothing visible references them.
6. Parameters: check each `description` says what the parameter does in the
   script body, and that `aiDescription` text quoting a default matches
   `_parameters.ts`. Check `select` option labels are meaningful to a person
   choosing among them.
7. Tier 2 pass: for each preset, check caption, subCaption and footnote
   numbers and method statements against step 3. Do not review style.
8. For the fields the AI reads in the instance language (label,
   variantLabel, importantNotes), check that the fr and pt versions carry
   the same meaning and the same ids as the en. Do not review the fr and pt
   of `aiDescription`; the AI never reads them.
9. Consistency within the module: sibling metrics across admin levels
   should have parallel labels, variants, and text that differs only where
   the computation differs.
10. Confirm or reject the known items for `{MODULE}` from section 7 of the
    briefing.
11. Write `{OUT}`, including the "One-line truth per metric" section. Lane
    C relies on that section's exact heading.

---

## Lane C: the catalog as a whole

Scope: `{AI_VIEW}/catalog.md`, the full listing exactly as an agent receives
it from `get_available_metrics`, judged as one document. Runs after the
module lanes, and reads the "One-line truth per metric" section of every
`{REPORTS}/lane_m*.md`.

The question for this lane: given a realistic user request, can an agent
pick the right metric, values, disaggregations and preset from the catalog
alone, and does the summary tell the truth about what the number is?

Procedure:

1. Read the briefing sections 3 and 5, then read `catalog.md` end to end
   once, as an agent would, before checking anything.
2. Build a table: for each metric, the one-line truth from the module lane
   beside the catalog's `summary`. Flag any summary that disagrees with the
   truth (S1), omits grain or unit where that changes interpretation (S2),
   or merely restates the label (S3).
3. Sibling and near-duplicate test. For each group that shares a label, or
   whose labels differ only in wording, state what an agent would need to
   know to choose between them and whether the catalog line provides it.
   Groups to examine at least: m3-02, m3-03, m3-04, m3-05 (four admin
   levels of the same three metrics); `m4-*`, `m4a-*` and `m6-*` (three
   modules that each compute coverage); m10-01 and m10-03 (observed versus
   carry-forward, indicator versus variant); m11 versus m3-02 (two
   disruption methods); m12-01-01 versus m3-01-01 and `m6-*` (indicator
   values versus volume versus coverage).
4. Scenario test. For each of these user requests, write which metric and
   preset an agent should choose and whether the catalog text leads there:
   "how complete is reporting by district last year", "which indicators had
   outliers", "show ANC1 coverage trend nationally", "which regions had
   service disruptions in 2025", "compare adjusted and unadjusted volumes",
   "what is the neonatal mortality rate", "what did the facility survey say
   about stock-outs", "is the urban-rural gap in immunisation widening".
   Add two more of your own. A request the catalog cannot route correctly is
   a finding against the metrics involved.
5. Preset lines. For each preset, judge whether `id` plus `label.en` plus
   the derived notes tell an agent the figure type, grain, what varies, and
   when to use it instead of its siblings. Check every metric-level and
   preset-level `importantNotes` that appears in the catalog for validity on
   both surfaces (at the time of writing, only the metric note on m3-01-01).
6. Terminology and format across the whole catalog: `Admin Area N` versus
   "region" or "district", consistent naming of the count variables,
   consistent bracket conventions in variants, formatAs consistent with
   what the summary implies.
7. Length. The whole catalog is sent to the agent in one message, so length
   costs. Flag summaries longer than one sentence and notes longer than
   needed.
8. Write `{OUT}`. Include the scenario table in full.

---

## Lane X: cross-cutting consistency and authoring docs

Scope: things that only show when modules are compared, plus the authoring
documents that future authors and agents will follow. Runs after the module
lanes; step 7 reads their reports in `{REPORTS}`.

Procedure:

1. Read the briefing, then `{AI_VIEW}/catalog.md` and `{AI_VIEW}/modules.md`.
2. Parameter descriptions across modules. `SELECTEDCOUNT` (m003, m011,
   m012) and `SELECTED_COUNT_VARIABLE` (m004, m005) name the same choice;
   check their descriptions and option labels agree, and that the option
   labels are understandable by someone who has not read m002. Same for
   any other parameter that recurs.
3. `valueLabelReplacements` across modules: the same column name
   (`count_final_*`, `coverage_cov`, `value`, `estimate`) should carry the
   same label wherever it means the same thing, and different labels where
   it does not.
4. Cross-metric references. List every metric id mentioned inside any text
   field of any module (`grep -o 'm[0-9]\+a\?-[0-9]\+-[0-9]\+'` over the
   `_metrics` directories, then filter to prose fields). Check each target
   exists, is visible, and is in the same module or the text says which
   module it is in.
5. Terminology. Grep all en text for "region", "district", "province",
   "sub-district", "service delivery volume", the em-dash character, and
   "Admin Area" with inconsistent casing. Report as S4 with locations, in
   one finding per term.
6. Language coverage of the fields the AI reads in the instance language
   (label, variantLabel, importantNotes, module label, parameter
   descriptions): list every missing `pt` and any fr or pt that carries
   different ids or numbers than the en. m011 is known to lack pt entirely.
7. Shared text presets (`_shared/text_presets.ts`): check the footnotes'
   method statements against what m004 and m006 actually compute, using
   the "One-line truth per metric" sections of `{REPORTS}/lane_m004.md`
   and `lane_m006.md` rather than re-reading the R. Note which presets use
   the shared text and which have hand-written variants that have drifted.
8. Authoring documents. Check `DOC_METRICS_REFERENCE.md`,
   `DOC_VIZPRESET_STANDARDS.md`, `DOC_MODULES.md`, and `_metrics_template.ts`
   against the schema in `.validation/_module_definition_github.ts`, the
   build checks in `build_definitions.ts`, and section 3 of the briefing.
   Report every statement that is out of date or that would lead an author
   to write text the AI never sees or that the build rejects. Propose, as
   one finding, a short "writing for the AI" section for
   `DOC_METRICS_REFERENCE.md` that states which fields the AI reads and the
   per-field criteria from the briefing, so future authors get it right at
   authoring time.
9. Write `{OUT}`.

---

## Lane T (optional): translation fidelity

Run only if the orchestrator was told to include translations. Scope: the
fr and pt of every field an AI or a person reads, for all modules. Read the
`_metrics/*.ts`, `_parameters.ts` and `_core.ts` sources directly;
`{AI_VIEW}` is not needed.

Procedure: for each metric, preset and parameter, compare fr and pt against
en for meaning, ids, numbers and completeness. Report omissions (a sentence
in en absent from fr), contradictions, and ids or numbers that differ. Do
not report style. Group findings by module. Use severity S1 for a
contradiction or a wrong id, S2 for an omitted instruction or number, S3
for an omitted explanatory sentence. Write `{OUT}`.
