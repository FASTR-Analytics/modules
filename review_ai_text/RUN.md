# Run the AI-text review

This directory is a durable, repeatable review of the human-authored text in
this repo that AI agents consume. To run it, point an agent at this file and
say "run this". The result is a report. A run adds files under
`review_ai_text/reports/` and leaves every tracked file as it found it.

Files:

- `RUN.md` (this file): the orchestrator's steps.
- `BRIEFING.md`: context every lane must read. Which fields reach an AI, how,
  the quality criteria, the severity scale, known items, pitfalls.
- `LANES.md`: the prompt for each lane agent, with placeholders.
- `render_ai_view.ts`: renders, from each `definition.json`, the text an AI
  receives. Lanes review that output rather than reconstructing it.
- `reports/<date>/`: created by a run. Per-lane reports plus the
  consolidated `REPORT.md`.

## Placeholders used below and in LANES.md

| Placeholder | Meaning |
| --- | --- |
| `REPO` | Absolute path of the modules repo root, the parent of this directory |
| `APP` | Absolute path of the wb-fastr repo root |
| `REPO_HEAD`, `APP_HEAD` | Commit hashes recorded in step 0.4 |
| `REPORTS` | Absolute path of `REPO/review_ai_text/reports/<date>` |
| `AI_VIEW` | `REPORTS/ai_view` |
| `DELTA` | `REPORTS/briefing_delta.md`, or the word `none` |
| `MODULE` | A module directory name such as `m003` |
| `OUT` | The lane's report file, always inside `REPORTS` |

## Orchestrator steps

You are the orchestrator. You launch lanes, you do not do lane work
yourself, and you create or change files only inside `REPORTS`. The one
command that touches tracked files is the build in step 0.5, and it must
leave them unchanged.

### 0. Prepare

1. Read `BRIEFING.md` in full.
2. Set `REPO`. Set `APP` to the sibling directory `wb-fastr`. If it is not
   there, ask the user for its path. If the user says it is unavailable,
   skip the `APP` parts of steps 4 and 6, pass `APP` as the word `none`, and
   record "APP not available" in the report header.
3. Confirm the working tree is clean: `git -C REPO status --porcelain` must
   print nothing except untracked paths under `review_ai_text/reports/`.
   Otherwise stop and tell the user.
4. Record `REPO_HEAD` = `git -C REPO rev-parse HEAD` and `APP_HEAD` =
   `git -C APP rev-parse HEAD`.
5. Check the definitions are fresh: run `deno task build` in `REPO`, then
   `git -C REPO status --porcelain`. If any `definition.json` changed, run
   `git -C REPO checkout -- '*/definition.json'` to restore them, then stop
   and tell the user to rebuild and commit before running the review. The
   review only ever runs on committed definitions that match their sources.
6. Create `REPORTS` (ISO date; add a suffix if it exists).
7. Check the briefing's account of the consumer still holds. Open the
   wb-fastr files listed in section 9 of the briefing and compare with
   section 3. If a field's treatment has changed (for example
   `vizPreset.description` is now rendered, or fr and pt are now used),
   write the differences to `REPORTS/briefing_delta.md` and set `DELTA` to
   that path; otherwise set `DELTA` to `none`. Do not edit the briefing
   during a run; propose the update in the final report.
8. Render the AI view:

   ```sh
   cd REPO && deno run --allow-read --allow-write review_ai_text/render_ai_view.ts REPORTS/ai_view
   ```

   Read `AI_VIEW/catalog.md` once yourself so you can judge the lane reports.

### 1. Launch the module lanes

One agent per active module directory: every `mNNN` not listed in
`_frozen_modules.ts` (m001 to m006 and m009 to m012 at the time of
writing). Use the general-purpose agent type, run in the background, all at
once. For each, the prompt is the common preamble from `LANES.md` followed
by the Lane M section, with every placeholder filled. `OUT` is
`REPORTS/lane_<MODULE>.md`.

If the harness limits parallel agents, run in batches of five.

### 2. Launch Lane C and Lane X

When every module lane has written its report, launch Lane C (`OUT` =
`REPORTS/lane_catalog.md`) and Lane X (`OUT` = `REPORTS/lane_crosscut.md`)
together. Both read the module lane reports in `REPORTS`.

Launch Lane T (`OUT` = `REPORTS/lane_translations.md`) only if the user
asked for translations. It can run alongside the module lanes.

### 3. Consolidate

First, spot-check every S1 finding yourself: open the evidence it cites and
confirm the code says what the finding claims. Do not confirm by re-reading
the lane's prose. Mark each S1 as "orchestrator-confirmed" or "not
confirmed".

Then write `REPORTS/REPORT.md`:

1. **Header.** Date, `REPO_HEAD`, `APP_HEAD`, which lanes ran, whether a
   briefing delta existed, whether `APP` was available.
2. **Summary.** Counts by severity, and the five to ten findings that most
   change what an agent would say to a user, each in one sentence with its
   finding id.
3. **Findings by module.** For each module, the module lane's findings
   merged with anything Lane C or Lane X raised about that module. When two
   lanes raise the same issue, keep one entry and cite both finding ids.
   Keep the full finding record (severity, confidence, location, current,
   problem, evidence, proposed, translations affected).
4. **Cross-cutting findings.** Lane X's items that span modules, and Lane
   C's scenario table.
5. **Authoring-doc findings.** Lane X step 8, including the proposed
   "writing for the AI" section.
6. **Coverage.** A table of metric by field marking checked-accurate,
   finding, or not verifiable, built from the lanes' "Verified accurate"
   and "Not verifiable" sections. State every gap in coverage.
7. **Observations about the consumer.** Merged and deduplicated from all
   lanes. Labelled as facts about wb-fastr, not module findings.
8. **Known items disposition.** Each row of briefing section 7 with
   confirmed, rejected, or not examined, and the finding id if confirmed.
9. **Proposed briefing updates.** Anything a lane or you found that should
   change `BRIEFING.md` before the next run (a rule that no longer holds, a
   new known item, a pitfall that produced false findings).

### 4. Hand over

Tell the user: the path of `REPORT.md`, the severity counts, the top
findings in one line each, and whether any lane failed or any coverage gap
exists. Do not commit. Do not edit module files. The user decides what to
apply.

## Rerunning

Each run is a fresh `reports/<date>/` directory. The lanes do not read
previous reports, so a rerun is a clean review. If the user wants a delta
instead, add to every lane prompt: "Read `<previous REPORT.md>`. For each
earlier finding in your scope, mark it fixed, unchanged, or superseded, in
addition to your normal output."

Between runs the user applies the "Proposed briefing updates" from the last
report to `BRIEFING.md`, in its own commit.
