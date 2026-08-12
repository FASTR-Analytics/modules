# PLAN: M4/M5/M6 admin3 outputs must carry `admin_area_2`

Status: proposed 2026-08-12, not started. Raised by an adversarial review of
the app-side per-family structure split (wb-fastr commit `3d320bb4`).

## Problem

Four results objects report at admin-area-3 granularity but drop the parent
`admin_area_2` column:

| Module | Results object |
| --- | --- |
| m004 | `M4_coverage_estimation_admin_area_3.csv` |
| m005 | `M5_denominators_admin3.csv` |
| m005 | `M5_combined_results_admin3.csv` |
| m006 | `M6_coverage_estimation_admin3.csv` |

The app scopes a project to a single admin-area-2 (national projects are
unaffected). When a results object carries `admin_area_2`, the scope filter is
applied **directly** to that column. When it carries only `admin_area_3`, the
app must instead derive the set of AA3 names belonging to the project's AA2
from the run's facilities parquet, and filter by that derived list.

That derivation needs to know which facility registry (HMIS or HFA) the
results object belongs to — and for m004/m005/m006 it cannot, because every
`dataSource` on those modules is `sourceType: "results_object"` (they read
upstream module outputs, never a dataset directly), so the app's
`getDatasetFamily` returns `undefined`.

Two consequences on the app side, neither fixable in the modules repo alone:

1. The app now **fails closed** — those four results objects render empty in
   an AA2-scoped project rather than silently showing every district in the
   country. Correct, but the data is unavailable.
2. Even where derivation does work, it matches AA3 areas **by name**, so an
   instance with duplicate district names across regions folds the twin's
   numbers in (a known accepted latent, measured nil in production today).

## The change

Add `admin_area_2` to those four outputs, alongside the existing
`admin_area_3`:

- `script.R` in each module: carry the parent AA2 through the aggregation
  rather than dropping it. **Derive it from the input data, never hardcode**
  — the input CSV carries admin columns only down to the family's configured
  depth, and that depth is now per-registry.
- `_results_objects.ts` in each module:
  `createTableStatementPossibleColumns` gains `admin_area_2: "TEXT NOT NULL"`
  for those results objects. The app hard-errors on any CSV header that is
  not declared, so the declaration and the script must land together.
- `deno task build` to regenerate each `definition.json`.

## Why this is the durable fix

With `admin_area_2` present, those results objects take the direct-filter
path. The family no longer needs deriving, the name-matched AA3 derivation is
retired for them, and the duplicate-district-name latent disappears for these
outputs. It also removes the only case where an admin-bearing results object
cannot be scoped at all.

## Rollout notes

- **Packages are immutable.** This fixes results objects in packages
  generated *after* it ships; existing attached packages keep the fail-closed
  behaviour until their project is repointed to a newly generated package.
  That is expected, and needs no app change.
- Ships in lockstep with the app the usual way (deploy the app first, then
  push the modules repo), though nothing in the app has to change for this:
  the app already handles both column shapes.
- No `inputKey` implications for other modules — the change is additive to
  these results objects' schemas only.

## Verification

- Run each module and confirm the four CSVs carry a non-null `admin_area_2`
  on every row, consistent with the `admin_area_3` value's parent in the
  input facilities data.
- In the app, attach a newly generated package to an AA2-scoped project and
  confirm the four visualizations render that AA2's districts only.
