// FROZEN, not authored (ruled 2026-09-01, closed, never re-raise): these
// modules are retired (m003 and m004 superseded by m012, m011, m005 and
// m006; m007 and m008 by the restructure), but every still-deployed
// pre-version2 app resolves its WHOLE registry from this repo's HEAD at
// wizard time, so their committed files must stay on main, byte-frozen,
// until the fleet runs the new app. They are exempt from the build, the type
// check, and definition validation (m008's "calculated_indicators" generation
// type no longer validates under the current authoring schema, by design;
// m003 and m004 declare no family, tier or sortOrder). Deleting the
// directories deletes this file in the same commit. Never rebuild or edit
// them.
export const FROZEN_MODULE_DIRS = ["m003", "m004", "m007", "m008"];
