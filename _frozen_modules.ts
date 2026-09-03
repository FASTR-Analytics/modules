// FROZEN, not authored (ruled 2026-09-01, closed — never re-raise): m007 and m008 are
// retired, but every still-deployed pre-restructure app resolves its WHOLE
// registry — m007 and m008 included — from this repo's HEAD at wizard time,
// so their committed files must stay on main, byte-frozen, until the fleet
// runs the new app. They are exempt from the build, the type check, and
// definition validation (m008's "calculated_indicators" generation type no
// longer validates under the current authoring schema, by design). PLAN_1e
// (module cleanup) deletes the directories AND this file in the same
// commit. Never rebuild or edit them.
export const FROZEN_MODULE_DIRS = ["m007", "m008"];
