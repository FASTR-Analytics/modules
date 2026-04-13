# Git Hooks

This directory contains shared git hooks for the repo.

## Setup (one-time, per developer)

```sh
git config core.hooksPath .githooks
```

Run this once after cloning or pulling this directory.

## What the hooks do

### `pre-push`

Before any push, runs the `definition.json` validator. If errors are found, they are displayed and you are asked whether to push anyway:

```
OK    m001/definition.json
OK    m002/definition.json
FAIL  m003/definition.json
      labels[0].id: Required
...
-------------------------------------------------------
  VALIDATION ERRORS
-------------------------------------------------------
FAIL  m003/definition.json
      labels[0].id: Required
-------------------------------------------------------

Push anyway? (y/N)
```

- Press `y` to push despite the errors (CI will also flag them)
- Press `N` or Enter to cancel the push and fix the errors first

Requires [Deno](https://deno.com) to be installed. If Deno is not found, the hook skips validation with a warning.
