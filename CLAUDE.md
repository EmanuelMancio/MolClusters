<!--
SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>

SPDX-License-Identifier: GPL-3.0-only
-->

# CLAUDE.md

## Commands

- Run all pre-commit hooks (ruff check --fix, ruff format, REUSE license lint, pytest, and
  `cz check` commit-message lint): `uv run pre-commit run --all-files`. All hooks are `local`
  and run the uv-installed tools, so their versions come from `uv.lock` (bumped weekly by
  Dependabot, `.github/dependabot.yml`); don't use `pre-commit autoupdate`.

Note: pytest config lives solely in `pytest.ini` (`addopts` already includes `--cov=src
--cov-report html`), so a plain `uv run pytest` produces an HTML coverage report under `htmlcov/`.

## Design notes

- LAMMPS topologies carry no residue names: the config's `lammps_resnames` (name -> LAMMPS
  molecule id or `"first-last"` id range) supplies them, applied in `main._apply_lammps_resnames`.
- `MolClusters` keeps cluster identity stable across formation/merge/split events with a
  two-step assignment in `__update_clusters`: each previous cluster picks its best connected group
  (`__best_groups`: most molecules, then purest), then each group continues its largest
  contributor (`__get_older_cluster`: then oldest id) or becomes a new cluster.
- The analysis is frame-oriented: `MolClusters.run()` mutates in-memory cluster state frame by
  frame and only writes results to disk once the full run finishes.
- The config's `distance_backend` ("serial"/"OpenMP") only accelerates `ConnectionTable`'s "cm"
  rule (`capped_distance`/`self_capped_distance`); MDAnalysis silently ignores it under its
  auto-selected "nsgrid" method (typical when the cutoff is much smaller than the box), and "hb"
  rules always run serially since `HydrogenBondAnalysis` doesn't expose a backend option.

## Conventions

- Commit messages must follow Conventional Commits using the Angular convention (types `feat`,
  `fix`, `docs`, `style`, `refactor`, `perf`, `test`, `build`, `ci`, `chore`, `revert`, with an
  optional scope, e.g. `fix(release): ...`), enforced by pre-commit + commitizen; commitizen uses
  them to drive semver bumps and `CHANGELOG.md` generation.
- Every file needs a REUSE/SPDX license header (see `REUSE.toml`, `LICENSES/`), enforced by
  pre-commit's `reuse-lint-file` hook.
- A change to clustering, id assignment, or any other computed output is a breaking change
  (`type(scope)!:` + `BREAKING CHANGE:` footer, even when `type` is `fix`) only if it changes a
  result that was previously *valid/expected* — e.g. the dominance algorithm changing which
  plausible cluster-id assignment you get. Correcting output that was simply wrong (a bug nobody
  would have expected or relied on, e.g. deriving correct atom names for LAMMPS topologies) is a
  normal `fix`, even though the stored values change.
