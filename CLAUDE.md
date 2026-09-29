<!--
SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>

SPDX-License-Identifier: LGPL-3.0-or-later
-->

# CLAUDE.md

## Commands

- Run all pre-commit hooks (ruff check --fix, ruff format, REUSE license lint, pytest, and
  `cz check` commit-message lint): `uv run pre-commit run --all-files`. The pytest hook runs
  only the tests affected by what changed (`--testmon --no-cov`, blind to data files); the
  whole suite is `uv run pytest`. All hooks are `local`
  and run the uv-installed tools, so their versions come from `uv.lock` (bumped weekly by
  Dependabot, `.github/dependabot.yml`); don't use `pre-commit autoupdate`.
- Docs (MkDocs Material + mkdocstrings, `docs/`, `mkdocs.yml`, `docs` dependency group):
  `uv run --group docs mkdocs serve` to preview, `uv run --group docs mkdocs build --strict`
  as CI does (`.github/workflows/docs.yml`, which deploys `main` to GitHub Pages). A change to
  a config key, the CLI help, an output file, the report format, a property or id assignment
  updates its user-guide page too (table in `docs/developer-guide/docs.md`).

Note: pytest config lives solely in `pytest.ini` (`addopts` already includes `--cov=src
--cov-report html`), so a plain `uv run pytest` produces an HTML coverage report under `htmlcov/`.

## Design notes

The architecture is documented in `docs/developer-guide/` (`architecture.md`, `tracking.md`,
`builtin-analyses.md`, `report-format.md`); read the relevant page before changing a part.
What's easy to get wrong:

- `ClusterTracker` (`tracker.py`) owns clusters and ids; `MolClusters` only runs it and the
  analyses. Analyses (`analysis/`, one `FrameAnalysis` per output) see the run through `Run`
  and each frame through `Frame`, never the tracker. Clusters and groups are read-only outside
  the tracker (`Cluster._update` is its only mutator; returned arrays are read-only).
- Every file write goes through `Run.output` (`RunOutput`), and each analysis declares its files
  in `outputs` (`OutputFile`), which drives the earlier-run checks and the end-of-run summary.
- A new built-in is one line in `analysis/builtins.py`'s `BUILTINS`, in order; `JsonReport`
  runs `last` and knows no other analysis: others add to the report through
  `report_frame`/`report_cluster`. Changing the report schema (in `report.py`'s docstring)
  bumps `VERSION`; floats are rounded to `DECIMALS` where they're computed.
- LAMMPS: residue names come from the config's `lammps_resnames`; dumps record steps, so
  `lammps_timestep` is passed as `dt`, which is the MD timestep, not the frame spacing.
- `--traj-memory` must not use MDAnalysis' `in_memory`/`transfer_to_memory` (they retime frame
  i as i x dt); `main._load_into_memory` copies frames and times into a `_TimedMemoryReader`.
- `distance_backend` "OpenMP" gives no speedup on the PyPI Windows wheels (built without
  OpenMP, `USED_OPENMP is False`); `ConnectionTable` warns about it.

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
- Units (amu, e, D, Å, Å³, g/cm³, degrees, ps) are pinned in `tests/test_units.py`, each against
  a physical fact rather than a copy of its formula; a new reported quantity or output file gets
  a test there too.
