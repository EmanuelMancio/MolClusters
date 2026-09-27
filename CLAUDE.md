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
- `ClusterTracker` (`tracker.py`) owns the clusters and their ids, with no analysis or file
  output; `MolClusters` runs the analyses on top of it (`MolClusters.tracker`). The tracker keeps
  cluster identity stable across formation/merge/split events with a two-step assignment in
  `update`: each previous cluster picks its best connected group (`_best_groups`: most molecules,
  then purest), then each group continues its largest contributor (`_get_older_cluster`: then
  oldest id) or becomes a new cluster.
- The analyses on top of the tracker are being moved into the `analysis/` package, one
  `FrameAnalysis` subclass per output: `prepare(run)`, `analyse(frame)` on every frame including
  frame 0, `finish(run)`, all called in list order from `MolClusters.run()`. They see the run
  through `Run` (universe, config, `output`, `analysis(Type)` lookup of earlier analyses) and each
  frame through `Frame` (index, time, read-only clusters, `find`), never the tracker itself. All
  file writes go through `Run.output` (`output.py`'s `RunOutput`: `path(name)` for whole files,
  buffered `append` for per-frame ones), and each analysis declares its files in `outputs`
  (`OutputFile`), which drives the earlier-run check and the end-of-run summary. Library users add
  their own analyses with `MolClusters(..., analyses=[...])`, which run after the config-enabled
  built-ins (there is no registry yet); `Run` dispatches every hook and notes which analysis raised
  an error (`_blame`), and the CLI prints those notes. So far only
  `SizeEvolution` (evo.txt) has moved; solute-solvent, coordinates, nucleus and the JSON report
  are still methods of `MolClusters`.
- The analysis is frame-oriented: `MolClusters.run()` calls `tracker.update()` frame by frame,
  mutating in-memory cluster state, and only writes results to disk once the full run finishes.
- The config's `distance_backend` ("serial"/"OpenMP") only accelerates `ConnectionTable`'s "cm"
  rule (`capped_distance`/`self_capped_distance`); MDAnalysis silently ignores it under its
  auto-selected "nsgrid" method (typical when the cutoff is much smaller than the box), and "hb"
  rules always run serially since `HydrogenBondAnalysis` doesn't expose a backend option. Even
  when the backend is actually used, benchmarking found it delivers zero speedup on the official
  PyPI Windows wheel (confirmed on 2.10.0): its `c_distances_openmp` extension imports without
  error but was compiled with OpenMP disabled (`MDAnalysis.lib.distances.USED_OPENMP is False`),
  so it silently runs the same serial code. `ConnectionTable` warns (via `_warn_if_openmp_unavailable`)
  when that's detected.

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
