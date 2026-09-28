<!--
SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>

SPDX-License-Identifier: LGPL-3.0-or-later
-->

# CLAUDE.md

## Commands

- Run all pre-commit hooks (ruff check --fix, ruff format, REUSE license lint, pytest, and
  `cz check` commit-message lint): `uv run pre-commit run --all-files`. All hooks are `local`
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

- LAMMPS topologies carry no residue names: the config's `lammps_resnames` (name -> LAMMPS
  molecule id or `"first-last"` id range) supplies them, applied in `main._apply_lammps_resnames`.
- LAMMPS dumps record step numbers, not times, and MDAnalysis times them as step x dt (dt = 1
  unless given). The optional `lammps_timestep` (e.g. `"2 fs"`) is passed as `dt` by
  `main._lammps_dump_timestep`, which warns when a dump is read without it. A dump reader's `dt`
  is the MD timestep, not the time between frames, so `main._log_system` works from the frames'
  own times instead.
- `--traj-memory` doesn't use MDAnalysis' `in_memory`/`transfer_to_memory`, whose in-memory
  reader times frame i as i x dt from 0 (losing start times, uneven spacing, and a dump's
  spacing): `main._load_into_memory` copies the frames, times included, into a
  `main._TimedMemoryReader`. MDAnalysis' MemoryReader applies transformations once to the stored
  frames, not per read, so a transformation can't restore times.
- `ClusterTracker` (`tracker.py`) owns the clusters and their ids, with no analysis or file
  output; `MolClusters` is only the runner on top of it. The tracker keeps
  cluster identity stable across formation/merge/split events with a two-step assignment in
  `update`: each previous cluster picks its best connected group (`_best_groups`: most molecules,
  then purest), then each group continues its largest contributor (`_get_older_cluster`: then
  oldest id) or becomes a new cluster. Each `update` (and the first frame) leaves a
  `Transition` (`tracker.transition`, `Frame.transition` for analyses): `flows` (previous id,
  current id) -> molecules, id 0 for none, `born`, `merged` (absorbed -> survivor) and
  `dissolved`. `Lineage` turns it into events stamped with the first frame that shows them, so a
  cluster's DeathTime is the first frame it's gone from and its lifetime a whole number of frame
  spacings (a one-frame cluster lived one spacing); lifetimes born at the start or alive at the
  end are censored.
- The analyses on top of the tracker live in the `analysis/` package, one `FrameAnalysis`
  subclass per output: `prepare(run)`, `analyse(frame)` on every frame including
  frame 0, `finish(run)`, all called in list order from `MolClusters.run()`. They see the run
  through `Run` (universe, config, `output`, `analysis(Type)` lookup of earlier analyses) and each
  frame through `Frame` (index, time, read-only clusters, `transition`, `find`, `output`), never the tracker itself. All
  file writes go through `Run.output` (`output.py`'s `RunOutput`: `path(name)` for whole files,
  buffered `append` for per-frame ones, overwritten on their first flush of a run so runs never
  mix; names may include folders, created on demand, all
  under the run's output directory, `--output-dir`), and each analysis declares its files in `outputs`
  (`OutputFile`), which drives the earlier-run checks (what gets overwritten; files of a family
  like `cls-id<id>.gro` that this run didn't rewrite) and the end-of-run summary. Library users add
  their own analyses with `MolClusters(..., analyses=[...])`, which run after the config-enabled
  built-ins but before the report (`Builtin.last`, placed by `build_builtins`); `Run` dispatches every hook through `Run._hook`, which
  notes which analysis raised an error (the CLI prints those notes), tags what's logged meanwhile
  with the analysis' name (`logger.contextualize(analysis=...)`, shown as a `[Name]` prefix by
  `log._formatter`) and adds up its time in `Run.durations`, logged at the end of the run next to the tracker's.
  Analyses log with loguru's global `logger` (no logging service on `Run`): count per frame,
  summarise once in `finish` (see the `FrameAnalysis` docstring). Built-ins take their options as
  constructor arguments; the registry `analysis/builtins.py`'s `BUILTINS` says how each is built
  from the config and which options it `needs`, and `MolClusters.__init__` builds the ones that
  run (`build_builtins`): those whose options are set, unless the config's `analyses` (class name
  -> bool) turns them off. `MolClsConfig` checks `analyses` against the registry (a lazy import,
  since the analysis package imports the config): unknown names, or one turned on without an
  option it `needs`, are errors. A new built-in is one `BUILTINS` line, in this order:
  `SizeEvolution` (evo.txt), `Lineage` (cluster_events.csv, cluster_lifetimes.csv; adds each
  frame's transition to the report), `SoluteSolvent` (solute_solvent.csv), `ClusterCoordinates`
  (coordinates/cls-n/cls-id/solute-*.gro), `Nucleus` (nucleus_data.csv; its `nuclei` per cluster id are for
  later analyses), `JsonReport` (molclusters.jsonl.zst, see below; `last`).
- The report (`report.py` holds its format, encoders and reader; `analysis/report.py` the
  `JsonReport` that writes it): JSON Lines, a header line then one line per frame, appended
  through `RunOutput` as each frame is analysed, so nothing accumulates in memory and an
  interrupted run leaves its frames readable. Encoded with orjson (NaN becomes null; the stdlib
  wrote invalid bare NaN), every float rounded to `DECIMALS` (6) places where it's computed, and
  compressed per the config's `report_compression` (zstd default, gzip, none): `RunOutput`
  compresses `.zst`/`.gz` names into one stream per run, ending a block on every flush.
  `report.read_report` reads it back frame by frame (compression told from the magic bytes),
  warning when it ends before the header's `n_frames`, or as pandas tables. Its schema is in
  `report.py`'s docstring; changing it bumps `VERSION`. Derived fields (diameter, volume,
  density) are kept on purpose, for users reading the file without the reader. Any analysis adds
  to it by overriding `FrameAnalysis.report_frame`/`report_cluster`: `JsonReport` finds those
  that do through `Run.analyses` (it runs last, so it sees them all), files what they return
  under their class name (two of one name are an error), rounded by `report.rounded`, and lists
  them in the header's `contributors` (`analyses` lists every analysis of the run). `Nucleus` adds its nuclei that way; the report itself
  knows no other analysis.
- What analyses read (`cluster.py`): `MolGroup` is any set of residues (a nucleus, say), with
  properties computed on the group made whole across PBC (`whole()`; never
  moves the shared Universe for good; residues are placed along a spanning tree, a `Cluster`'s
  own graph or a plain group's minimum spanning tree, so groups wider than half the box stay
  whole, and a cluster wrapping around the box is warned about once). The whole positions and
  the geometric properties (`@_per_frame`: `radius`, Rg, sphericity, dipole moment, shape,
  center of mass, dipole, `bsphere`) are computed once per frame and residue set
  (`_frame_cache`), since analyses read them repeatedly (`diameter`/`volume`/`density` all go
  through `radius`). A group can't be changed through what it hands out: every array it
  returns is read-only (`_readonly`; callers change a `.copy()`), `residues` is a new group
  over a copy of its indices each time (a group's `ix` can't be made read-only: MDAnalysis'
  Cython code needs it writable), and its `universe` is a read-only property. Its `radius`
  (and so `diameter`, `volume`, `density`)
  is an equivalent sphere's: sqrt(5/3)·Rg (a uniform sphere's) plus `radius_buffer`, half the
  atoms' mean van der Waals radius by element, for the atoms' size; `radius_of_gyration` is Rg
  itself. `Cluster(MolGroup)` adds the id, a frozen `graph` (`_frozen_graph`: `nx.freeze` plus
  read-only views of the edge, node and graph attribute dicts),
  `birth_time`/`age`, `neighbors` and `distance`. Clusters are read-only outside the tracker:
  its only mutator is `Cluster._update`, which `ClusterTracker.update` calls to swap in the new
  frame's frozen graph. The objects are live, so they describe the current frame only.
- `MolClusters` (`molclusters.py`) is the runner: the constructor checks resnames and builds the
  analyses (the config fills in its own defaults, e.g. `solvent`: the rules' residue names
  other than the solute, since only those can cluster);
  `run()` rewinds the trajectory, creates a fresh `ClusterTracker` (`MolClusters.tracker`, None
  before the first run), calls `tracker.update()` frame by frame and drives the analyses, so
  running again gives the same results (analyses start over in `prepare`). Results are read from
  the analyses, found with `MolClusters.analysis(Type)`; per-frame files are buffered, and the
  whole-run ones are written when the run finishes.
- `ConnectionTable`'s "cm" rule compares each molecule's center of mass made whole by minimum
  image around its first atom (`_whole_centers_of_mass`, no bonds needed): trajectories split
  molecules across PBC, and a split molecule's plain COM can land half a box away.
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
- Units (amu, e, D, Å, Å³, g/cm³, degrees, ps) are pinned in `tests/test_units.py`, each against
  a physical fact rather than a copy of its formula; a new reported quantity or output file gets
  a test there too.
