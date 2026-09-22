<!--
SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>

SPDX-License-Identifier: GPL-3.0-only
-->

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

MolClusters analyzes the formation and lifetime of molecular clusters from a molecular dynamics
simulation trajectory, built on MDAnalysis and NetworkX. It ships as a CLI (`molclusters`) whose
entry point is `src/molclusters/main.py:main`.

## Commands

Dependencies are managed with `uv` (see `uv.lock`, `uv_build` backend).

- Install/sync deps: `uv sync`
- Run the CLI: `uv run molclusters <trajectory> <topology> <input.yaml> [--traj-memory] [--in-memory-step N]`
- Run tests: `uv run pytest`
- Run a single test: `uv run pytest tests/test_cluster.py::TestCluster::test_cluster_creation`
- Lint: `uv run ruff check .`
- Format: `uv run ruff format .`
- Run all pre-commit hooks (ruff check --fix, ruff format, conventional-commit message lint, REUSE license lint): `uv run pre-commit run --all-files`
- Version bump / changelog (commitizen, semver, tag format `v$version`): `uv run cz bump`

Note: pytest config lives solely in `pytest.ini` (`addopts` already includes `--cov=src
--cov-report html`), so a plain `uv run pytest` produces an HTML coverage report under `htmlcov/`.

## Architecture

Core pipeline, all under `src/molclusters/`:

- `main.py` — CLI entry point (argparse). Builds an `MDAnalysis.Universe` from the topology +
  trajectory, applies `lammps_resnames` (if set) via `_apply_lammps_resnames` since LAMMPS
  topologies carry no residue names, fills in missing vdW radii per atom, reads the analysis
  config, then drives `MolClusters.run()`.
- `config.py` — `MolClsConfig` (pydantic-settings, env prefix `MOLCLS_`), parses the YAML/JSON/TOML
  input file. Per-residue-pair `rules` (`"cm <dist>"` for a center-of-mass distance cutoff, or
  `"hb [d <dist>] [a <angle>]"` for hydrogen bonding) are validated and compiled into a
  `SymmetricDict[(resname, resname), Rule]`. Also owns `solute`/`solvent`/`nucleus`/`follow`/
  `ignore_composition`, including expanding the `"solute"` keyword into the configured solute
  resnames, and `lammps_resnames` (name -> LAMMPS molecule id or `"first-last"` id range),
  compiled into a resid -> resname lookup (`resname_for_resid`) for LAMMPS inputs, which have no
  residue names of their own.
- `conntable.py` — `ConnectionTable` builds a per-frame NetworkX graph of molecule-molecule
  connections from the compiled rules, using MDAnalysis capped-distance search for `cm` rules and
  `HydrogenBondAnalysis` for `hb` rules. `subconntables()` yields each connected component
  (`_SubConnTable`) as a candidate cluster for the current frame.
- `cluster.py` — `MDAResidueGroupAnalyzer` wraps an MDAnalysis `ResidueGroup` with derived physical
  properties (radius of gyration, sphericity, dipole moment, density, shape parameter, ...),
  including periodic-image "make whole" centering. `Cluster` extends it with a NetworkX graph of
  intra-cluster connectivity and add/remove/merge/separate operations.
- `molclusters.py` — `MolClusters` orchestrates the per-frame loop. Each frame it rebuilds the
  connectivity table, then reconciles the new connected components against the previous frame's
  clusters with a dominance algorithm (`__construct_dominance` / `__check_dominance`) that decides
  whether a component continues an existing cluster, is a merge, or is a new cluster — this is what
  keeps cluster identity stable across formation/merge/split events. It optionally runs solute-solvent
  and "nucleus" sub-cluster analyses. `MolClustersData` serializes per-frame cluster state to
  `molclusters.json`; `run()` also writes `evo.txt`, `solute_solvent.csv`, `nucleus_data.csv`, and
  per-cluster `.gro`/`.ndx` files.
- `symdict.py` — `SymmetricDict`, an order-independent `(a, b) == (b, a)` mapping used for pairwise
  rules and hydrogen-bond analyzers.
- `log.py` — loguru-based logging setup (`start_logging`) plus a `_logger_wraps` decorator for
  function entry/exit/exception tracing.

The analysis is frame-oriented: `MolClusters.run()` walks the trajectory frame by frame, mutating
in-memory cluster state, and only writes results to disk once the full run finishes.

## Conventions

- Requires Python >=3.12.
- Ruff (`.ruff.toml`) enforces numpy-style docstrings, double quotes, 88-char lines, and the
  bandit (`S`), annotations (`ANN`), pathlib (`PTH`), and perf (`PERF`) rule sets; files under
  `tests/` are exempt from `S101`, `ANN201`, and docstring rules.
- Commit messages must follow Conventional Commits (enforced by pre-commit + commitizen); commitizen
  uses them to drive semver bumps and `CHANGELOG.md` generation.
- Every file needs a REUSE/SPDX license header (see `REUSE.toml`, `LICENSES/`), enforced by
  pre-commit's `reuse-lint-file` hook.
