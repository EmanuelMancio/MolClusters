# Development setup

MolClusters is managed with [uv](https://docs.astral.sh/uv/), which installs the
right Python, the dependencies pinned in `uv.lock` and the package itself in
editable mode.

```bash
git clone https://github.com/EmanuelMancio/MolClusters.git
cd MolClusters
git switch dev
uv sync                       # runtime + dev dependencies (pytest, ruff, ...)
uv run pre-commit install     # pre-commit and commit-msg hooks
uv run molclusters --version
```

Work happens on `dev`, or on branches off it; `main` only receives releases (see
[Commits and releases](releases.md)).

## Checks

Everything CI checks runs locally through pre-commit:

```bash
uv run pre-commit run --all-files
```

| Hook | What it runs |
| --- | --- |
| `ruff-check` | `ruff check --fix` (lint rules in `.ruff.toml`) |
| `ruff-format` | `ruff format` |
| `reuse-lint-file` | [REUSE](https://reuse.software/) license headers |
| `pytest` | the tests affected by your changes (see below) |
| `commitizen` | commit message format (on `commit-msg`) |

All hooks are `local` and run the tools uv installed, so their versions come from
`uv.lock`, which Dependabot bumps weekly (`.github/dependabot.yml`). Don't use
`pre-commit autoupdate`.

The `pytest` hook runs `pytest --testmon --no-cov`:
[pytest-testmon](https://testmon.org/) records which code each test runs (in
`.testmondata`, git-ignored) and reruns only the tests whose code changed since its
last run, so a commit usually takes seconds; the first run, and any change to the
installed packages, runs them all. It doesn't see changes to data files (such as
`tests/data/`), so after changing those, or to check everything, run the whole suite
as CI does, with coverage:

```bash
uv run pytest
```

## Conventions

- **Docstrings** are [NumPy style](https://numpydoc.readthedocs.io/en/latest/format.html),
  on every public module, class and function; the [API reference](../reference/index.md)
  is generated from them.
- **License headers**: every file needs a REUSE/SPDX header, as a comment at the
  top or, for files that can't carry one (data, docs pages), an annotation in
  `REUSE.toml`. Code is `LGPL-3.0-or-later`; configuration and data files are `CC0-1.0`.
- **Line endings**: every text file is stored and checked out with LF
  (`.gitattributes`), so a tool that rewrites a file as CRLF on Windows, as
  `cz bump` does with `pyproject.toml`, still commits it with LF.
- **Logging**: use loguru's `logger`. The package disables its own logger on import
  (`molclusters/__init__.py`); `start_logging` enables it for a run. Detail the
  terminal doesn't need goes through `logger.bind(**FILE_ONLY)`.

## Layout

```text
src/molclusters/
├── main.py          # the `molclusters` command: arguments, loading, errors
├── molclusters.py   # MolClusters, the runner
├── tracker.py       # ClusterTracker and Transition: clusters and their ids
├── conntable.py     # ConnectionTable: which molecules connect, frame by frame
├── cluster.py       # MolGroup and Cluster: groups of molecules and their properties
├── config.py        # MolClsConfig: the configuration, parsed and checked
├── analysis/        # FrameAnalysis, Run, Frame, and the built-in analyses
├── output.py        # OutputFile and RunOutput: the files a run writes
├── report.py        # the JSON report's format, encoders and reader
├── log.py           # logging setup
└── symdict.py       # SymmetricDict, for per-residue-pair rules
tests/               # pytest suite, with small data sets in tests/data
docs/                # this site (mkdocs.yml at the root)
scripts/             # release helpers
```
