# Testing

```bash
uv run pytest                      # the whole suite, with an HTML coverage report
uv run pytest tests/test_tracker.py -k split
```

pytest's configuration is in `pytest.ini` only; its `addopts` already include
`--cov=src --cov-report html`, so every run writes a coverage report to `htmlcov/`.
CI (`.github/workflows/tests.yml`) runs the suite on Python 3.12, 3.13 and 3.14 on
every push and pull request to `main` and `dev`.

## The suite

| File | Covers |
| --- | --- |
| `test_config.py` | parsing and validation of every configuration key |
| `test_conntable.py` | `cm` and `hb` connections, PBC, whole-molecule centers of mass |
| `test_tracker.py` | id assignment: formation, growth, split, merge, dissolution, ties, transitions |
| `test_cluster.py` | `MolGroup`/`Cluster`, making groups whole, the per-frame cache |
| `test_properties.py` | properties and outputs of a real system (`tests/data/met-mal`) |
| `test_analysis.py` | `FrameAnalysis`, `Run` and `Frame`, each built-in analysis, errors and logging |
| `test_molclusters.py` | the runner: configuration defaults, runs, user analyses, built-ins |
| `test_output.py` | `RunOutput` buffering, overwriting, compression |
| `test_report.py` | the report format and the reader |
| `test_main.py` | the command line, LAMMPS handling, in-memory loading |
| `test_log.py`, `test_symdict.py` | logging, `SymmetricDict` |
| `test_units.py` | the units of every reported quantity |

## Fixtures and data

`tests/conftest.py` builds **small synthetic systems whose clustering is known by
construction**: `make_universe(frames, n_res, resnames)` returns an in-memory
Universe of two-atom residues, each frame given as the groups of resids that must
form one cluster under a `cm 3.0` rule. Residues of a group sit 2 Å apart in a line,
groups 50 Å apart, in a 2000 Å box, so nothing connects across groups. Prefer these
to real trajectories: the expected result is written in the test.

`tests/data/` holds small real inputs: `met-mal/` (malic acid `MAL` and
2-methylerythritol `MOL`, GROMACS) and `lammps_mini.*` (a LAMMPS DATA file and dump).
In `met-mal.tpr`, MOL has all-zero charges, so only MAL can form hydrogen bonds, and
only at looser angles such as `a 120`.

`captured_logs` collects what the package logs, for tests on messages (the
package's logger is disabled by default and loguru writes to its own stderr
reference, so `capsys` doesn't see it).

## Units

Every unit MolClusters reports (amu, e, D, Å, Å³, g/cm³, degrees, ps) is pinned in
`tests/test_units.py`, **each against a physical fact rather than a copy of its
formula**: a model's literature value, a textbook case, a geometry built to size. A
wrong conversion factor then fails there. A new reported quantity or output file
gets a test there too, and `TestOutputUnits` checks that the files carry the same
units as the properties.
