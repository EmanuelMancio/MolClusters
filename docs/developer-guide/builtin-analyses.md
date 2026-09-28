# Adding a built-in analysis

A built-in analysis is a [`FrameAnalysis`](../reference/analysis.md) that the
configuration turns on. Writing the class is the same as for a
[user analysis](../user-guide/custom-analyses.md); making it a built-in takes a line
in the registry.

## 1. Write the class

Put it in its own module under `src/molclusters/analysis/`, named after what it
produces, and export it from `analysis/__init__.py`.

```python
"""Provides `Contacts`, the contacts between solutes and solvents (contacts.csv)."""

from collections.abc import Iterable

from loguru import logger

from ..output import OutputFile
from .base import Frame, FrameAnalysis, Run


class Contacts(FrameAnalysis):
    """Counts the solute-solvent connections of every frame.

    Attributes
    ----------
    counts : list[int]
        The number of such connections, per frame.
    """

    outputs = (OutputFile("contacts.csv"),)

    def __init__(self, solutes: Iterable[str]) -> None:
        self.solutes = set(solutes)
        self.counts: list[int] = []

    def prepare(self, run: Run) -> None:
        self.counts = []  # a run can be repeated: start over here

    def analyse(self, frame: Frame) -> None:
        ...

    def finish(self, run: Run) -> None:
        logger.info(f"{sum(self.counts)} solute-solvent contacts in total")
        ...  # write run.output.path("contacts.csv")
```

Conventions every built-in follows:

- **Options are constructor arguments**, not reads of `run.config`, so the class
  works on its own from Python and in tests.
- **Results are attributes** (arrays, tables, mappings), readable after `run()`
  through `MolClusters.analysis(Type)`; files are written from them.
- **Start over in `prepare`**, not only in `__init__`.
- **Declare every file** in `outputs`; set it in `__init__` if it depends on options.
- **Log a summary in `finish`**, not per frame: count what's worth reporting in
  `analyse`, then log one INFO line or one warning saying how to fix the problem.
  Per-frame messages are DEBUG/TRACE and use loguru's own formatting
  (`logger.debug("cluster {}: ...", cls.id)`), so a filtered message is never built.
  A warning that could repeat is logged once.
- **Read other analyses' results** with `run.analysis(Type)` (they run earlier, so
  on each frame their results are up to date), e.g. `Nucleus.nuclei`.
- **Add to the report** through `report_frame`/`report_cluster` if the analysis has
  per-frame or per-cluster data worth keeping; the report itself knows no analysis.

## 2. Register it

Add a line to `BUILTINS` in `analysis/builtins.py`, in the order it should run:

```python
Builtin(Contacts, lambda c: Contacts(c.solute), needs=("solute",)),
```

- `build` builds it from the configuration;
- `needs` lists the configuration options it can't run without: it runs only when
  they are set, and turning it on explicitly without them is a configuration error;
- `last=True` runs it after the user's analyses too (only the report does).

The configuration then accepts `analyses: {Contacts: false}`, `describe()` lists it,
and `MolClusters` builds it. The current order is `SizeEvolution`, `Lineage`,
`SoluteSolvent`, `ClusterCoordinates`, `Nucleus`, then `JsonReport` (last).

## 3. Test and document it

- Tests in `tests/test_analysis.py`, on the small synthetic systems of
  `tests/conftest.py`, where the expected clusters are known by construction.
- A new quantity or output file gets a units test in `tests/test_units.py` (see
  [Testing](testing.md#units)).
- A section in [Output files](../user-guide/outputs.md), a row in the
  [`analyses` table](../user-guide/configuration.md#analyses), and a line in the
  results table of [Using it from Python](../user-guide/python.md#the-analyses-and-their-results).
- If it adds to the report, the report schema changes: see
  [Changing the report format](report-format.md).
