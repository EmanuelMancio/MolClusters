# Changing the report format

The JSON report is read by people outside MolClusters, often without its reader, so
its format is versioned and documented in one place.

## Where things live

- `report.py`: the format (`FORMAT`, `VERSION`, `DECIMALS`, `UNITS`,
  `CONNECTION_COLUMNS`), the record encoders (`cluster_record`, `group_record`,
  `connection_rows`, `rounded`, `dumps`) and the reader (`read_report`, `Report`).
  **Its module docstring is the schema.**
- `analysis/report.py`: `JsonReport`, which writes the header in `prepare` and one
  line per frame in `analyse`, appended through `RunOutput`, so nothing accumulates
  in memory and an interrupted run leaves its frames readable.
- The user documentation: [The JSON report](../user-guide/report.md).

## Rules

- **Change `VERSION`** whenever the format changes, and say what changed in the
  comment next to it. The reader refuses reports newer than the version it knows.
- **Update the schema** in `report.py`'s docstring and the user page.
- **Round every float where it is computed** to `DECIMALS` places (`_round` for
  built-in fields; contributions go through `rounded`). NaN is written as `null`
  (orjson does this; the stdlib `json` wrote invalid bare `NaN`).
- **Keep derived fields** (diameter, volume, density) even though they follow from
  others: people read the file without the reader.
- **Units** go in `UNITS`, and get a test in `tests/test_units.py`.

## Adding fields from an analysis

An analysis adds fields by overriding `FrameAnalysis.report_frame` and/or
`report_cluster`. `JsonReport.prepare` finds the analyses that do (through
`Run.analyses`; it runs last, so it sees all of them), files what they return under
their class name, and lists them in the header's `contributors`. Two contributors of
the same class name are an error, since their fields would collide.

A built-in's contribution is part of the format: document it in the schema and
bump `VERSION`. A user analysis' contribution isn't.

## Compression

`report_compression` picks the suffix (`report_name`), and `RunOutput` compresses
any `.zst`/`.gz` name into one stream per run, ending a block on every flush. The
reader tells the compression from the magic bytes, not the name.
