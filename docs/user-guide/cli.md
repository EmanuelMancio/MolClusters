# Command line

```bash
molclusters TRAJ TOP INP [options]
```

```text
usage: molclusters [-h] [--traj-memory] [--in-memory-step N]
                   [--log-level {TRACE,DEBUG,INFO,WARNING,ERROR}]
                   [--output-dir OUTPUT_DIR] [--version]
                   traj top inp
```

## Arguments

`traj`
:   The trajectory, in any format MDAnalysis reads (told from its extension). Files
    ending in `.lammpsdump`, `.lammpstrj` or `.dump` (compressed or not) are read as
    LAMMPS dumps.

`top`
:   The topology, with residue names, masses, partial charges and bonds (e.g. a
    GROMACS `.tpr`). See [Quick start](quickstart.md#1-what-you-need).

`inp`
:   The [configuration file](configuration.md) (YAML, JSON or TOML).

!!! note "Argument order"
    The trajectory comes first, then the topology.

## Options

`--output-dir DIR`
:   Directory for the results and the log, created if missing. Default: the current
    directory. Files of an earlier run in it are overwritten (the log says which);
    files of an earlier run that this one didn't write, e.g. `cls-id<id>.gro` files
    for ids this run doesn't have, are left and warned about.

`--traj-memory`
:   Load the whole trajectory into memory before the analysis. It can make the
    analysis faster, notably for compressed or slow-to-read formats, but needs
    memory for every frame's coordinates (about 12 bytes per atom per frame).

`--in-memory-step N`
:   Keep only every N-th frame when loading the trajectory into memory (requires
    `--traj-memory`); the whole analysis then runs on those frames, with their
    original times. Default: 1 (every frame). Useful for a quick first look at a long
    trajectory; remember that events between the kept frames are not seen (see
    [How clusters are tracked](concepts.md#events)).

`--log-level LEVEL`
:   Minimum level of the messages logged: `TRACE`, `DEBUG`, `INFO` (default),
    `WARNING` or `ERROR`. `DEBUG` logs every lineage event and other per-frame detail.

`--version`
:   Print the version and exit.

`-h`, `--help`
:   Print the help and exit.

## The log

Every run writes `molclusters_<YYYYMMDD>_<HHMM>.log` to the output directory, plus a
machine-readable copy with a `.json` suffix (one JSON record per message). Besides
what the terminal shows, the log file has the module and line of each message, the
progress every 10% of the frames and, when a run fails, the full traceback.

The first lines record the MolClusters, Python and library versions, the command
line, the system and the effective configuration, so a log is enough to know how a
result was produced.

Messages from an analysis are prefixed by its name, e.g. `[Lineage]`.

## Exit status

| Code | Meaning |
| --- | --- |
| 0 | the analysis finished |
| 1 | an error stopped it; the message says what went wrong (and which analysis, if one raised it), the log file has the traceback |
| 130 | interrupted with ++ctrl+c++ |

An interrupted or failed run keeps what it had written so far: the per-frame files
(the JSON report, `cluster_events.csv`, the `.gro` files) hold the frames analysed
before it stopped. The tables written at the end (`evo.txt`, `cluster_lifetimes.csv`,
`solute_solvent.csv`, `nucleus_data.csv`) are not written.
