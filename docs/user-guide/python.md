# Using it from Python

Everything the command line does is available from Python, e.g. in a notebook, and
the results can be read straight from the analyses, as NumPy arrays and pandas
tables, without going through the files.

```python
import MDAnalysis as mda

from molclusters import MolClusters, start_logging
from molclusters.analysis import Lineage, Nucleus, SizeEvolution
from molclusters.config import MolClsConfig, read_config

start_logging("INFO", filename="results/molclusters.log")  # optional (1)

u = mda.Universe("topol.tpr", "traj.xtc")
u.guess_TopologyAttrs(to_guess=["elements"])  # a .tpr has no elements (2)

config = MolClsConfig(
    rules={"MAL": {"MAL": "hb d 3.5 a 120", "MOL": "cm 6.0"}},
    solute=["MAL"],
    nucleus=["solute"],
    analyses={"ClusterCoordinates": False},
)
# or: config = read_config("config.yml")

mc = MolClusters(u, config)
mc.run(output_dir="results")

evo = mc.analysis(SizeEvolution).data          # Time, NClusters, Min/Avg/MaxSize per frame
lifetimes = mc.analysis(Lineage).lifetimes()    # the cluster_lifetimes.csv table
```

1.  MolClusters logs nothing unless logging is started. `start_logging` sets up
    the same terminal output and log files as the command line; call it once.
2.  The command line does this, and the other preparation steps below, for you.

## What the command line does besides

The `molclusters` command prepares the Universe before running the analysis. From
Python, do what your system needs of it:

- **Elements**: needed for the cluster radius. Guess them from the atom names when
  the topology has none: `u.guess_TopologyAttrs(to_guess=["elements"])`.
- **LAMMPS**: the command line applies [`lammps_resnames`](lammps.md#residue-names),
  reads the dump's `element` column and passes
  [`lammps_timestep`](lammps.md#times) to the reader. From Python, give the residue
  names and elements to the Universe yourself (`u.add_TopologyAttr("resnames", ...)`,
  `u.add_TopologyAttr("elements", ...)`) and the timestep as the reader's `dt`, in ps:
  `mda.Universe("system.data", "dump.lammpstrj", format="LAMMPSDUMP", dt=0.002)`.
- **Frames in memory**: `--traj-memory` loads the frames into memory keeping each
  frame's time. MDAnalysis' own `u.transfer_to_memory()` renumbers the times as
  0, dt, 2dt, ..., which loses a trajectory's start time; slice the trajectory
  only if that doesn't matter to you.

The residues must be numbered 1 to N in topology order; `MolClusters` checks it
and says how to renumber them otherwise.

## Configuration objects

`MolClsConfig` takes the same keys as the [configuration file](configuration.md),
as keyword arguments, and checks them the same way (a `pydantic.ValidationError`
lists every problem). `read_config(path)` reads a file into one.

`config.describe()` returns the *Effective configuration* text of the log.

## The analyses and their results

`mc.analyses` lists the analyses of the run, in the order they run. Find one by
class with `mc.analysis(Type)` (None when it doesn't run). Their results, once
`run()` returns:

| Analysis | Results |
| --- | --- |
| `SizeEvolution` | `data`: array, one row per frame: time, number of clusters, min, average and max size |
| `Lineage` | `lifetimes()`: the lifetimes table (pandas); `lives`: id → what is known of each cluster |
| `SoluteSolvent` | `data`: array, one row per frame, the columns of `solute_solvent.csv` |
| `Nucleus` | `data`: array, the columns of `nucleus_data.csv`; `nuclei`: id → the last frame's nuclei |
| `ClusterCoordinates` | `solute_ids`: the solutes' resids; `follow_skipped`: frames not followed, by cluster id |
| `JsonReport` | `name`: the report's file name, to [read it back](report.md#reading-it) |

`run()` can be called again: every run starts from the first frame with fresh ids
and results. `mc.tracker` is the [`ClusterTracker`](../reference/tracker.md) of the
last run, left at the last frame.

## Tracking only

To follow the clusters without any analysis, e.g. to inspect them frame by frame
yourself, use the tracker directly:

```python
from molclusters import ClusterTracker

u.trajectory[0]
tracker = ClusterTracker(u, config)       # the clusters of the current frame
for ts in u.trajectory[1:]:
    tracker.update()                      # reconcile with the new frame
    for cid, cluster in tracker.clusters.items():
        print(ts.time, cid, cluster.size, cluster.radius, cluster.age)
    print(tracker.transition.merged)      # {absorbed id: id it merged into}
```

Clusters ([`Cluster`](../reference/cluster.md)) are live objects: their properties
are those of the frame the Universe is on. See the
[API reference](../reference/index.md) for everything they offer.

To run analyses of your own on the clusters, rather than a loop like this one, see
[Writing your own analysis](custom-analyses.md): they get the same bookkeeping,
files and report as the built-ins.
