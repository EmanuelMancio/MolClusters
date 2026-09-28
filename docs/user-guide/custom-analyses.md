# Writing your own analysis

An analysis is a subclass of [`FrameAnalysis`](../reference/analysis.md) that
MolClusters calls on every frame, with the clusters up to date. Pass instances to
`MolClusters` and they run after the built-ins (but before the JSON report, so they
can add to it):

```python
from molclusters import FrameAnalysis, MolClusters, OutputFile


class LargestCluster(FrameAnalysis):
    """Writes the largest cluster of every frame to largest.csv."""

    outputs = (OutputFile("largest.csv"),)  # (1)

    def prepare(self, run):  # (2)
        self.sizes = []
        run.output.append("largest.csv", "Time,Id,Size,Age\n")

    def analyse(self, frame):  # (3)
        if not frame.clusters:
            self.sizes.append(0)
            return
        largest = max(frame.clusters.values(), key=lambda c: c.size)
        self.sizes.append(largest.size)
        frame.output.append(
            "largest.csv", f"{frame.time},{largest.id},{largest.size},{largest.age}\n"
        )

    def report_cluster(self, frame, cluster):  # (4)
        return {"age": cluster.age}


largest = LargestCluster()
mc = MolClusters(u, config, analyses=[largest])
mc.run(output_dir="results")
print(largest.sizes)
```

1.  The files the analysis writes, so the run can report them and tell which files
    of an earlier run are overwritten or left behind.
2.  Called once before the first frame. Start the results over here, not only in
    `__init__`: the same analysis may be run more than once.
3.  Called for every frame, the first one included.
4.  Optional: fields to add to each cluster's record in the JSON report.

## The hooks

A run calls, for every analysis in order:

| Hook | When | Receives |
| --- | --- | --- |
| `prepare(run)` | once, before the first frame | a `Run` |
| `analyse(frame)` | on every frame, the first included (required) | a `Frame` |
| `finish(run)` | once, after the last frame | a `Run` |
| `report_frame(frame)` | on every frame, by the report | a `Frame` |
| `report_cluster(frame, cluster)` | on every cluster of every frame, by the report | a `Frame` and a `Cluster` |

Analyses are called in order, so on a given frame an analysis sees the results of
those before it. The built-ins run first, in the order of the
[`analyses` table](configuration.md#analyses), then yours in the order given, then
the report.

## What an analysis sees

**`Run`**, the run as a whole:

- `universe`: the MDAnalysis Universe;
- `config`: the configuration (a `MolClsConfig`);
- `n_frames`: the number of frames;
- `output`: where to write files (see [below](#writing-files));
- `analysis(Type)`: the first analysis of that type in the run, to use its results
  (while preparing, only the analyses before this one are found).

**`Frame`**, the current frame:

- `index`: 0 for the run's first frame, then 1, 2, ...;
- `time`: its time (ps);
- `universe`: the Universe, at this frame;
- `clusters`: a read-only mapping of cluster id → [`Cluster`](../reference/cluster.md);
- `transition`: how the clusters came from the previous frame's
  ([`Transition`](../reference/tracker.md): `flows`, `born`, `merged`, `dissolved`,
  and `sources(id)`, `destinations(id)`, `absorbed(id)`);
- `find(resid)`: the id of the cluster a molecule is in, or None;
- `output`: the same as `Run.output`.

A **`Cluster`** has an `id`, `birth_time` and `age`, its molecules (`resids`,
`resnames`, `composition`, `size`, `residues`, `atoms`), the `graph` of its
connections (a frozen NetworkX graph whose nodes are resids and whose edges carry
`distance`, and for hydrogen bonds `angle` and `n_hbonds`), `neighbors(resid)`,
and every [property](properties.md) (`radius`, `density`, `dipole_moment`,
`sphericity`, ...). `with cluster.whole() as atoms:` gives its atoms made whole
across the periodic boundaries, restored on exit.

Clusters are read-only. Every array they return (`resids`, `center_of_mass`,
`dipole`, ...) is read-only, so changing one in place raises a `ValueError`: work on
`array.copy()`. The graph's attributes can't be written either; `nx.Graph(cluster.graph)`
gives a copy you can change. `residues` returns a new group each time, and the
Universe is `cluster.universe`.

Clusters are also live: their properties are those of the current frame. To keep
data about a cluster across frames, key it by `cluster.id`, not by the object.

Any group of molecules can get the same properties as a `MolGroup`:
`MolGroup(frame.universe, [resid, ...])`.

## Writing files

Files go through `run.output` (or `frame.output`), always relative to the run's
output directory:

- `output.path("name.csv")` gives the path to write a whole file to, e.g. in
  `finish`;
- `output.append("name.csv", text)` appends to a file frame by frame. Appends are
  buffered (opening a file every frame is slow), and every file starts over in a
  new run. A name ending in `.gz` or `.zst` is compressed as it is written.

Names may include folders (`"mine/per-frame.txt"`), created as needed. Declare every
file in `outputs`; a family of files takes a pattern with an integer placeholder,
e.g. `OutputFile("mine/cluster-<id>.txt")`. Set `outputs` in `__init__` when the
files depend on the analysis' options.

## Adding to the report

Override `report_frame` and/or `report_cluster` to return a mapping of fields (JSON
types, NumPy scalars or arrays). They go under the analysis' class name in each
frame's or cluster's record, their floats rounded like the report's, and the header's
`contributors` lists the analysis. With the example above, every cluster's record
gets `"LargestCluster": {"age": 960.0}`. Return None to add nothing for a given
frame or cluster.

Two analyses of the same class name can't both add to the report.

## Logging

Log with loguru's `logger` (`from loguru import logger`); what an analysis logs is
prefixed with its name. To keep the log readable over long trajectories, count
things per frame and log one summary in `finish`, keep per-frame messages at DEBUG
level, and log a repeated warning only once.

## Errors

An error raised by an analysis stops the run; the message names the analysis, the
hook and the frame (`Raised by LargestCluster.analyse() on frame 12`). Files
appended to so far are kept.
