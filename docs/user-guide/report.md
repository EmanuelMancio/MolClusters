# The JSON report

The report, `molclusters.jsonl.zst`, records every cluster of every frame: its
molecules, composition, properties and connections, plus whatever other analyses add
(the nuclei of each cluster, how the clusters of each frame came from the previous
frame's). It is the most complete output, and the one to use for any analysis the
tables don't cover.

## Reading it

From Python, `molclusters.report.read_report` reads it back, whatever its
compression:

```python
from molclusters.report import read_report

report = read_report("results/molclusters.jsonl.zst")
report.header["units"]            # {'time': 'ps', 'mass': 'amu', ...}
report.header["config"]           # the configuration the run used

clusters = report.clusters()      # one row per cluster per frame (pandas)
connections = report.connections()  # one row per connection per frame
```

`clusters()` gives the columns `frame`, `time`, `id`, `birth_time`, `size`, `resids`,
`mass`, `volume`, `radius`, `radius_of_gyration`, `diameter`, `density`, `charge`,
`dipole_moment`,
`sphericity`, `shape_parameter`, one `composition.<resname>` column per residue name
(NaN when absent), and one `<Analysis>.<field>` column per field another analysis
added, e.g. `Nucleus.nuclei`.

`connections()` gives `frame`, `time`, `cluster` (its id), `i`, `j` (the two
resids), `distance` (Å), and for hydrogen bonds `angle` (degrees) and `n_hbonds`
(`<NA>` for `cm` rules).

For reports too large for memory, go through the frames one at a time; each is a
`dict` (see [the format](#format)):

```python
for frame in report:          # or report.frames()
    largest = max(frame["clusters"], key=lambda c: c["size"], default=None)
    if largest:
        print(frame["time"], largest["id"], largest["size"])
```

A report cut short (the run was interrupted or killed) is read up to its last
complete frame, with a warning.

Without MolClusters, any JSON Lines reader works: decompress the file
(`zstd -d molclusters.jsonl.zst`, or `gunzip` for a `.gz` report) and read one JSON
object per line, keeping in mind that the first line is the [header](#header), not a
frame. `pandas.read_json(path, lines=True)` also reads the compressed file directly,
with the header as its first row.

## Format

The report is [JSON Lines](https://jsonlines.org/): one JSON object per line,
written frame by frame, so it never has to fit in memory and an interrupted run
leaves every frame it analysed readable. It is compressed as the
[`report_compression`](configuration.md#report_compression) setting says: zstd
(`.jsonl.zst`, default), gzip (`.jsonl.gz`) or none (`.jsonl`).

Floats are rounded to 6 decimal places. A value that doesn't exist (NaN) is `null`.

### Header

The first line describes the run:

```json
{
  "format": "molclusters-report",
  "version": 4,
  "software": "MolClusters 0.6.0",
  "trajectory": "/abs/path/traj.trr",
  "topology": "/abs/path/topol.tpr",
  "n_frames": 53,
  "decimals": 6,
  "units": {"time": "ps", "birth_time": "ps", "mass": "amu", "volume": "angstrom^3",
            "radius": "angstrom", "radius_of_gyration": "angstrom",
            "diameter": "angstrom", "density": "g/cm^3", "charge": "e",
            "dipole_moment": "D", "distance": "angstrom", "angle": "degrees"},
  "connection_columns": ["i", "j", "distance", "angle", "n_hbonds"],
  "analyses": ["SizeEvolution", "Lineage", "SoluteSolvent", "ClusterCoordinates",
               "Nucleus", "JsonReport"],
  "contributors": ["Lineage", "Nucleus"],
  "config": {"rules": {"MAL": {"MAL": "hb d 3.5 a 120", "MOL": "cm 6.0"}}, "...": "..."}
}
```

(shown on several lines here; in the file it is a single line). `trajectory` and
`topology` are `null` for a trajectory built in memory. `analyses` lists every
analysis of the run, and `contributors` those that added fields to the records.

### Frames

Every other line is a frame:

```json
{"frame": 1, "time": 960.0, "clusters": [...], "Lineage": {...}}
```

`frame` is the trajectory frame number, `time` its time in ps.

### Clusters

Each cluster of a frame is:

```json
{
  "id": 3,
  "birth_time": 120.0,
  "size": 3,
  "resids": [107, 169, 171],
  "composition": {"MAL": 3},
  "mass": 402.26099,
  "volume": 834.55576,
  "radius": 5.840575,
  "radius_of_gyration": 3.981872,
  "diameter": 11.681149,
  "density": 0.80039,
  "charge": 0.0,
  "dipole_moment": 15.696198,
  "sphericity": 0.73869,
  "shape_parameter": 0.138796,
  "connections": [[169, 107, 2.684737, 171.434846, 1], [171, 169, 3.12, 158.2, 2]],
  "Nucleus": {"nuclei": [...], "combined_dipole_moment": 15.696198}
}
```

`birth_time` is the time of the frame the cluster first appeared in (the first
frame's for clusters already there). Each connection is a row of the header's
`connection_columns`: the two resids, the distance (center-of-mass distance for `cm`
rules, donor-acceptor distance for `hb` rules), and for hydrogen bonds the D-H-A angle
and the number of hydrogen bonds between the two molecules, of which the angle and
distance are the shortest one's (`null` for `cm` rules). The properties are described
in [Cluster properties](properties.md).

### What other analyses add

Analyses add their own fields under their class name, in a frame's or a cluster's
record:

`Nucleus` (per cluster)
:   `nuclei`: one record per nucleus of the cluster, with the cluster's fields but
    `id`, `birth_time` and `connections`; `combined_dipole_moment`: the dipole
    moment of all its nuclei together (`null` without nuclei).

`Lineage` (per frame)
:   how the frame's clusters came from the previous frame's:

    ```json
    "Lineage": {"flows": [[0, 9, 3], [1, 0, 1], [1, 1, 2], [2, 0, 1], [2, 11, 1],
                          [3, 3, 2], [6, 3, 2], ...],
                "born": [9, 10, 11, ...], "merged": [[6, 3]],
                "dissolved": [2, 5, 7, 8]}
    ```

    Each flow is `[previous id, current id, molecules]`, id 0 standing for no
    cluster (free molecules); a cluster that carries on has a flow to itself.
    `born` are the clusters new in this frame, `merged` pairs
    `[absorbed id, id it merged into]`, and `dissolved` the clusters that ended
    without merging. On the first frame every cluster is born.

    Here, three free molecules formed cluster 9; cluster 1 kept two molecules and
    lost one; cluster 2 dissolved, one molecule going free and one into the new
    cluster 11; and cluster 6 merged into cluster 3, bringing two molecules.

[Your own analyses](custom-analyses.md#adding-to-the-report) can add fields the same
way.

### Versions

The `version` field changes whenever the format does:

| Version | Change |
| --- | --- |
| 4 | added `radius_of_gyration` to clusters and nuclei |
| 3 | added `birth_time` and the `Lineage` fields |
| 2 | the first JSON Lines report |

A report of an earlier version simply lacks the later fields; `read_report` reads
every version from 2 on. MolClusters 0.6
and earlier wrote `molclusters.json`, a single JSON document, which `read_report`
does not read: load it with Python's `json` module.
