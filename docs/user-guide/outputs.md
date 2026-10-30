# Output files

Every file goes to the output directory (`--output-dir`, or the `output_dir` of
`MolClusters.run` from Python). Which files a run writes depends on the analyses the
[configuration](configuration.md#analyses) runs. Units are ps, Å, amu, e, D and
g/cm³ throughout (see [Cluster properties](properties.md)).

| File | Analysis | Written |
| --- | --- | --- |
| [`evo.txt`](#evotxt) | `SizeEvolution` | at the end |
| [`cluster_events.csv`](#cluster_eventscsv) | `Lineage` | frame by frame |
| [`cluster_lifetimes.csv`](#cluster_lifetimescsv) | `Lineage` | at the end |
| [`solute_solvent.csv`](#solute_solventcsv) | `SoluteSolvent` | at the end |
| [`coordinates/`](#coordinates) | `ClusterCoordinates` | frame by frame |
| [`nucleus_data.csv`](#nucleus_datacsv) | `Nucleus` | at the end |
| [`molclusters.jsonl.zst`](report.md) | `JsonReport` | frame by frame |
| [`molclusters_<date>_<time>.log`](cli.md#the-log) | — | as it runs |

## `evo.txt`

The number of clusters and their sizes (in molecules), one row per frame, as a
whitespace-separated table with a `#` header line (`numpy.loadtxt` reads it
directly):

```text
# Time NClusters MinSize AvgSize MaxSize
0.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00 ...
9.600000000000000000e+02 8.000000000000000000e+00 2.000000000000000000e+00 ...
```

| Column | Meaning |
| --- | --- |
| `Time` | time of the frame (ps) |
| `NClusters` | number of clusters |
| `MinSize`, `AvgSize`, `MaxSize` | smallest, average and largest cluster size; 0 in a frame without clusters |

## `cluster_events.csv`

Every event of the run, one row each (see [Events](concepts.md#events)):

```text
Frame,Time,Event,Cluster,Other,NMols
1,48.0,formation,1,,2
2,96.0,dissolution,1,,2
8,384.0,split,22,18,1
44,2112.0,merge,296,300,2
```

Read these rows as:

- `formation,1,,2`: cluster 1 formed of 2 free molecules;
- `dissolution,1,,2`: cluster 1, which held 2 molecules, dissolved;
- `split,22,18,1`: the new cluster 22 split off cluster 18, taking 1 of its
  molecules, while cluster 18 carried on;
- `merge,296,300,2`: cluster 300 merged into cluster 296, bringing 2 molecules.

In `split` and `merge` rows, the molecules move from `Other` to `Cluster`: from the
cluster that split to the new one, and from the absorbed cluster to the one that
absorbed it.

| Column | Meaning |
| --- | --- |
| `Frame` | trajectory frame number of the first frame showing the event |
| `Time` | its time (ps) |
| `Event` | `formation`, `split`, `merge` or `dissolution` |
| `Cluster` | the cluster the event is about (see below) |
| `Other` | the other cluster involved, if any |
| `NMols` | the number of molecules concerned |

| `Event` | `Cluster` | `Other` | `NMols` |
| --- | --- | --- | --- |
| `formation` | the new cluster, made only of free molecules | — | its size |
| `split` | the new cluster | the cluster it took molecules from | molecules taken from `Other` |
| `merge` | the cluster that absorbed `Other` | the cluster that ended | molecules of `Other` it received |
| `dissolution` | the cluster that ended without merging | — | its last size |

A cluster made of pieces of several clusters has a `split` row for each of them, and
a merge of several clusters has a `merge` row for each one absorbed. The clusters of
the first frame have no event: whatever formed them happened before the run.

A new cluster with even a single molecule from an earlier cluster counts as a
`split` of it. For instance, when a dimer breaks up and one of its molecules pairs
with a free molecule, the new pair gets a `split` row with `NMols` 1, and the dimer
a `dissolution` row, in the same frame. Such rows are common. To keep only the
splits where a piece of two or more molecules broke off, filter on `NMols >= 2`.

## `cluster_lifetimes.csv`

One row per cluster seen in the run:

| Column | Meaning |
| --- | --- |
| `Id` | the cluster id |
| `BirthTime` | time of the first frame it was in (ps) |
| `DeathTime` | time of the first frame it was gone from (ps); empty if alive at the end |
| `Lifetime` | `DeathTime − BirthTime`, or up to the last frame if alive at the end (ps) |
| `NFrames` | number of frames it was in |
| `BornAtStart` | it was there at the first frame: its lifetime is a lower bound |
| `AliveAtEnd` | it was still there at the last frame: its lifetime is a lower bound |
| `Origin` | `initial` (there at the first frame), `formation` or `split` |
| `Parents` | the clusters it took molecules from at birth, `;`-separated, most molecules first |
| `Fate` | `merge`, `dissolution` or `alive` |
| `MergedInto` | the cluster that absorbed it, for `Fate` = `merge` |
| `BirthSize`, `MaxSize`, `LastSize` | its size at birth, largest and last size |

For a lifetime distribution free of censoring, keep the rows where both
`BornAtStart` and `AliveAtEnd` are false:

```python
import pandas as pd

lives = pd.read_csv("cluster_lifetimes.csv")
whole = lives[~lives.BornAtStart & ~lives.AliveAtEnd]
print(whole.Lifetime.describe())
```

Their distribution, one bar per multiple of the frame spacing (see
[How precise the times are](concepts.md#how-precise-the-times-are)), on a log scale
since most clusters are short-lived:

```python
import matplotlib.pyplot as plt

counts = whole.Lifetime.value_counts().sort_index()
dt = counts.index.min()  # the frame spacing, the shortest lifetime
plt.bar(counts.index / 1000, counts, width=0.8 * dt / 1000)
plt.yscale("log")
plt.xlabel("Lifetime (ns)")
plt.show()
```

![Bar chart of the number of clusters by lifetime, on a log scale, for the quick start's run analysed frame by frame (48 ps apart). Of the 13,380 clusters born and ended within the run, 9,905 lived a single frame; the counts fall steeply, to about ten per bar around 1 ns, and the longest-lived cluster lasted 3.5 ns.](../assets/outputs/lifetimes-light.svg#only-light)
![Bar chart of the number of clusters by lifetime, on a log scale, for the quick start's run analysed frame by frame (48 ps apart). Of the 13,380 clusters born and ended within the run, 9,905 lived a single frame; the counts fall steeply, to about ten per bar around 1 ns, and the longest-lived cluster lasted 3.5 ns.](../assets/outputs/lifetimes-dark.svg#only-dark)

In the quick start's system, analysed frame by frame (48 ps apart), three clusters in
four lived a single frame.

## `solute_solvent.csv`

The clusters that hold at least one solute and at least one solvent molecule,
averaged over those clusters, one row per frame:

| Column | Meaning |
| --- | --- |
| `Time` | time of the frame (ps) |
| `NCls` | number of such clusters |
| `NSolt`, `NSolv` | their average number of solute and solvent molecules (0 without such clusters) |
| `Radius` | their average [radius](properties.md#size) (Å) |
| `Density` | average density (g/cm³) |
| `Charge` | average total charge (e) |
| `Dipole` | average dipole moment (D) |
| `Spher` | average [sphericity](properties.md#shape) |
| `Shape` | average [shape parameter](properties.md#shape) |

The averages are empty (NaN) in frames without such clusters.

For instance, the solute and solvent molecules of these clusters over time, smoothed
over 1 ns:

```python
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("results/solute_solvent.csv")
mean = data.rolling(21, center=True).mean()  # 21 frames: 1 ns at 48 ps per frame
plt.plot(data.Time / 1000, mean.NSolt, label="solute")
plt.plot(data.Time / 1000, mean.NSolv, label="solvent")
plt.xlabel("Time (ns)")
plt.legend()
plt.show()
```

![Two lines over 50 ns, the average number of malic acid (solute) and of 2-methylerythritol (solvent) molecules in the clusters holding both, as a 1 ns running mean. Both rise over the first 20 ns, malic acid from about 1 to 4 and 2-methylerythritol from about 1 to 3, peak around 22 ns at about 6.4 and 4.7, and from 25 ns on hold between about 4 and 6 and between 2 and 4. Malic acid is above 2-methylerythritol almost all the time.](../assets/outputs/solute-solvent-light.svg#only-light)
![Two lines over 50 ns, the average number of malic acid (solute) and of 2-methylerythritol (solvent) molecules in the clusters holding both, as a 1 ns running mean. Both rise over the first 20 ns, malic acid from about 1 to 4 and 2-methylerythritol from about 1 to 3, peak around 22 ns at about 6.4 and 4.7, and from 25 ns on hold between about 4 and 6 and between 2 and 4. Malic acid is above 2-methylerythritol almost all the time.](../assets/outputs/solute-solvent-dark.svg#only-dark)

## `nucleus_data.csv`

The [nuclei](configuration.md#nucleus) inside the clusters, one row per frame:

| Column | Meaning |
| --- | --- |
| `Time` | time of the frame (ps) |
| `NNuc` | the average number of nuclei of the clusters that have any (0 without nuclei) |
| `Size` | the nuclei's average size (molecules) |
| `Radius`, `Density`, `Charge`, `Dipole`, `Spher`, `Shape` | the nuclei's average properties, as in `solute_solvent.csv` |

The averages are over every nucleus of the frame, and empty (NaN) in frames without
nuclei. Each cluster's own nuclei are in the [report](report.md).

## `coordinates/`

The clusters holding at least one solute, made whole across the periodic boundaries,
as GROMACS `.gro` files (one frame per time the cluster was seen; the title line
says which cluster and when):

| File | Holds |
| --- | --- |
| `cls-n<size>.gro` | every such cluster of that size, whichever its id: an ensemble of what an N-mer looks like |
| `cls-id<id>.gro` | one cluster's own trajectory, whatever its size |
| `solute-<resid>.gro` | with [`follow: [solute]`](configuration.md#follow), the cluster of that solute molecule, in the frames where it holds no other solute |

```text
Cluster-5 - Time = 144.0
   30
  181MAL    CM0    1  16.158  16.040  16.224
  ...
```

These files can be many: one per cluster id holding a solute. Turn the analysis off
with `analyses: {ClusterCoordinates: false}` when you don't need them.

!!! note "Atom counts change between frames"
    A cluster's size changes over its life, so the frames of a `cls-id<id>.gro` or
    `solute-<resid>.gro` file don't all have the same number of atoms. Most
    trajectory readers expect a constant atom count; read such files frame by
    frame, or split them.

A `cls-id<id>.gro` file follows one cluster through its life. Here is cluster 2872
of the quick start's run, analysed frame by frame, at four of its 72 frames: it split
off another cluster with 4 molecules, grew to 17 within 0.6 ns, and was down to 4
malic acids when it merged into another cluster 3.4 ns after its birth.

![One cluster drawn at four times, side by side at the same scale, malic acid in blue and 2-methylerythritol in grey. At 11.33 ns it holds 4 molecules, two of each, just split off another cluster; at 11.66 ns, 13 molecules; at 11.90 ns, 17, its largest; at 14.74 ns, its last frame, 4 malic acids.](../assets/outputs/cls-id-light.svg#only-light)
![One cluster drawn at four times, side by side at the same scale, malic acid in blue and 2-methylerythritol in grey. At 11.33 ns it holds 4 molecules, two of each, just split off another cluster; at 11.66 ns, 13 molecules; at 11.90 ns, 17, its largest; at 14.74 ns, its last frame, 4 malic acids.](../assets/outputs/cls-id-dark.svg#only-dark)

The picture was drawn with VMD from four of the file's frames, split as the note
above says.

## Re-running in the same directory

A new run overwrites the files it writes. Families of files (the `.gro` files) need
not have the same members from one run to the next, so a file of an earlier run
that this run didn't write is left in place, and the log warns about it: delete
those to keep only the latest run's results, or use a fresh `--output-dir`.
