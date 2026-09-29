# Configuration

The analysis settings are read from a YAML (`.yml`, `.yaml`), JSON (`.json`) or TOML
(`.toml`) file, told apart by its extension. Only `rules` is required.

=== "YAML"

    ```yaml
    rules:
      MAL:
        MAL: hb d 3.5 a 120
        MOL: cm 6.0
      MOL:
        MOL: cm 6.0
    solute: [MAL]
    nucleus: [solute]
    follow: [solute]
    report_compression: gzip
    analyses:
      ClusterCoordinates: false
    ```

=== "TOML"

    ```toml
    solute = ["MAL"]
    nucleus = ["solute"]
    follow = ["solute"]
    report_compression = "gzip"

    [rules.MAL]
    MAL = "hb d 3.5 a 120"
    MOL = "cm 6.0"

    [rules.MOL]
    MOL = "cm 6.0"

    [analyses]
    ClusterCoordinates = false
    ```

=== "JSON"

    ```json
    {
      "rules": {
        "MAL": {"MAL": "hb d 3.5 a 120", "MOL": "cm 6.0"},
        "MOL": {"MOL": "cm 6.0"}
      },
      "solute": ["MAL"],
      "nucleus": ["solute"],
      "follow": ["solute"],
      "report_compression": "gzip",
      "analyses": {"ClusterCoordinates": false}
    }
    ```

Residue names are case-sensitive and must match the topology's. A key MolClusters
doesn't know is ignored with a warning (check it for typos), and a residue name the
topology doesn't have is warned about, since it would silently select nothing.

| Key | Default | Summary |
| --- | --- | --- |
| [`rules`](#rules) | *required* | which residue pairs connect, and how |
| [`solute`](#solute) | none | residue names of the solutes |
| [`solvent`](#solvent) | the rules' other residue names | residue names of the solvents |
| [`nucleus`](#nucleus) | none | residue names that make up nuclei |
| [`follow`](#follow) | none | write one coordinate file per solute molecule |
| [`ignore_composition`](#ignore_composition) | none | cluster compositions not to count as clusters |
| [`analyses`](#analyses) | all that can run | turn built-in analyses on or off |
| [`report_compression`](#report_compression) | `zstd` | compression of the JSON report |
| [`distance_backend`](#distance_backend) | `serial` | MDAnalysis backend for `cm` distances |
| [`flush_threads`](#flush_threads) | `4` | output files written at once |
| [`lammps_resnames`](#lammps_resnames) | none | residue names of LAMMPS molecule ids |
| [`lammps_timestep`](#lammps_timestep) | none | timestep of a LAMMPS run, for dump times |

## `rules`

A nested mapping `first residue name → second residue name → rule`. Each rule says
when two molecules of those residue names are connected in a frame. Residue-name
pairs without a rule never connect, so a rule is needed for every pair that can be
part of the same cluster, including a species with itself.

Rules are symmetric: `MAL: {MOL: ...}` and `MOL: {MAL: ...}` are the same rule.
Write each pair once. If both are written, they must be the same rule: different
rules for the same pair are an error.

### `cm`: center-of-mass distance

```yaml
MOL: cm 6.0
```

Connects two molecules whose centers of mass are at most the given distance apart,
in ångströms. Each molecule is made whole across the periodic boundaries before its
center of mass is computed, and the distance follows the minimum-image convention.
Finding these connections needs only masses; the clusters' charge and dipole moment
need partial charges as well (see [Quick start](quickstart.md#1-what-you-need)).

### `hb`: hydrogen bond

```yaml
MAL: hb d 3.5 a 120
```

Connects two molecules that share at least one hydrogen bond, as found by
MDAnalysis' [`HydrogenBondAnalysis`](https://docs.mdanalysis.org/stable/documentation_pages/analysis/hydrogenbonds.html).
The rule takes any of these `flag value` pairs, in any order, each optional:

| Flag | Default | Meaning |
| --- | --- | --- |
| `d` | 3.5 | donor-acceptor distance cutoff (Å) |
| `a` | 150 | donor-hydrogen-acceptor angle cutoff (degrees): a straight bond is 180° |
| `hmin` | 0.9 | hydrogens have a mass (amu) above this... |
| `hmax` | 1.1 | ...and below this |
| `hq` | 0.3 | hydrogens have a partial charge (e) above this |
| `aq` | -0.5 | acceptors have a partial charge (e) below this |

So `hb` alone is `hb d 3.5 a 150 hmin 0.9 hmax 1.1 hq 0.3 aq -0.5`, MDAnalysis'
own defaults for guessing hydrogens and acceptors. Donors are the atoms the
hydrogens are bonded to. This rule therefore needs **masses, partial charges and
bonds** in the topology (a GROMACS `.tpr` or a CHARMM `.psf` has them). Hydrogen
bonds within a molecule are ignored.

!!! tip
    United-atom or coarse-grained models, or topologies whose charges are all zero
    for a species, find no hydrogens or acceptors with the defaults. Check the
    *Connections per frame, by rule* line at the end of the log: a rule that never
    connected anything is warned about. Loosen `a` (e.g. `a 120`) if bonds you
    expect are missed at the default 150°.

A pair of molecules joined by several hydrogen bonds (e.g. a carboxylic acid dimer)
is one connection; the report keeps the geometry of the shortest bond and the
number of bonds.

## `solute`

A list of residue names: the species of interest. Setting it turns on the
[`SoluteSolvent`](outputs.md#solute_solventcsv) and
[`ClusterCoordinates`](outputs.md#coordinates) analyses, and defines the `solute`
keyword: in `nucleus`, `follow` and `ignore_composition`, the word `solute` stands
for every residue name listed here.

`solute` itself cannot be used in `rules` or inside `solute`.

## `solvent`

A list of residue names: the solvents, for [`SoluteSolvent`](outputs.md#solute_solventcsv).
When `solute` is set and `solvent` isn't, every residue name used in
[`rules`](#rules) that isn't a solute is a solvent. Residues no rule names can't
join a cluster, so they aren't solvents here (ions, say). The log's effective
configuration lists the solvents.

## `nucleus`

A list of residue names, or `solute`. A **nucleus** is a connected group of
molecules of these residue names inside a cluster: the cluster's graph restricted to
them, split into its connected pieces. A cluster can hold several nuclei. Setting it
turns on the [`Nucleus`](outputs.md#nucleus_datacsv) analysis, which also adds each
cluster's nuclei to the [report](report.md).

## `follow`

```yaml
follow: [solute]
```

Writes, for each solute molecule, a `coordinates/solute-<resid>.gro` file with the
cluster it belongs to in every frame where that cluster holds no other solute. Only
the `solute` keyword has an effect so far; other entries are warned about and
ignored.

## `ignore_composition`

A list of compositions, each a list of residue names (or `solute`). A connected
group whose set of residue names is exactly one of these is not treated as a
cluster: it gets no id and appears in no output, and its molecules count as free.
Only which residue names are present matters, not how many molecules of each.

```yaml
ignore_composition:
  - [MOL]         # groups of MOL only
  - [MOL, SOL]    # groups made of MOL and SOL, and nothing else
```

## `analyses`

Turns built-in analyses on or off by class name. An analysis left out runs when
the options it needs are set:

| Analysis | Needs | Writes |
| --- | --- | --- |
| `SizeEvolution` | — | [`evo.txt`](outputs.md#evotxt) |
| `Lineage` | — | [`cluster_events.csv`, `cluster_lifetimes.csv`](outputs.md#cluster_eventscsv) |
| `SoluteSolvent` | `solute` | [`solute_solvent.csv`](outputs.md#solute_solventcsv) |
| `ClusterCoordinates` | `solute` | [`coordinates/*.gro`](outputs.md#coordinates) |
| `Nucleus` | `nucleus` | [`nucleus_data.csv`](outputs.md#nucleus_datacsv) |
| `JsonReport` | — | [`molclusters.jsonl.zst`](report.md) |

```yaml
analyses:
  ClusterCoordinates: false   # skip the .gro files, which can be many
  JsonReport: false
```

An unknown name, or an analysis turned on (`true`) without an option it needs, is
an error.

## `report_compression`

How the [JSON report](report.md) is compressed: `zstd` (the default, fastest and
smallest, `molclusters.jsonl.zst`), `gzip` (readable by more tools,
`molclusters.jsonl.gz`) or `none` (`molclusters.jsonl`).

## `distance_backend`

`serial` (default) or `OpenMP`: the MDAnalysis acceleration backend for the `cm`
rules' distance searches. It often makes no difference: MDAnalysis ignores it when
it picks its grid search (typical when the cutoff is much smaller than the box), and
some MDAnalysis builds, such as the PyPI wheels for Windows, have OpenMP compiled
out (MolClusters warns when it detects this). `hb` rules always run serially.

## `flush_threads`

How many output files are written at once (default `4`, at least `1`). The files are
kept in memory and written when the buffer fills up and at the end of the run, where
a run with one coordinate file per cluster can have thousands of them. Writing
several at once overlaps the wait for the file system on each. On a local disk the
gain is small: on Windows, 4 is about as fast as it gets. On a network file system,
such as an HPC cluster's shared scratch space, where every file costs a round trip
to a server, more may help; `1` writes the files one after another. It changes only
how long writing takes, not what's written.

## `lammps_resnames`

LAMMPS topologies have molecule ids but no residue names. This maps each name used
elsewhere in the file to the LAMMPS molecule id(s) it stands for: a single id, an
inclusive range `"first-last"`, or a list of them.

```yaml
lammps_resnames:
  SOL: "1-500"
  NA: 501
  CL: [502, "510-519"]
```

Every molecule id of the topology must be covered. See [LAMMPS systems](lammps.md).

## `lammps_timestep`

The timestep of the LAMMPS run, with its unit (`fs`, `ps` or `ns`), e.g. `"2 fs"`
for `units real` or `"0.001 ps"` for `units metal`. LAMMPS dumps store step numbers
only; with this, their times are in ps. See [LAMMPS systems](lammps.md#times).

## Environment variables

Any key can also be set by an environment variable named `MOLCLS_` followed by the
key in upper case, e.g. `MOLCLS_REPORT_COMPRESSION=gzip`; list and mapping values
are written as JSON. The configuration file takes precedence over the environment.
