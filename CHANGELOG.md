<!--
SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>

SPDX-License-Identifier: LGPL-3.0-or-later
-->

# CHANGELOG

## v0.7.0 (2026-09-28)

### Highlights

- **Errors in earlier results fixed.** v0.6.0 and earlier have errors in dipole
  moments, volumes, densities and, in some cases, cluster membership. See
  [Known issues in earlier versions](https://emanuelmancio.github.io/MolClusters/user-guide/known-issues/)
  for which results are affected and how to correct them.
- **Relicensed** under the GNU LGPL v3 or later (was GPL-3.0-only), so MolClusters can
  be used as a library from software under any license.
- **Published on PyPI** (`pip install molclusters`), with a Zenodo DOI for each release,
  and a [documentation site](https://emanuelmancio.github.io/MolClusters/).

### BREAKING CHANGE

- arrays returned by MolGroup/Cluster are read-only: changing
one in place raises ValueError, so work on `array.copy()`. A Cluster's graph
attributes can't be written (TypeError); `nx.Graph(cluster.graph)` gives a
modifiable copy. `MolGroup.uni` is now the read-only `MolGroup.universe`, and
`residues` returns a new group each call.
- molclusters.json is replaced by molclusters.jsonl.zst
(or .jsonl.gz / .jsonl), in a new format: JSON Lines with snake_case keys
named as the MolGroup properties ("dipole_moment", "shape_parameter", ...),
the composition as a {resname: count} map, connections as
[i, j, distance, angle, n_hbonds] rows, the nuclei under
"Nucleus": {"nuclei", "combined_dipole_moment"}, and floats rounded to 6
decimal places. JsonReport's `data` attribute and encode_* methods are
gone; read the report with `molclusters.report.read_report`.
- Radius, Diameter, Volume and Density in molclusters.json, and
the Radius and Density columns of solute_solvent.csv and nucleus_data.csv,
include the buffer. Reading MolGroup.radius, .diameter, .volume or .density on a
Universe without elements now raises ValueError.
- Radius and Diameter in molclusters.json, and the Radius columns
of solute_solvent.csv and nucleus_data.csv, are now the equivalent sphere's
(sqrt(5/3) = 1.291x the radius of gyration), not the radius of gyration.
- find() returns None instead of False for a molecule
that is in no cluster; check it with `is None`.
- cls-n<size>.gro, cls-id<id>.gro and solute-<resid>.gro
are written to coordinates/ in the output directory instead of the
output directory itself.
- MolClusters' clusters, mol_clt, find(),
clusters_size_evo, size_evolution and nucleus are removed. Use
molcls.tracker.clusters / .mol_clt / .find() after run(), and
molcls.analysis(SizeEvolution).data or molcls.analysis(Nucleus).
- library code using the removed or renamed Cluster and
MDAResidueGroupAnalyzer members must move to the new names (MolGroup,
residues, graph, birth_time, age, distance); constructing a Cluster
now requires a connected group and cluster_id; and two Cluster objects
with the same graph no longer compare equal.
- cluster ids can differ from earlier releases for the same
trajectory and config. A cluster that loses a merge no longer continues in a
remnant elsewhere (the remnant gets a new id), and a pure fragment now wins a
tie against a mixed fragment at any size, not only as a dimer. Outputs keyed by
cluster id (cls-id*.gro, ids in molclusters.json, id-based lifetimes) are
therefore not reproducible against earlier versions.

### Feat

- **cluster**: make everything a group hands out read-only
- **report**: record each group's radius of gyration
- **report**: record the clusters' birth time and lineage in the report
- **analysis**: add Lineage, the clusters' events and lifetimes
- **tracker**: record how each frame's clusters came from the previous frame's
- **analysis**: list the run's analyses in the report header
- **analysis**: let analyses add their own fields to the report
- **analysis**: stream the report as compressed JSON Lines
- **output**: append bytes, and compress .gz and .zst files as they're written
- **analysis**: turn built-in analyses on or off from the config
- **analysis**: tag what each analysis logs with its name and time it
- **config**: give LAMMPS dump times in ps with an optional lammps_timestep
- **cluster**: grow the equivalent sphere by half the atoms' vdW radius
- **cluster**: report the equivalent sphere's radius as Radius and Diameter
- **analysis**: give each frame the run's output, to append to files
- **analysis**: write the cluster coordinates in a coordinates/ folder
- **cli**: write the results to a chosen directory with --output-dir
- **cluster**: give MolGroup and Cluster a read-only public API
- **analysis**: let library users add their own analyses
- **config**: make hb rule hydrogen/acceptor mass and charge criteria configurable
- **cli**: show the version and a description in the help
- **log**: log crashes, third-party warnings, run timing and rule usage
- **log**: log the effective configuration and clearer run progress
- **conntable**: warn when the OpenMP backend can't actually run in parallel
- **conntable**: add OpenMP backend support for cm-rule distances
- **main**: read .dump and .lammpstrj trajectories as LAMMPS dumps

### Fix

- **config**: default solvent to the rules' non-solute residues
- **config**: reject conflicting rules for a pair given both ways
- **output**: start every appended file over on each run
- **log**: shorten Windows warning paths on any OS
- **cli**: keep every frame's time with --traj-memory
- **conntable**: merge the H-bonds a pair of molecules shares into one connection
- **cluster**: keep groups reaching past half the box whole
- **cluster**: convert densities from amu/Å³ to g/cm³ with the right factor
- **cluster**: give the volume of the uniform sphere with the group's Rg
- **conntable**: stop counting H-bonds within a molecule as connections
- **conntable**: take whole molecules' centers of mass for "cm" rules
- **cluster**: convert dipole moments from e·Å, not atomic units, to Debye
- **conntable**: refuse residues not numbered 1 to N instead of mixing them up
- **molcls**: start every run from the first frame with a new tracker
- **molcls**: default the solvent on a copy of the config, not the caller's
- **analysis**: write the clusters' coordinates for the first frame too
- **tracker**: number each tracker's clusters from 1
- **main**: derive atom names from elements for LAMMPS topologies
- **main**: read elements from the LAMMPS dump trajectory
- **molcls**: assign cluster ids in two order-independent steps
- **cluster**: make groups whole without moving the shared positions
- **molcls**: correct older clusters calculation
- **molcls**: correct dominance for dimer split
- **main**: guess elements when the topology does not provide them
- **main**: store the default solvent list as a list, not a set
- **config**: name the offending rule in rule validation errors
- **cluster**: update the ResidueGroup when merging clusters
- **cluster**: rebuild the ResidueGroup of clusters split off by separate()
- replace __dict__ in MolClusters.__slots__ with the real slot list
- **release**: keep uv.lock in sync and use GitHub's exact skip-ci token

### Refactor

- **tracker**: rename _gen_origin_cluster_counter to _count_origins
- **tracker**: return None from find() for a molecule in no cluster
- **molcls**: make MolClusters a runner, found analyses by type
- **molcls**: delete the unused radius_evolution and cluster index writer
- **analysis**: write the JSON report as a FrameAnalysis
- **analysis**: find the nuclei as a FrameAnalysis
- **analysis**: write the cluster coordinates as a FrameAnalysis
- **analysis**: run the solute-solvent analysis as a FrameAnalysis
- **analysis**: give analyses a run context, a base class and an output service
- **analysis**: run the cluster size evolution as a FrameAnalysis
- **molcls**: extract cluster tracking into ClusterTracker
- **output**: write cluster .gro frames via an in-memory stream
- explicit unhashability, guard private mda API, isolate tests

### Perf

- **cluster**: compute scalar properties once per frame
- **output**: buffer cluster .gro frames instead of reopening files per frame
- **conntable**: vectorize connection building and restrict hb searches

## v0.6.0 (2026-09-22)

### Feat

- **config**: warn when lammps_resnames ranges overlap for the same name
- **config**: add lammps_resnames mapping for LAMMPS topologies
- also write per-cluster-identity coordinate output

### Fix

- **ci**: pin setup-uv to an exact tag, v10 doesn't exist
- **pre-commit**: scope reuse-lint-file hook to the pre-commit stage
- temp-file hygiene, logging consistency, and TypeVar bound
- **config**: store lammps_resnames as ranges instead of per-id entries
- correct dominance-adjacent correctness bugs from code review
- correct config drift and untrack test fixtures gitignore rule
- restrict reuse-lint-file hook to pre-commit stage

### Refactor

- **config**: simplify rule and solute-keyword parsing

## v0.5.0 (2026-09-22)

- build system change: `poetry` -> `uv`

## v0.4.0 (2026-05-09)

### BREAKING CHANGE

- dropped python 3.11 support.

### Feat

- **config**: implement MolClsConfig class to handle configuration
- **log**: allow string to be passed to logger initializer
- **log**: change log decorator so no parentheses are needed
- **log**: add logging tools and configuration
- drop python 3.11 and add new dependencies

### Fix

- correct package name in version

### Refactor

- update type hints to use built-in types (e.g., list, dict) across multiple files

## v0.3.1 (2026-03-03)

### Fix

- update config keys from 'solute' to 'nucleus' and adjust JSON indentation
- use upper element name for vdwradii

## v0.3.0 (2025-04-11)

### Feat

- **main**: add trajectory loading options for memory management

### Fix

- **main**: correct element extraction from atom name

## v0.2.0 (2025-04-11)

### Feat

- **MolClusters**: add support for ignore_composition in cluster analysis
- **SubConnTable**: add resnames property

## v0.1.2 (2025-04-11)

### Fix

- **ConnectionTable**: show correct values in attributes

### Refactor

- **ConnectionTable**: stop holding HB info after usage

## v0.1.1 (2025-04-11)

### Fix

- resolve numpy RuntimeWarning
- correctly encodes connection information
- correctly open temporary file in MolClusters
- update import path for vdwradii
- resolve type hints errors

## v0.1.0 (2025-04-11)

### Feat

- add attributes to connections

## v0.1.0-dev.0 (2025-04-11)

### Feat

- add radii information to topology
- nucleus analysis
- add clusters nucleus analysis
- add resids properties to Cluster
- add software version and absolute paths to json output
- change dipole to debye
- add json export for molclusters
- add version option to command
- add properties analysis of clusters
- implements initial HB rule algorithm
- initial solute-solvent analysis algorithm
- use YAML file as settings for analysis
- add MolCluster Class and dominance algorithm
- add init files so test can import MolClusters
- finish implementation of ConnTable
- change code to support ConnTable as class
- add generator of sub graphs of ConnTable
- change conn_tab to class ConnTable
- add cluster size evolution analysis
- read large chunks of trajectory at once
- add size property to Cluster
- add trajectory analysis

### Fix

- add encode_connections method to MolClustersData for connection encoding
- avoid error when input does not have solute information
- solve clusters of one molecule
- remove unused variables
- unify use of coord instead of trajectory
- add guards to eliminate errors in evolution calculations
- correct MolClusters to use instance universe
- free command line interface
- import Counter
- consider the reciprocal in input file
- ConnTab raises error if access key do not exist
- change handling of zero distance
- save evolution data to file
- correct position of data in clusters_evo_info
- correct merging of multiple clusters_index
- correctly update conn_tab during trajectory
- correct identification of new clusters
- solve dict size change runtime error
- solve merge and separation problem
- correct identification of clusters

### Refactor

- avoid unnecessary nucleus checks
- move cluster ag and uni to holder class
- break up update_cluster in smaller functions
- modify Cluster to use SubConnTables
- add SubConnTable to ConnTable
- changes ConnTab construction
