# User guide

This guide is for people analysing their simulations with MolClusters, from the
command line or from Python.

- **[Installation](installation.md)**: install the `molclusters` command.
- **[Quick start](quickstart.md)**: a first analysis, start to finish.
- **[How clusters are tracked](concepts.md)**: what a connection, a cluster and a
  cluster id are, and what the formation, split, merge and dissolution events mean.
  Worth reading before interpreting any result.
- **[Configuration](configuration.md)**: every key of the configuration file.
- **[Command line](cli.md)**: the `molclusters` command's arguments and options.
- **[Output files](outputs.md)**: what each output file holds, column by column.
- **[The JSON report](report.md)**: the full per-frame, per-cluster record, and how
  to read it back.
- **[Cluster properties](properties.md)**: how radius, density, dipole moment,
  sphericity and the other properties are computed, and their units.
- **[LAMMPS systems](lammps.md)**: residue names, elements and times for LAMMPS
  topologies and dumps.
- **[Using it from Python](python.md)**: run an analysis from a script or notebook
  and get the results as arrays and tables.
- **[Writing your own analysis](custom-analyses.md)**: add per-frame analyses of
  your own.
- **[Troubleshooting](troubleshooting.md)**: common errors and warnings, and what
  to do about them.
