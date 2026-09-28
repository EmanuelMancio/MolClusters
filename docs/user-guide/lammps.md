# LAMMPS systems

MolClusters reads a LAMMPS DATA file as the topology and a LAMMPS dump as the
trajectory. LAMMPS files lack some information the analysis needs, which the
configuration and the dump supply.

```bash
molclusters dump.lammpstrj system.data config.yml
```

Dumps are recognised by the extensions `.lammpsdump`, `.lammpstrj` and `.dump`
(compressed or not). A dump must have the atom `id` and coordinates; MDAnalysis
reads unwrapped (`xu yu zu`), scaled (`xs ys zs`) or plain (`x y z`) coordinates.

## Residue names

A DATA file has a molecule id per atom (the `mol` column of `atom_style full` or
`molecular`), which MDAnalysis exposes as the residue id, but no residue names.
[`lammps_resnames`](configuration.md#lammps_resnames) names them: each name used in
the rest of the configuration maps to the molecule ids it stands for.

```yaml
lammps_resnames:
  SOL: "1-500"          # molecules 1 to 500
  NA: 501
  CL: [502, "510-519"]

rules:
  SOL:
    SOL: cm 3.2
    NA: cm 3.5
```

Every molecule id of the topology must be covered (the error lists the missing
ones), and an id can't be given two names.

Molecule ids must be numbered 1 to N in the DATA file's order, as for any topology
(see [Troubleshooting](troubleshooting.md#residues-numbered-otherwise)).

## Elements and atom names

A DATA file has neither elements nor atom names, and MolClusters needs elements for
the [cluster radius](properties.md#size). Write them to the dump with an `element`
column:

```text
dump        traj all custom 1000 dump.lammpstrj id mol type element x y z
dump_modify traj element C H O N
```

MolClusters reads that column from the dump's first frame and assigns the elements
to the topology's atoms by id. The atom names are then set to the elements, so the
`.gro` coordinate files have meaningful names. Without an `element` column, the run
stops with an error saying so.

## Charges and hydrogen bonds

`hb` rules need partial charges and bonds: use an `atom_style full` DATA file
(charges in the `Atoms` section, and a `Bonds` section), and check that the
hydrogen mass and charge criteria of the [rule](configuration.md#hb-hydrogen-bond)
fit your force field.

## Times

A dump records step numbers, not times, so without more information every time in
the results (the `Time` columns, birth times, lifetimes) is a step number, and the
log warns about it. Give the run's timestep, with its unit, and times are in ps:

```yaml
lammps_timestep: "2 fs"      # units real
# lammps_timestep: "0.001 ps"  # units metal
```

The unit is required, since LAMMPS' own time unit depends on the run's `units`
style. The log's *Trajectory* line shows the times as read, e.g.
`Trajectory: 501 frame(s) from 0 to 1000 ps, every 2 ps`: check it.
