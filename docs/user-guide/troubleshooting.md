# Troubleshooting

When something goes wrong, the terminal shows what; the log file
(`molclusters_<date>_<time>.log` in the output directory) has the full traceback and
the versions of everything involved, which is what a bug report needs.

## Errors

### Residues numbered otherwise

> MolClusters identifies molecules by residue id, and needs the residues numbered 1
> to N in topology order, but residue 1 has resid 0 ...

Molecules are identified by residue id, which must be 1, 2, ..., N in the order of
the topology: no offset, gaps or repeats. GROMACS `.gro`/`.pdb` files of large
systems wrap resids at 99999, and some tools start at 0. Renumber the residues in
the topology, or, from Python, before creating `MolClusters`:

```python
import numpy as np
u.residues.resids = np.arange(1, len(u.residues) + 1)
```

### Invalid configuration

```text
Invalid configuration:
  Invalid rule MAL:MAL ('hb d 3.5 a'): 'hb' flags must be given as '<flag> <number>' pairs.
```

Every problem of the configuration file is listed, with the key it is about. See
[Configuration](configuration.md) for the syntax of each key.

### No van der Waals radius for an element

> No van der Waals radius is known for element(s) ['X'] ...

The cluster radius needs every atom's element. When the topology has no elements,
they are guessed from the letters of the atom names, which fails for names like
`OW1` → `OW` or for virtual sites. Use a topology with elements (e.g. a `.pdb` or
`.psf` with an element column), or rename the atoms.

### The topology has neither elements nor atom names

A LAMMPS DATA file paired with a dump without an `element` column. Add it to the
dump, see [LAMMPS systems](lammps.md#elements-and-atom-names).

### 'lammps_resnames' does not cover LAMMPS molecule id(s) ...

Every molecule id of a LAMMPS topology needs a name; add the listed ids to
[`lammps_resnames`](configuration.md#lammps_resnames).

### Cannot assign donor-hydrogen pairs ... no bond information

`hb` rules find donors through the hydrogens' bonds. Use a topology with bonds, or
switch the rule to `cm`.

## Warnings

**`'rules' names residue(s) [...] that are not in the topology`**
:   A misspelled or absent residue name, which selects nothing. Names are
    case-sensitive.

**`Rule(s) ... never connected any molecules`**
:   A rule whose cutoff is too tight, in the wrong unit (distances are in Å, not
    nm), or, for `hb` rules, whose hydrogen/acceptor criteria match no atoms (zero
    charges, united atoms). See [`hb`](configuration.md#hb-hydrogen-bond).

**`Ignoring unknown config key(s)`**
:   A key MolClusters doesn't know, probably a typo; the warning lists the known
    keys.

**`LAMMPS dumps record step numbers, not times`**
:   Set [`lammps_timestep`](lammps.md#times).

**`The trajectory has no box`**
:   The system is analysed as non-periodic: distances are computed without periodic
    boundary conditions and molecules are taken as they are. Right for a system that
    isn't periodic (e.g. a cluster in vacuum); for a periodic simulation whose file
    lost its box (e.g. an `.xyz`, or a `.pdb` without `CRYST1`), molecules split
    across the boundaries stay split, connections across them are missed and the
    properties of their clusters are wrong: dipole moments by far the most, since
    each atom on the wrong side adds its charge times the box length. Use a
    trajectory with the box.

**`The topology has no partial charges`**
:   A coordinate-only topology (`.gro`, `.pdb`): the charge and dipole moment of
    clusters and nuclei are empty (NaN, `null` in the report), and `hb` rules can't
    run. The other properties don't need charges. Use a topology with charges, such
    as a GROMACS `.tpr`, a CHARMM `.psf` or a LAMMPS DATA file with
    `atom_style full`, to get them.

**`Molecule X N reaches ...% of half the box from its first atom`**
:   `cm` rules, and every property in a topology without bonds, make molecules
    whole by taking each atom at its nearest image to the molecule's first atom,
    which fails for a molecule reaching half the box or more (a long polymer or
    surfactant in a small box): it is torn apart, and its center of mass and the
    properties of the groups holding it are wrong. Warned from 80% of that limit,
    once. Use a larger box, or, for the properties, a topology with bonds.

**`cluster N wraps around the periodic box`**
:   A cluster spanning the box (percolating) is connected to its own periodic image
    and can't be made whole, so its size, density, shape and dipole mean little
    while it does. Reported once per cluster.

**`Solutes were not followed in N frame(s)`**
:   With `follow: [solute]`, a cluster holding more than one solute is not written
    to the `solute-<resid>.gro` files in those frames.

**`N file(s) of an earlier run are left in ...`**
:   The output directory holds files of an earlier run that this one didn't write,
    such as `.gro` files of cluster ids this run doesn't have. Delete them, or use a
    fresh `--output-dir`.

**`distance_backend='OpenMP' was requested, but this MDAnalysis build ... was compiled without OpenMP support`**
:   `distance_backend: OpenMP` has no effect with this MDAnalysis build; see
    [`distance_backend`](configuration.md#distance_backend).

## Performance

- Hydrogen-bond rules are much slower than `cm` rules. When a center-of-mass cutoff
  describes your clusters well enough, prefer it.
- Every cluster's properties are computed every frame for the report and the
  tables; turn off the analyses you don't need with
  [`analyses`](configuration.md#analyses). `ClusterCoordinates` in particular writes
  many files. Writing them on a network file system (an HPC cluster's shared scratch
  space) may go faster with a higher [`flush_threads`](configuration.md#flush_threads),
  or write to a local disk and copy the results afterwards.
- For a first look at a long trajectory, analyse every N-th frame with
  `--traj-memory --in-memory-step N`.
- The end of the log says how long the tracking and each analysis took.

## Reporting a bug

Open an issue at
[github.com/EmanuelMancio/MolClusters/issues](https://github.com/EmanuelMancio/MolClusters/issues)
with the log file of the run (it records the versions, the command and the effective
configuration) and, if you can, a small system that reproduces the problem.
