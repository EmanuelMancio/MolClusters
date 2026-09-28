# Known issues in earlier versions

!!! warning "Results from MolClusters v0.6.0 or earlier"

    Every release up to and including **v0.6.0** has the errors below, all fixed in
    **v0.7.0**. Don't use those versions for new work, and check results you
    obtained with them against this page: some can be corrected by a constant
    factor, others need the analysis run again.

## Summary

| Issue | Affected outputs | Fix for old results |
| --- | --- | --- |
| [Dipole moments 1.89× too small](#dipole-moments) | `Dipole` columns, JSON `Dipole Moment`, `NucleiDipole` | multiply by 1.890 |
| [Volumes 2.152× too small](#volumes-and-densities) | JSON `Volume` | multiply by 2.152 |
| [Densities 0.780× too small](#volumes-and-densities) | `Density` columns, JSON `Density` | multiply by 1.2815 |
| [Cluster properties wrong with nuclei](#cluster-properties-with-nucleus-analysis) | JSON cluster `Radius`, `Diameter`, `Volume`, `Density`, `Dipole Moment`, `Sphericity`, `Shape`, `NucleiDipole` | run again |
| [Clusters wider than half the box torn apart](#clusters-wider-than-half-the-box) | every geometric property of those clusters | run again |
| [Molecules split across the box in "cm" rules](#molecules-split-across-the-box-cm-rules) | which molecules cluster, so every output | run again, if the trajectory has split molecules |
| [Residues not numbered 1 to N](#residues-not-numbered-1-to-n) | which molecules cluster, so every output | renumber and run again |
| [Times with `--traj-memory`](#times-with-traj-memory) | every `Time` column | run again, or correct the times |
| [`.gro` files across runs](#coordinate-files) | `cls-n*.gro`, `cls-id*.gro`, `solute-*.gro` | run again in an empty directory |
| [H-bond connections](#h-bond-connections) | JSON `Connections` | run again, if you use them |
| [Cluster ids in some merges and splits](#cluster-ids) | cluster ids | none needed for sizes and membership |

The factors correct an old value to what that version meant to report. v0.7.0 also
[changed some definitions](#other-changes-in-v070-not-errors), so a corrected old
value and a v0.7.0 value of the same property can still differ.

## Dipole moments

MDAnalysis gives dipoles in e·Å, but they were converted to Debye with the factor for
e·a₀ (atomic units), so every dipole moment was 1/a₀ = 1.890 times too small: an
SPC/E water came out at 1.24 D instead of 2.35 D.

**Affected:** the `Dipole` column of `solute_solvent.csv` and `nucleus_data.csv`, and
`Dipole Moment` and `NucleiDipole` in `molclusters.json`.
**Fix:** multiply by 1.890, except for the
[clusters holding nuclei](#cluster-properties-with-nucleus-analysis), which need a new
run.

## Volumes and densities

The volume was that of a sphere whose radius is the radius of gyration Rg, but a
uniform sphere of radius R has Rg = √(3/5) R, so volumes were (5/3)<sup>3/2</sup> = 2.152
times too small. Densities are mass over that volume, and were also converted from
amu/Å³ to g/cm³ with the inverse of the right factor (×0.6022 instead of ×1.6605);
together, densities came out 0.780 times their intended value.

**Affected:** `Volume` and `Density` in `molclusters.json`, and the `Density` column of
`solute_solvent.csv` and `nucleus_data.csv`. `Radius` and `Diameter` were the radius
of gyration and twice it, as documented at the time, and aren't affected.
**Fix:** multiply volumes by 2.152 and densities by 1.2815 (averages too, since the
factors are constant).

## Cluster properties with nucleus analysis

With `nucleus` configured, analysing a nucleus moved its atoms in place in the shared
coordinates, recentred on the nucleus' first molecule, for the rest of the frame. The
cluster's geometric properties, computed afterwards, used those shifted atoms and could
be far off (a cluster radius of 779 Å instead of 2.3 Å in one test). In the systems
checked, every cluster with two or more nuclei was affected, with dipoles inflated to
around 100 D and `NucleiDipole` collapsed to around 5 D, while clusters with a single
nucleus came out right.

**Affected:** in `molclusters.json`, the geometric properties (`Radius`, `Diameter`,
`Volume`, `Density`, `Dipole Moment`, `Sphericity`, `Shape`) and `NucleiDipole` of
clusters holding nuclei. Cluster membership and ids, `solute_solvent.csv` and
`nucleus_data.csv` were not affected (beyond the unit factors above).
**Fix:** run again.

## Clusters wider than half the box

To make a cluster whole across periodic boundaries, its first molecule was put at the
box centre and every other molecule wrapped into the box. A cluster reaching more than
half a box length from that molecule (a long chain, or any large cluster in a small
box) had its far end placed on the wrong side, and its geometric properties (radius,
volume, density, shape, sphericity, dipole, center of mass) followed.

**Affected:** only clusters that large; smaller clusters are unchanged.
**Fix:** run again if your clusters can span half the box.

## Molecules split across the box ("cm" rules)

"cm" rules compared each molecule's plain center of mass. Trajectories often split
molecules across periodic boundaries (raw GROMACS `.xtc`/`.trr`, wrapped LAMMPS
dumps), and a split molecule's center of mass lands between its pieces, up to half a
box away, so contacts were missed or invented. In one test trajectory, 3.8% of
molecule-frames were split, and the clusters of 29 of 105 frames changed.

**Affected:** runs with "cm" rules on trajectories with split molecules; "hb" rules
aren't affected. Trajectories made whole beforehand (e.g. `gmx trjconv -pbc mol`) are
fine.
**Fix:** run again.

## Residues not numbered 1 to N

Molecules were found by residue number, as the n-th residue for number n. A topology
numbered otherwise (an offset, gaps or repeats: LAMMPS molecule ids with gaps or a
molecule 0, multi-chain PDB files, `.gro` numbers wrapping at 99999) silently put
molecules in the wrong clusters and computed their properties on the wrong atoms.
v0.7.0 refuses such topologies, saying how to renumber them.

**Affected:** only such topologies; those numbered 1 to N are unaffected.
**Fix:** renumber the residues and run again (see
[Troubleshooting](troubleshooting.md#residues-numbered-otherwise)).

## Times with `--traj-memory`

(v0.3.0 to v0.6.0.) Loading the trajectory into memory timed frame i as i × dt from 0.
A trajectory that doesn't start at t = 0 (a continuation run) had all its times
shifted, and one saved at uneven intervals (a LAMMPS dump with changing dump
intervals) had them wrong outright.

**Affected:** the `Time` of every output, only with `--traj-memory` and such a
trajectory.
**Fix:** run again, or shift the times by the trajectory's start time if it was only
offset.

## Coordinate files

The `cls-n<size>.gro`, `solute-<resid>.gro` and (v0.6.0) `cls-id<id>.gro` files were
appended to, so running twice in the same directory added the second run's frames after
the first run's. They also never held the trajectory's first frame.

**Affected:** those files, when a directory held more than one run.
**Fix:** run again in an empty directory.

## H-bond connections

An "hb" rule between a residue name and itself counted H-bonds within a molecule as
connections from that molecule to itself, and a pair of molecules sharing several
H-bonds kept the distance and angle of an arbitrary one.

**Affected:** the `Connections` of `molclusters.json` (entries `(m, m)`, and distances
and angles) and the connection counts in the log. Which molecules cluster was not
affected.
**Fix:** run again if you use the connections.

## Cluster ids

In some merges and splits (ties between equally large pieces, a split into dimers) the
id could go to the wrong piece. v0.7.0 also assigns ids with a new, order-independent
algorithm (see [How clusters are tracked](concepts.md)), so ids can differ from earlier
versions anyway. Cluster membership and sizes are not affected; per-id results
(`cls-id*.gro`, anything following a cluster by id) can be.

## Other changes in v0.7.0 (not errors)

These change results between versions without the old ones being wrong:

- `Radius` and `Diameter` are now those of the equivalent sphere, √(5/3) Rg plus half
  the atoms' mean van der Waals radius, rather than the radius of gyration (reported
  as `radius_of_gyration` in the report); `Volume` and `Density` follow. See
  [Cluster properties](properties.md).
- The JSON report is now `molclusters.jsonl.zst`, one line per frame; see
  [The JSON report](report.md).
- The coordinate files are written in a `coordinates/` folder, and all outputs in the
  directory given by `--output-dir`.

The [changelog](../changelog.md) lists every change.
