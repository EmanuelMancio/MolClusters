# Cluster properties

The same properties are computed for clusters and nuclei (any group of molecules),
in the tables and the [JSON report](report.md). All of them are computed on the
group **made whole** across the periodic boundaries: its molecules are placed next
to each other along the cluster's connections (or, for a nucleus, a minimum spanning
tree of its molecules' centers), so a cluster wider than half the box stays whole.
A cluster connected to its own periodic image (e.g. a slab spanning the box) cannot
be made whole; the log warns once per such cluster.

| Property | Report field | Table column | Unit |
| --- | --- | --- | --- |
| number of molecules | `size` | `Size` | — |
| mass | `mass` | — | amu (g/mol) |
| equivalent radius | `radius` | `Radius` | Å |
| radius of gyration | `radius_of_gyration` | — | Å |
| equivalent diameter | `diameter` | — | Å |
| equivalent volume | `volume` | — | Å³ |
| density | `density` | `Density` | g/cm³ |
| total charge | `charge` | `Charge` | e |
| dipole moment | `dipole_moment` | `Dipole` | D (debye) |
| sphericity | `sphericity` | `Spher` | — (0 to 1) |
| shape parameter | `shape_parameter` | `Shape` | — (−0.25 to 2) |

## Mass and charge

The mass is the sum of the atoms' masses, and the charge the sum of their partial
charges, both from the topology.

## Size

A cluster is not a sphere, but its size is summarised by the **equivalent sphere**:
the uniform sphere with the same mass-weighted radius of gyration \(R_g\), whose
radius is \(\sqrt{5/3}\,R_g\). The radius of gyration only sees the atoms' centers,
so the sphere would end at the outer atoms' centers; half the atoms' mean van der
Waals radius \(\bar r_\text{vdW}\) (by element, from MDAnalysis' table) is added to
account for the atoms' own size:

\[
R = \sqrt{\tfrac{5}{3}}\,R_g + \tfrac{1}{2}\,\bar r_\text{vdW},
\qquad
d = 2R,
\qquad
V = \tfrac{4}{3}\pi R^3,
\qquad
\rho = \frac{m}{V}.
\]

Half, not the full van der Waals radius, because bonded atoms' van der Waals spheres
overlap. On methanol-malic acid mixtures the full radius puts a lone malic acid at
0.70 g/cm³ and half at 1.21 g/cm³ (the solid is 1.61 g/cm³), and the choice barely
changes large clusters. The density of small clusters is an estimate; compare it
between clusters rather than with bulk densities.

The report also gives \(R_g\) itself, as `radius_of_gyration`, for analyses that
need it without the equivalent sphere's assumptions.

This needs the atoms' elements: the command line guesses them from the atom names
when the topology has none (e.g. a `.tpr`).

## Dipole moment

The magnitude of the dipole vector \(\boldsymbol\mu = \sum_i q_i(\mathbf r_i -
\mathbf r_\text{COM})\) about the group's center of mass, in debye
(1 D = 0.2081943 e·Å). For a neutral group the dipole does not depend on the
reference point; for a charged group it does, so dipoles of charged clusters are
only comparable with that convention in mind.

The report's `Nucleus.combined_dipole_moment` is the dipole of all a cluster's
nuclei taken together, as one group.

## Shape

Both shape descriptors come from the eigenvalues \(\lambda_1, \lambda_2, \lambda_3\)
of the mass-weighted gyration tensor, with \(\bar\lambda\) their mean, as defined
by Dima and Thirumalai.[^dima2004]

[^dima2004]: R. I. Dima and D. Thirumalai, "Asymmetry in the Shapes of Folded and
    Denatured States of Proteins", *J. Phys. Chem. B* **108**(21), 6564–6570 (2004).
    [doi:10.1021/jp037128y](https://doi.org/10.1021/jp037128y)

**Sphericity** is one minus the asphericity \(\Delta\):

\[
\text{sphericity} = 1 - \Delta,
\qquad
\Delta = \frac{3}{2}\,\frac{\sum_i (\lambda_i - \bar\lambda)^2}{\left(\sum_i \lambda_i\right)^2}.
\]

It is 1 when the three moments are equal (a sphere, but also e.g. a cube), 0 for a
rod (all mass on a line), and 0.75 for a flat disk or ring. It is not Wadell's
surface-area sphericity.

**Shape parameter**:

\[
S = 27\,\frac{\prod_i (\lambda_i - \bar\lambda)}{\left(\sum_i \lambda_i\right)^3},
\]

0 for a sphere, negative for oblate (flattened) shapes down to −0.25 for a disk, and
positive for prolate (elongated) ones up to 2 for a rod.

## Connections

For each connection, the report records the distance between the two molecules:
their centers of mass' distance for `cm` rules, and the donor-acceptor distance of
their shortest hydrogen bond for `hb` rules, with its donor-hydrogen-acceptor angle
and the number of hydrogen bonds between them.
