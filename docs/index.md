# MolClusters

**MolClusters** follows molecular clusters through a molecular dynamics trajectory:
which molecules are connected to which, how the clusters they form are born, grow,
split, merge and dissolve, how long they live, and what they look like — size,
radius, density, charge, dipole moment and shape — frame by frame.

It reads the simulations [MDAnalysis](https://www.mdanalysis.org/) can (GROMACS,
LAMMPS, CHARMM/NAMD, AMBER, ...), takes the rules that connect molecules
from a small YAML, JSON or TOML file, and writes plain CSV/text tables, per-cluster
coordinate files and a compressed JSON Lines report of every cluster of every frame.

!!! warning "Results from v0.6.0 or earlier?"

    Those versions have errors in dipoles, volumes, densities and, in some cases,
    cluster membership, fixed in v0.7.0. See
    [Known issues in earlier versions](user-guide/known-issues.md) for which results
    are affected and how to correct them.

<div class="grid cards" markdown>

-   :material-rocket-launch:{ .lg .middle } **Get started**

    ---

    Install MolClusters and analyse your first trajectory in a few minutes.

    [:octicons-arrow-right-24: Installation](user-guide/installation.md) ·
    [Quick start](user-guide/quickstart.md)

-   :material-cog:{ .lg .middle } **Configure**

    ---

    Connect molecules by center-of-mass distance or hydrogen bonds, pick solutes,
    nuclei and the analyses to run.

    [:octicons-arrow-right-24: Configuration](user-guide/configuration.md)

-   :material-file-table:{ .lg .middle } **Read the results**

    ---

    What every output file holds, its columns and units, and how to load the JSON
    report with pandas.

    [:octicons-arrow-right-24: Output files](user-guide/outputs.md) ·
    [The JSON report](user-guide/report.md)

-   :material-language-python:{ .lg .middle } **Extend it**

    ---

    Drive MolClusters from Python and add your own per-frame analyses, which can
    write files and add fields to the report.

    [:octicons-arrow-right-24: Python](user-guide/python.md) ·
    [Custom analyses](user-guide/custom-analyses.md)

</div>

## What it does

For every frame of the trajectory, MolClusters

1. **connects molecules** that satisfy one of your *rules*: their centers of mass
   are closer than a cutoff (`cm`), or they share a hydrogen bond (`hb`);
2. **finds the clusters**: groups of two or more molecules linked to each other,
   directly or through other molecules;
3. **keeps each cluster's id** from one frame to the next, so a cluster that grows,
   shrinks, splits or merges is still recognisable as the same cluster (see
   [How clusters are tracked](user-guide/concepts.md));
4. **runs the analyses** on them: cluster sizes over time, births and deaths
   (lineage and lifetimes), solute-solvent clusters, nuclei inside clusters,
   coordinates of the clusters holding solutes, and a full JSON report.

```bash
molclusters traj.xtc topol.tpr config.yml --output-dir results
```

## Authorship and acknowledgments

MolClusters was written by Emanuel Fernandes Dias Mancio[:fontawesome-brands-orcid:](https://orcid.org/0000-0002-0262-711X){ .profile title="ORCID" aria-label="ORCID" }[:academicons-lattes:](http://lattes.cnpq.br/3118069372734394){ .profile title="Lattes CV" aria-label="Lattes CV" },
with the important contribution of Prof. Kaline Coutinho[:fontawesome-brands-orcid:](https://orcid.org/0000-0002-7586-3324){ .profile title="ORCID" aria-label="ORCID" }[:simple-clarivate:](https://www.webofscience.com/wos/author/record/C-2104-2012){ .profile title="Web of Science" aria-label="Web of Science" }[:academicons-lattes:](http://lattes.cnpq.br/9205662588542783){ .profile title="Lattes CV" aria-label="Lattes CV" },
who supervised the work and suggested improvements to it.

We thank the Brazilian funding agencies CAPES, CNPq and FAPESP for the fellowships
and approved research projects. This work was done under FAPESP fellowships for
Emanuel Fernandes Dias Mancio (grants: #2022/01284-1, #2025/15166-9), listed in the
[FAPESP research database](https://bv.fapesp.br/pt/pesquisador/709594/emanuel-fernandes-dias-mancio).

MolClusters is free software, licensed under the
[GNU LGPL v3](https://github.com/EmanuelMancio/MolClusters/blob/main/LICENSE) or any
later version: you can use it as a library from software under any license, and changes
to MolClusters itself that you distribute stay under the LGPL.
