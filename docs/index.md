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

## Citing

If you use MolClusters in your research, please cite it. Every release is archived on
Zenodo with its own DOI: cite the version you used, listed on the
[Zenodo record](https://doi.org/10.5281/zenodo.23023320), or
[10.5281/zenodo.23023320](https://doi.org/10.5281/zenodo.23023320) for MolClusters as a
whole (it resolves to the latest version). The repository's
[CITATION.cff](https://github.com/EmanuelMancio/MolClusters/blob/main/CITATION.cff) has
the details, and GitHub's "Cite this repository" button formats them.

### Software and methods it builds on

If you can, please also cite the work MolClusters relies on:

- **MDAnalysis**, which reads the trajectories and computes the distances,
  hydrogen bonds and group properties. MDAnalysis asks for both of its papers:
    - N. Michaud-Agrawal, E. J. Denning, T. B. Woolf and O. Beckstein, "MDAnalysis:
      A Toolkit for the Analysis of Molecular Dynamics Simulations", *J. Comput.
      Chem.* **32**, 2319–2327 (2011).
      [doi:10.1002/jcc.21787](https://doi.org/10.1002/jcc.21787)
    - R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
      J. Domański, D. L. Dotson, S. Buchoux, I. M. Kenney and O. Beckstein,
      "MDAnalysis: A Python Package for the Rapid Analysis of Molecular Dynamics
      Simulations", *Proc. 15th Python in Science Conf.*, 98–105 (2016).
      [doi:10.25080/Majora-629e541a-00e](https://doi.org/10.25080/Majora-629e541a-00e)
- **With `hb` rules**, MDAnalysis' hydrogen bond analysis: P. Smith, R. M. Ziolek,
  E. Gazzarrini, D. M. Owen and C. D. Lorenz, "On the interaction of hyaluronic acid
  with synovial fluid lipid membranes", *Phys. Chem. Chem. Phys.* **21**, 9845–9857
  (2019). [doi:10.1039/C9CP01532A](https://doi.org/10.1039/C9CP01532A)
- **For sphericity or the shape parameter**: R. I. Dima and D. Thirumalai, "Asymmetry
  in the Shapes of Folded and Denatured States of Proteins", *J. Phys. Chem. B*
  **108**, 6564–6570 (2004). [doi:10.1021/jp037128y](https://doi.org/10.1021/jp037128y)
- **For radius, diameter, volume or density**, the sources of the van der Waals radii
  (MDAnalysis' table), for the elements in your system: A. Bondi, *J. Phys. Chem.*
  **68**, 441–451 (1964), [doi:10.1021/j100785a001](https://doi.org/10.1021/j100785a001);
  R. S. Rowland and R. Taylor, *J. Phys. Chem.* **100**, 7384–7391 (1996),
  [doi:10.1021/jp953141+](https://doi.org/10.1021/jp953141%2B); M. Mantina,
  A. C. Chamberlin, R. Valero, C. J. Cramer and D. G. Truhlar, *J. Phys. Chem. A*
  **113**, 5806–5812 (2009), [doi:10.1021/jp8111556](https://doi.org/10.1021/jp8111556).
- **NetworkX**, for the connection graphs and clusters: A. A. Hagberg, D. A. Schult and
  P. J. Swart, "Exploring Network Structure, Dynamics, and Function using NetworkX",
  *Proc. 7th Python in Science Conf.*, 11–15 (2008).
  [doi:10.25080/TCWV9851](https://doi.org/10.25080/TCWV9851)
- **NumPy**: C. R. Harris *et al.*, "Array programming with NumPy", *Nature* **585**,
  357–362 (2020). [doi:10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2)

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
