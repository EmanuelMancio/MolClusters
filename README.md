<!--
SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>

SPDX-License-Identifier: LGPL-3.0-or-later
-->

# MolClusters

[![PyPI](https://img.shields.io/pypi/v/molclusters.svg)](https://pypi.org/project/molclusters/)
[![Python versions](https://img.shields.io/pypi/pyversions/molclusters.svg)](https://pypi.org/project/molclusters/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23023320.svg)](https://doi.org/10.5281/zenodo.23023320)
[![Tests](https://github.com/EmanuelMancio/MolClusters/actions/workflows/tests.yml/badge.svg?branch=dev)](https://github.com/EmanuelMancio/MolClusters/actions/workflows/tests.yml)
[![Coverage](https://codecov.io/gh/EmanuelMancio/MolClusters/branch/dev/graph/badge.svg)](https://codecov.io/gh/EmanuelMancio/MolClusters)
[![Docs](https://github.com/EmanuelMancio/MolClusters/actions/workflows/docs.yml/badge.svg)](https://emanuelmancio.github.io/MolClusters/)
[![License: LGPL v3+](https://img.shields.io/badge/license-LGPL--3.0--or--later-blue.svg)](https://github.com/EmanuelMancio/MolClusters/blob/main/LICENSE)
[![REUSE compliant](https://img.shields.io/badge/REUSE-compliant-green.svg)](https://reuse.software/)

Track the formation, growth, splitting, merging and lifetime of molecular clusters in
molecular dynamics trajectories (GROMACS, LAMMPS and any other format
[MDAnalysis](https://www.mdanalysis.org/) reads).

**Documentation: <https://emanuelmancio.github.io/MolClusters/>**

> ⚠️ **Results from v0.6.0 or earlier?** Those versions have errors in dipoles, volumes,
> densities and, in some cases, cluster membership, fixed in v0.7.0. See
> [Known issues in earlier versions](https://emanuelmancio.github.io/MolClusters/latest/user-guide/known-issues/)
> for which results are affected and how to correct them.

## Installation

MolClusters needs Python 3.12 or newer, and is published on
[PyPI](https://pypi.org/project/molclusters/):

```bash
uv tool install molclusters   # the molclusters command, in an environment of its own
# or
pip install molclusters       # into your current environment, e.g. to use it from Python
```

The latest development version installs from this repository with
`pip install git+https://github.com/EmanuelMancio/MolClusters@dev`.

It depends on [MDAnalysis](https://www.mdanalysis.org/), [NetworkX](https://networkx.org/),
[NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/),
[Pydantic](https://docs.pydantic.dev/) and a few others, installed automatically.

## Usage

```bash
molclusters traj.xtc topol.tpr config.yml --output-dir results
```

The [user guide](https://emanuelmancio.github.io/MolClusters/latest/user-guide/) covers the
configuration file, the output files and using MolClusters from Python; `molclusters -h`
lists the command-line options.

## Contributing

Bug reports and pull requests are welcome; for major changes, please open an issue first to
discuss what you would like to change. See
[CONTRIBUTING.md](https://github.com/EmanuelMancio/MolClusters/blob/main/CONTRIBUTING.md) and the
[developer guide](https://emanuelmancio.github.io/MolClusters/latest/developer-guide/) for the
development setup and conventions.

## Citing

If you use MolClusters in your research, please cite it. Every release is archived on
Zenodo with its own DOI: cite the version you used, listed on the
[Zenodo record](https://doi.org/10.5281/zenodo.23023320), or
[10.5281/zenodo.23023320](https://doi.org/10.5281/zenodo.23023320) for MolClusters as a
whole (it resolves to the latest version).
[CITATION.cff](https://github.com/EmanuelMancio/MolClusters/blob/main/CITATION.cff) has the
details, and GitHub's "Cite this repository" button formats them.

If you can, please also cite the work MolClusters builds on:
[MDAnalysis](https://www.mdanalysis.org/citations/) (Michaud-Agrawal *et al.*,
[J. Comput. Chem. 2011](https://doi.org/10.1002/jcc.21787); Gowers *et al.*,
[SciPy 2016](https://doi.org/10.25080/Majora-629e541a-00e)), and, depending on the
features you use, the method papers listed in the
[documentation](https://emanuelmancio.github.io/MolClusters/latest/#software-and-methods-it-builds-on).

## License

MolClusters is free software, licensed under the
[GNU LGPL v3](https://choosealicense.com/licenses/lgpl-3.0/) or any later version; see
[LICENSE](https://github.com/EmanuelMancio/MolClusters/blob/main/LICENSE) and, for the GNU
GPL v3 it builds on, [COPYING](https://github.com/EmanuelMancio/MolClusters/blob/main/COPYING). You can use it
as a library from software under any license; changes to MolClusters itself that you
distribute must stay under the LGPL. Every file's license is declared per the
[REUSE](https://reuse.software/) specification.

## Authorship

This package was written by Emanuel Fernandes Dias Mancio ([ORCID](https://orcid.org/0000-0002-0262-711X), [Lattes](http://lattes.cnpq.br/3118069372734394)), with the important contribution of Prof. Kaline Coutinho ([ORCID](https://orcid.org/0000-0002-7586-3324), [Web of Science](https://www.webofscience.com/wos/author/record/C-2104-2012), [Lattes](http://lattes.cnpq.br/9205662588542783)) who supervised the work and gave suggestions to its improvement.

## Acknowledgments

We thank the Brazilian funding agencies CAPES, CNPq and FAPESP for the fellowships and approved research projects. This work was done under FAPESP fellowships for Emanuel Fernandes Dias Mancio (grants: #2022/01284-1, #2025/15166-9), listed in the [FAPESP research database](https://bv.fapesp.br/pt/pesquisador/709594/emanuel-fernandes-dias-mancio).
