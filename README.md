<!--
SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>

SPDX-License-Identifier: GPL-3.0-only
-->

# MolClusters

[![Tests](https://github.com/EmanuelMancio/MolClusters/actions/workflows/tests.yml/badge.svg?branch=dev)](https://github.com/EmanuelMancio/MolClusters/actions/workflows/tests.yml)
[![Docs](https://github.com/EmanuelMancio/MolClusters/actions/workflows/docs.yml/badge.svg)](https://emanuelmancio.github.io/MolClusters/)
[![License: GPL v3](https://img.shields.io/badge/license-GPL--3.0--only-blue.svg)](LICENSE)
[![REUSE compliant](https://img.shields.io/badge/REUSE-compliant-green.svg)](https://reuse.software/)

Track the formation, growth, splitting, merging and lifetime of molecular clusters in
molecular dynamics trajectories (GROMACS, LAMMPS and any other format
[MDAnalysis](https://www.mdanalysis.org/) reads).

**Documentation: <https://emanuelmancio.github.io/MolClusters/>**

## Installation

MolClusters needs Python 3.12 or newer, and is installed from this repository:

```bash
uv tool install git+https://github.com/EmanuelMancio/MolClusters
# or
pip install git+https://github.com/EmanuelMancio/MolClusters
```

It depends on [MDAnalysis](https://www.mdanalysis.org/), [NetworkX](https://networkx.org/),
[NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/),
[Pydantic](https://docs.pydantic.dev/) and a few others, installed automatically.

## Usage

```bash
molclusters traj.xtc topol.tpr config.yml --output-dir results
```

The [user guide](https://emanuelmancio.github.io/MolClusters/user-guide/) covers the
configuration file, the output files and using MolClusters from Python; `molclusters -h`
lists the command-line options.

## Contributing

Bug reports and pull requests are welcome; for major changes, please open an issue first to
discuss what you would like to change. See [CONTRIBUTING.md](CONTRIBUTING.md) and the
[developer guide](https://emanuelmancio.github.io/MolClusters/developer-guide/) for the
development setup and conventions.

## Citing

If you use MolClusters in your research, please cite it: [CITATION.cff](CITATION.cff) has the
details, and GitHub's "Cite this repository" button formats them.

## License

MolClusters is free software, licensed under the
[GNU GPL v3](https://choosealicense.com/licenses/gpl-3.0/) only; see [LICENSE](LICENSE).
Every file's license is declared per the [REUSE](https://reuse.software/) specification.

## Authorship

This package was written by Emanuel Fernandes Dias Mancio ([ORCID](https://orcid.org/0000-0002-0262-711X), [Lattes](http://lattes.cnpq.br/3118069372734394)), with the important contribution of Prof. Kaline Coutinho ([ORCID](https://orcid.org/0000-0002-7586-3324), [Web of Science](https://www.webofscience.com/wos/author/record/C-2104-2012), [Lattes](http://lattes.cnpq.br/9205662588542783)) who supervised the work and gave suggestions to its improvement.

## Acknowledgments

We thank the Brazilian funding agencies CAPES, CNPq and FAPESP for the fellowships and approved research projects. This work was done under FAPESP fellowships for Emanuel Fernandes Dias Mancio (grants: #2022/01284-1, #2025/15166-9), listed in the [FAPESP research database](https://bv.fapesp.br/pt/pesquisador/709594/emanuel-fernandes-dias-mancio).
