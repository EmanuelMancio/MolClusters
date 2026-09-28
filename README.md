<!--
SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>

SPDX-License-Identifier: GPL-3.0-only
-->

# Molecular Clusters

Module to analyze formation and life-time of molecular clusters from a molecular dynamics simulation

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

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
See the [developer guide](https://emanuelmancio.github.io/MolClusters/developer-guide/) for the
development setup and conventions.

## License

This script is licensed under the [GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/) license - see [LICENSE](LICENSE) for more details

## Authorship

This package was written by Emanuel Fernandes Dias Mancio ([ORCID](https://orcid.org/0000-0002-0262-711X), [Lattes](http://lattes.cnpq.br/3118069372734394)), with the important contribution of Prof. Kaline Coutinho ([ORCID](https://orcid.org/0000-0002-7586-3324), [Web of Science](https://www.webofscience.com/wos/author/record/C-2104-2012), [Lattes](http://lattes.cnpq.br/9205662588542783)) who supervised the work and gave suggestions to its improvement.

## Acknowledgments

We thank the Brazilian funding agencies CAPES, CNPq and FAPESP for the fellowships and approved research projects. This work was done under FAPESP fellowships for Emanuel Fernandes Dias Mancio (grants: #2022/01284-1, #2025/15166-9), listed in the [FAPESP research database](https://bv.fapesp.br/pt/pesquisador/709594/emanuel-fernandes-dias-mancio).
