# Installation

MolClusters needs **Python 3.12 or newer**. It is installed from its GitHub
repository (it isn't published on PyPI).

!!! danger "Don't use cluster properties from MolClusters 0.6.0 or earlier"
    Versions up to 0.6.0 computed several cluster properties wrongly. Use a later
    version, and recompute any results obtained with those versions:

    - **Dipole moments** were 1.89 times too small: they were converted to debye
      with the factor for e·a₀ instead of e·Å.
    - **Volumes** were 2.15 times too small (a sphere whose radius is the radius of
      gyration, not the uniform sphere with that radius of gyration).
    - **Densities** used that volume and an inverted unit factor (amu/Å³ to g/cm³,
      off by 2.76 times).
    - **Clusters wider than half the box** could be split apart when made whole,
      which distorted all their geometric properties.

    The **radius and diameter** also changed meaning: they were the radius of
    gyration and twice it, and are now the [equivalent sphere's](properties.md#size),
    so they can't be compared across versions.

    Which molecules form a cluster could differ too: `cm` rules used the plain
    center of mass of molecules split across the periodic boundaries, which could
    connect molecules that weren't close, or miss ones that were.

=== "uv (recommended)"

    [uv](https://docs.astral.sh/uv/) installs the `molclusters` command in an
    isolated environment of its own:

    ```bash
    uv tool install git+https://github.com/EmanuelMancio/MolClusters
    molclusters --version
    ```

    Upgrade it later with `uv tool upgrade molclusters`.

=== "pipx"

    ```bash
    pipx install git+https://github.com/EmanuelMancio/MolClusters
    molclusters --version
    ```

=== "pip"

    Into an environment you manage yourself (a virtual environment or a conda
    environment), for instance to use MolClusters from Python as well:

    ```bash
    python -m pip install git+https://github.com/EmanuelMancio/MolClusters
    molclusters --version
    ```

To install a given release rather than the latest code, add its tag,
`git+https://github.com/EmanuelMancio/MolClusters@v<version>`. The releases and what
changed in each are listed in the [changelog](../changelog.md).

## Dependencies

They are installed automatically:

| Package | Used for |
| --- | --- |
| [MDAnalysis](https://www.mdanalysis.org/) | reading topologies and trajectories, distances, hydrogen bonds |
| [NetworkX](https://networkx.org/) | the graph of connections between molecules |
| [NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/) | numbers and tables |
| [Pydantic](https://docs.pydantic.dev/), [PyYAML](https://pyyaml.org/) | reading and checking the configuration file |
| [orjson](https://github.com/ijl/orjson), [zstandard](https://github.com/indygreg/python-zstandard) | writing the JSON report |
| [loguru](https://github.com/Delgan/loguru), [tqdm](https://github.com/tqdm/tqdm) | logging and the progress bar |

## From a clone, for development

To work on MolClusters itself, see the
[development setup](../developer-guide/setup.md).
