# Installation

MolClusters needs **Python 3.12 or newer**, and is published on
[PyPI](https://pypi.org/project/molclusters/) from v0.7.0.

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

    [Known issues in earlier versions](known-issues.md) lists every error, the
    outputs it affects and how to correct old results.

=== "uv (recommended)"

    [uv](https://docs.astral.sh/uv/) installs the `molclusters` command in an
    isolated environment of its own:

    ```bash
    uv tool install molclusters
    molclusters --version
    ```

    Upgrade it later with `uv tool upgrade molclusters`.

=== "pipx"

    ```bash
    pipx install molclusters
    molclusters --version
    ```

=== "pip"

    Into an environment you manage yourself (a virtual environment or a conda
    environment), for instance to use MolClusters from Python as well:

    ```bash
    python -m pip install molclusters
    molclusters --version
    ```

To install a given release, pin its version, e.g. `molclusters==0.7.0`; the releases
and what changed in each are listed in the [changelog](../changelog.md). The latest
development version, not yet released, installs from the repository's `dev` branch:
`git+https://github.com/EmanuelMancio/MolClusters@dev` in place of `molclusters`.

## Container image

Each release from the one after v0.8.0 is also published as a container image on the
GitHub Container Registry, `ghcr.io/emanuelmancio/molclusters`, for computers where
installing Python packages is awkward, such as HPC clusters. It's a standard
[OCI](https://opencontainers.org/) image (linux/amd64 and linux/arm64), so Docker,
Podman, Apptainer and Singularity all run it. It holds only MolClusters and its
dependencies, at the exact versions it was tested with, and runs the `molclusters`
command with the arguments you give it.

Tags: `X.Y.Z` for a release (e.g. `0.9.0`), `X.Y` for the latest release of a minor
version, and `latest`. Pin a release (`X.Y.Z`) in anything you need to reproduce, or
better its digest (`ghcr.io/emanuelmancio/molclusters@sha256:...`), which `docker
inspect` or the package page on GitHub shows and which can't change.

=== "Apptainer / Singularity"

    Convert the image into a single `.sif` file once, then run it like a command.
    Apptainer runs it as you, in the current directory, with your home directory
    and the current directory available, so paths work as usual:

    ```bash
    apptainer pull molclusters.sif docker://ghcr.io/emanuelmancio/molclusters:X.Y.Z
    apptainer run molclusters.sif traj.xtc topol.tpr input.yaml
    ```

    `singularity` takes the same commands. Data outside your home and the current
    directory (e.g. a scratch file system) needs binding: `apptainer run --bind
    /scratch molclusters.sif ...`. In a batch job, put the `apptainer run` line
    where the `molclusters` command would go.

=== "Docker / Podman"

    The image's working directory is `/data`: mount the directory with your files
    there, and run as yourself, so that MolClusters can write its results there and
    they are yours:

    ```bash
    docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/data" \
        ghcr.io/emanuelmancio/molclusters:X.Y.Z traj.xtc topol.tpr input.yaml
    ```

    Without `--user` it runs as uid 1000, never as root. With rootless Podman, use
    `--userns=keep-id` in place of `--user`:

    ```bash
    podman run --rm --userns=keep-id -v "$PWD:/data" \
        ghcr.io/emanuelmancio/molclusters:X.Y.Z traj.xtc topol.tpr input.yaml
    ```

    Paths given to `molclusters` are paths inside the container: files outside the
    mounted directory need their own `-v` (`:ro` for inputs it only reads).

`molclusters --help` is the image's default command, and `--version` tells which
release it holds. The [command-line options](cli.md) are the same as for an
installed `molclusters`. Use `--traj-memory` only if the container may use the
memory it needs (Docker's `--memory`, or the job's limit).

Each image is signed with a build provenance attestation: [GitHub's
CLI](https://cli.github.com/) checks that it was built by MolClusters' own release
workflow, from the repository's sources, and not altered since:

```bash
gh attestation verify oci://ghcr.io/emanuelmancio/molclusters:X.Y.Z -R EmanuelMancio/MolClusters
```

To build the image yourself, e.g. from the `dev` branch, run `docker build -t
molclusters .` (or `podman build`) in a clone of the repository.

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
