# Developer guide

This guide is for people working on MolClusters itself: fixing bugs, adding
analyses, changing what it computes.

- **[Development setup](setup.md)**: clone, install, run the checks.
- **[Architecture](architecture.md)**: the modules, what each owns, and how a run
  flows through them.
- **[The tracking algorithm](tracking.md)**: how cluster ids are assigned, in the
  detail needed to change it.
- **[Adding a built-in analysis](builtin-analyses.md)**: the registry and the
  conventions an analysis follows.
- **[Changing the report format](report-format.md)**: versioning and the reader.
- **[Testing](testing.md)**: the test suite, its fixtures, and the units tests.
- **[Commits and releases](releases.md)**: commit messages, branches, what counts as
  a breaking change, and how a release is cut.
- **[Writing the docs](docs.md)**: this site, and how to build and publish it.

Pull requests are welcome; for larger changes, please
[open an issue](https://github.com/EmanuelMancio/MolClusters/issues) first to
discuss them.
