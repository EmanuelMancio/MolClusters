# Writing the docs

This site is built with [MkDocs](https://www.mkdocs.org/) and
[Material for MkDocs](https://squidfunk.github.io/mkdocs-material/); the API
reference is generated from the docstrings by
[mkdocstrings](https://mkdocstrings.github.io/). Its tools are in the `docs`
dependency group.

```bash
uv run --group docs mkdocs serve      # live preview at http://127.0.0.1:8000
uv run --group docs mkdocs build --strict
```

`--strict` turns warnings (broken links, missing pages, docstrings that can't be
parsed) into errors; CI builds with it.

## Layout

```text
mkdocs.yml                   # site configuration and navigation
docs/
├── index.md                 # home page
├── user-guide/              # for people analysing simulations
├── developer-guide/         # for people working on MolClusters
├── reference/               # API reference: one page per module, `::: module` directives
├── changelog.md             # includes CHANGELOG.md
├── assets/, javascripts/    # favicon, MathJax setup
```

A new page needs an entry in `nav` in `mkdocs.yml`. Pages are covered by the
`docs/**` annotation in `REUSE.toml`, so they need no license header of their own.

## What to update when

| When you change... | Update |
| --- | --- |
| a configuration key | [Configuration](../user-guide/configuration.md) |
| the command-line arguments or their help | [Command line](../user-guide/cli.md), which quotes `molclusters -h` |
| an output file or its columns | [Output files](../user-guide/outputs.md) |
| the report format | [The JSON report](../user-guide/report.md), and `report.py`'s docstring |
| a property's definition | [Cluster properties](../user-guide/properties.md) |
| id assignment | [How clusters are tracked](../user-guide/concepts.md), [The tracking algorithm](tracking.md) |
| a public class or function | its docstring (the reference follows); a new module gets a reference page |

## Docstrings

Docstrings use NumPy style and are rendered as they are, so the reference is only as
good as they are. Backticked names (`` `ClusterTracker` ``) render as code; to link
to another object in a docstring, use mkdocstrings' cross-reference syntax,
`` [`Cluster`][molclusters.cluster.Cluster] ``. Private names (a leading
underscore) are left out of the reference.

## Publishing

`.github/workflows/docs.yml` builds the site with `--strict` on every pull request and
push (to catch errors). On pushes to `main` and `dev` it also deploys it with
[mike](https://github.com/jimporter/mike) to the `gh-pages` branch: `main` as one
folder per minor version (`0.7`, `0.8`, ...) taken from `pyproject.toml`, copied to a
`latest` folder, and `dev` as a separate `dev` version (the unreleased docs). Pushing
docs-only changes to `main` therefore updates the current version's pages, and a
release of a new minor adds a folder. Every version's canonical URLs point at
`latest`, so search engines index that one. The site is at
<https://emanuelmancio.github.io/MolClusters/>; the version selector in the header
switches between versions. Links from before the site was versioned
(`.../MolClusters/user-guide/...`, as in the v0.7.0 README on PyPI and Zenodo) still
work: the workflow also puts `.github/gh-pages/404.html` at the branch's root, which
sends a missing path to the same page under `latest/`, so links in the README and
elsewhere should point at `latest/` directly. The workflow can also be run by hand
from the repository's *Actions* tab.

Set up once: *Settings → Pages → Build and deployment → Source* to **Deploy from a
branch**, `gh-pages` / root; and after the first deployment from `main` (which creates
`latest`) make the bare URL open it:

```bash
uv run --group docs mike set-default --push latest
```

To preview the versioned site locally, `uv run --group docs mike serve` (it serves
the `gh-pages` branch, so it shows nothing until a deployment exists). An old
version's pages are fixed once built; delete one with `mike delete --push <version>`.
