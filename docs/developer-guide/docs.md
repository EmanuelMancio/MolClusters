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

`.github/workflows/docs.yml` builds the site on every pull request (to catch
errors) and on pushes to `main`, and deploys the `main` build to GitHub Pages at
<https://emanuelmancio.github.io/MolClusters/>. It can also be run by hand from the
repository's *Actions* tab.

The repository's *Settings → Pages → Build and deployment → Source* must be set to
**GitHub Actions** once, for the deployment to be accepted.
