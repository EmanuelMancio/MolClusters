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
├── assets/, javascripts/    # favicon, figures, MathJax setup
```

A new page needs an entry in `nav` in `mkdocs.yml`. Pages are covered by the
`docs/**` annotation in `REUSE.toml`, so they need no license header of their own.

## What to update when

| When you change... | Update |
| --- | --- |
| a configuration key | [Configuration](../user-guide/configuration.md), and its figures for `rules` and `nucleus` |
| the command-line arguments or their help | [Command line](../user-guide/cli.md), which quotes `molclusters -h` |
| an output file or its columns | [Output files](../user-guide/outputs.md) |
| the report format | [The JSON report](../user-guide/report.md), and `report.py`'s docstring |
| a property's definition | [Cluster properties](../user-guide/properties.md) and its figures |
| id assignment | [How clusters are tracked](../user-guide/concepts.md) and its figures, [The tracking algorithm](tracking.md) |
| a public class or function | its docstring (the reference follows); a new module gets a reference page |

## Figures

Each figure is a pair of SVG files, `<name>-light.svg` and `<name>-dark.svg`, shown
with Material's `#only-light` and `#only-dark` and the same alt text, which
describes what the figure shows for readers who can't see it:

```markdown
![Alt text.](../assets/concepts/events-light.svg#only-light)
![Alt text.](../assets/concepts/events-dark.svg#only-dark)
```

The schematic figures of [How clusters are tracked](../user-guide/concepts.md),
[Configuration](../user-guide/configuration.md) (the rules, nuclei) and
[Cluster properties](../user-guide/properties.md) (size, shape) and the home page's
steps, all in
`docs/assets/concepts/`, are drawn by `scripts/docs_figures.py`; change them there
and rerun it rather than editing the SVG files:

```bash
uv run python scripts/docs_figures.py
```

The home page's picture (`docs/assets/home/droplet.webp`) is frame 990 of the
`tests/data/met-mal` run with the quick start's rules, ray-traced by VMD 1.9.3's
Tachyon with the scripts in `scripts/docs_snapshot/`. Tachyon writes no alpha
channel, so the scene is rendered on white and on black and `matte.py` recovers
the transparency from the two. VMD runs on Linux (or WSL); the rest anywhere:

```bash
uv run molclusters tests/data/met-mal/met-traj.trr tests/data/met-mal/met-mal.tpr \
    config.yml --output-dir run   # the quick start's config.yml, with report_compression: none
uv run python scripts/docs_snapshot/frame.py run/molclusters.jsonl 990 work/frame990.pdb
sh scripts/docs_snapshot/render.sh work/frame990.pdb work          # in Linux or WSL
uv run --with pillow python scripts/docs_snapshot/matte.py work/black.tga work/white.tga \
    docs/assets/home/droplet.webp
```

Which six clusters are coloured is set in `scene.tcl`, for this frame and view.

The plots of the quick start (`evo`) and of [Output files](../user-guide/outputs.md)
(lifetimes, solute-solvent clusters) come from the quick start's run on all frames
of the test data, with its `config.yml` as it stands (coordinates included):

```bash
uv run --with matplotlib python scripts/docs_plots.py run
```

The Output files page's picture of `cls-id2872.gro` is four of that file's frames,
side by side, labelled with their time and size:

```bash
uv run python scripts/docs_snapshot/cluster_row.py run/coordinates/cls-id2872.gro \
    work/row.pdb 40 11328 11664 11904 14736
vmd -dispdev text -e scripts/docs_snapshot/cluster.tcl -args work/row.pdb work/scene.dat
# Tachyon on white and on black as in render.sh, at -res 3600 3600, then:
uv run --with pillow python scripts/docs_snapshot/matte.py work/black.tga work/white.tga work/row.png
uv run --with pillow python scripts/docs_snapshot/label_row.py work/row.png docs/assets/outputs/cls-id \
    "11.33 ns|4 molecules, new" "11.66 ns|13 molecules" "11.90 ns|17, its largest" "14.74 ns|4, its last frame"
```

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
