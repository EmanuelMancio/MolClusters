# Commits and releases

## Commit messages

Commit messages follow [Conventional Commits](https://www.conventionalcommits.org/)
with the Angular types, checked by the `commitizen` commit-msg hook:

```text
<type>(<optional scope>): <description>

<optional body>

<optional footer>
```

Types: `feat`, `fix`, `docs`, `style`, `refactor`, `perf`, `test`, `build`, `ci`,
`chore`, `revert`. For example `fix(config): reject an empty rule` or
`feat(analysis): add Lineage, the clusters' events and lifetimes`. Commitizen reads
them to choose the next version and to write `CHANGELOG.md`: `fix` bumps the patch
version, `feat` the minor one, and a breaking change the major one (the minor one
while the version is 0.x).

## Breaking changes

A change to clustering, id assignment or any other computed output is a
**breaking change** (`type(scope)!:` and a `BREAKING CHANGE:` footer, even for a
`fix`) only if it changes a result that was previously valid and expected, e.g. a
change to the tie-breaking that picks a different, equally plausible, id assignment.

Correcting output that was simply wrong — a bug nobody would have expected or relied
on, such as deriving correct atom names for LAMMPS topologies — is a normal `fix`,
even though the stored values change.

Changes to the report format bump its own `VERSION` (see
[Changing the report format](report-format.md)); whether they are also a breaking
change for the package follows the same rule.

## Branches

- `dev` is where work happens; feature branches start from it and pull requests
  target it. Dependabot's updates also go to `dev`.
- `main` receives releases: `dev` is merged into it when a release is due.

## Cutting a release

Releases are cut by hand, when the changes on `dev` are worth one: nothing bumps the
version on its own, and commits simply collect under *Unreleased* until then. When
it's time, merge `dev` into `main` and, on `main` with the checks passing:

```bash
uv run cz bump --dry-run   # optional: see the version and changelog it would make
uv run cz bump
```

Commitizen computes the next version from the commits since the last tag, updates
`pyproject.toml`, `uv.lock`, `CITATION.cff` and `CHANGELOG.md`, and creates a
GPG-signed commit (`chore: release vX -> vY`) and an annotated tag `vY`. Its
pre-bump hooks refuse to run off `main` (`scripts/check_release_branch.py`), refresh
`uv.lock`, and set `CITATION.cff`'s `version` and `date-released`
(`scripts/update_citation.py`). Then push the branch and the tag:

```bash
git push origin main --follow-tags
```

and merge `main` back into `dev`. The push to `main` runs the tests and publishes this
documentation site (see [Writing the docs](docs.md)), and the tag runs the *Release*
workflow (`.github/workflows/release.yml`):

1. **build**: checks the tag matches the package version, runs the tests, builds the
   wheel and sdist and checks their metadata (`twine check`);
2. **pypi**: uploads them to [PyPI](https://pypi.org/project/molclusters/);
3. **github-release**: publishes the GitHub release, with the same files attached and
   the version's `CHANGELOG.md` section as its notes.

If it fails for a passing reason (a network error), re-run it from the Actions tab; if
the tagged code itself fails, don't move the tag: fix it on `dev` and cut the next
version. A version uploaded to PyPI can never be uploaded again, even after deleting
it, which is why the tests run before the upload.

## PyPI

The *pypi* job uses [trusted publishing](https://docs.pypi.org/trusted-publishers/):
PyPI trusts this repository's workflow directly, so no API token is stored anywhere.
It was set up once, before the first upload:

1. On PyPI, *Your account → Publishing → Add a new pending publisher*, with project
   name `molclusters`, owner `EmanuelMancio`, repository `MolClusters`, workflow
   `release.yml` and environment `pypi`. The first upload then creates the project.
2. On GitHub, *Settings → Environments → New environment* `pypi`; optionally with
   yourself as required reviewer, so every upload waits for your approval.

The PyPI page shows `README.md`, so its links must be absolute URLs: relative ones
work on GitHub but not on PyPI.

The release commit must not carry `[skip ci]`: GitHub would then skip every workflow
triggered by the push, the tag's included.

## DOIs (Zenodo)

Each published GitHub release, as the *Release* workflow makes, is archived on
[Zenodo](https://zenodo.org/), which gives it a DOI; a tag alone isn't. Zenodo takes the record's title, description, authors
(with their ORCID and affiliation), keywords and license from `CITATION.cff`, and the
version from the release's tag, so the authors in `CITATION.cff` are the ones credited.
A Zenodo record can't be deleted: check `CITATION.cff` before publishing a release.

The first archived release also creates a *concept DOI*, which always resolves to the
latest version. Add it to `CITATION.cff` as `doi:` and as a badge in the README, so that
citations can use the concept DOI or a specific version's one.
