<!--
SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>

SPDX-License-Identifier: LGPL-3.0-or-later
-->

# Contributing to MolClusters

Thanks for your interest in MolClusters! Bug reports, questions and pull requests
are all welcome.

## Reporting a bug or asking a question

[Open an issue](https://github.com/EmanuelMancio/MolClusters/issues) with:

- the MolClusters version (`molclusters --version`) and your Python version;
- the command you ran and your configuration file;
- the full error message or the log (`molclusters_<date>.log` in the output
  directory);
- if you can, a small trajectory and topology that reproduce the problem.

## Changing the code

For anything larger than a small fix, please open an issue first so we can agree
on the approach. Then:

1. Fork the repository and create a branch off **`dev`** (`main` only receives
   releases); open your pull request against `dev`.
2. Set up the development environment with [uv](https://docs.astral.sh/uv/):

   ```bash
   git clone https://github.com/<your-username>/MolClusters.git
   cd MolClusters
   git switch dev
   uv sync
   uv run pre-commit install
   ```

3. Make your change, with tests, and run every check CI runs:

   ```bash
   uv run pre-commit run --all-files
   ```

4. Write commit messages in the
   [Conventional Commits](https://www.conventionalcommits.org/) format (for example
   `fix(config): reject an empty rule`); the `commit-msg` hook checks them.

The [developer guide](https://emanuelmancio.github.io/MolClusters/developer-guide/)
covers the architecture, the conventions (docstrings, license headers, units) and
what counts as a breaking change.

## License

By contributing, you agree that your contributions are licensed under the
[GNU LGPL v3 or later](LICENSE), like the rest of the project. Every file carries a
[REUSE](https://reuse.software/) license header, which the pre-commit hooks check.
