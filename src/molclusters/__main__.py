# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""This is the entry point for the `molclusters` package.

The module imports and executes the `main` function from the `main` module,
serving as the starting point for the application when run as a script.

Usage:
    To execute the `molclusters` application, run this module directly.

Example:
    molclusters -h
"""

from .main import main

main()
