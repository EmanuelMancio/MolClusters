# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""MolClusters: A Python tool for analyzing cluster formation in molecular dynamics simulations.

This module provides functionality to study and analyze the formation of molecular clusters
from simulation data. It includes tools for processing and analyzing cluster
dynamics, enabling researchers to gain insights into molecular behavior of clusters.

Attributes
----------
    __version__ (str): The version of the MolClusters package.
    version (str): The version of the MolClusters package.

Modules
-------
    cluster: Submodule containing utilities for cluster analysis.
    molclusters: Core module defining the main `MolClusters` class.

Exports
-------
    MolClusters: The primary class for performing cluster analysis.
    start_logging: Function to set up molclusters logger.

Example
-------
    >>> from molclusters import MolClusters
    >>> analyzer = MolClusters(universe, config)
    >>> analyzer.run()
"""

from loguru import logger

logger.disable("molclusters")

from . import cluster as cluster
from .log import start_logging
from .molclusters import MolClusters
from .version import __version__, version

__all__ = ["version", "__version__", "start_logging", "MolClusters"]
