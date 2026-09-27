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
    tracker: The `ClusterTracker`, which keeps cluster ids stable over a trajectory.
    analysis: The analyses run on the tracked clusters, and their `FrameAnalysis` base.
    output: The files the analyses write.

Exports
-------
    MolClusters: The primary class for performing cluster analysis.
    ClusterTracker: Follows the clusters of a trajectory, without analysing them.
    FrameAnalysis: The base class of an analysis, to pass to `MolClusters`.
    Frame, Run: What an analysis sees of the current frame and of the run.
    OutputFile: The declaration of a file an analysis writes.
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
from .analysis import Frame, FrameAnalysis, Run
from .log import start_logging
from .molclusters import MolClusters
from .output import OutputFile
from .tracker import ClusterTracker
from .version import __version__, version

__all__ = [
    "version",
    "__version__",
    "start_logging",
    "MolClusters",
    "ClusterTracker",
    "FrameAnalysis",
    "Frame",
    "Run",
    "OutputFile",
]
