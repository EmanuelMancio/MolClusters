# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""The analyses run on the tracked clusters, each a `FrameAnalysis`.

Modules
-------
    base: `FrameAnalysis`, and the `Run` and `Frame` an analysis sees.
    builtins: `BUILTINS`, the built-in analyses in the order they run, and what
      each needs from the config.
    coordinates: `ClusterCoordinates`, the clusters holding solutes (.gro files).
    lineage: `Lineage`, where the clusters came from, what became of them and when
      (cluster_events.csv, cluster_lifetimes.csv).
    nucleus: `Nucleus`, the nuclei inside the clusters (nucleus_data.csv).
    report: `JsonReport`, every cluster of every frame (molclusters.jsonl).
    size: `SizeEvolution`, the number and sizes of the clusters over time (evo.txt).
    solute: `SoluteSolvent`, the clusters holding solutes and solvents
      (solute_solvent.csv).
"""

from .base import Frame, FrameAnalysis, Run
from .coordinates import ClusterCoordinates
from .lineage import Lineage
from .nucleus import Nucleus
from .report import JsonReport
from .size import SizeEvolution
from .solute import SoluteSolvent

__all__ = [
    "ClusterCoordinates",
    "Frame",
    "FrameAnalysis",
    "JsonReport",
    "Lineage",
    "Nucleus",
    "Run",
    "SizeEvolution",
    "SoluteSolvent",
]
