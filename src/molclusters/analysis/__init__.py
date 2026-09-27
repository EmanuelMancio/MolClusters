# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""The analyses run on the tracked clusters, each a `FrameAnalysis`.

Modules
-------
    base: `FrameAnalysis`, and the `Run` and `Frame` an analysis sees.
    coordinates: `ClusterCoordinates`, the clusters holding solutes (.gro files).
    size: `SizeEvolution`, the number and sizes of the clusters over time (evo.txt).
    solute: `SoluteSolvent`, the clusters holding solutes and solvents
      (solute_solvent.csv).
"""

from .base import Frame, FrameAnalysis, Run
from .coordinates import ClusterCoordinates
from .size import SizeEvolution
from .solute import SoluteSolvent

__all__ = [
    "ClusterCoordinates",
    "Frame",
    "FrameAnalysis",
    "Run",
    "SizeEvolution",
    "SoluteSolvent",
]
