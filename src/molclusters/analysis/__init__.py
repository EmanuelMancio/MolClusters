# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""The analyses run on the tracked clusters, each a `FrameAnalysis`.

Modules
-------
    base: `FrameAnalysis`, and the `Run` and `Frame` an analysis sees.
    size: `SizeEvolution`, the number and sizes of the clusters over time (evo.txt).
"""

from .base import Frame, FrameAnalysis, Run
from .size import SizeEvolution

__all__ = ["Frame", "FrameAnalysis", "Run", "SizeEvolution"]
