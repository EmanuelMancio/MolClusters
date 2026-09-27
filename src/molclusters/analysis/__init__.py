# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""The analyses run on the tracked clusters, each behind the `FrameAnalysis` interface.

Modules
-------
    base: The `FrameAnalysis` interface.
    size: `SizeEvolution`, the number and sizes of the clusters over time (evo.txt).
"""

from .base import FrameAnalysis
from .size import SizeEvolution

__all__ = ["FrameAnalysis", "SizeEvolution"]
