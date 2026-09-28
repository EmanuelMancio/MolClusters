# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Provides `SizeEvolution`, the number and sizes of the clusters over time (evo.txt)."""

import numpy as np
from loguru import logger

from ..output import OutputFile
from .base import Frame, FrameAnalysis, Run


class SizeEvolution(FrameAnalysis):
    """Records the number of clusters and their sizes, frame by frame.

    Attributes
    ----------
    data : np.ndarray
        One row per frame: time, number of clusters, and the minimum, average and
        maximum cluster size (all 0 in a frame without clusters).
    """

    __slots__ = ["data"]

    outputs = (OutputFile("evo.txt"),)

    def __init__(self) -> None:
        """Initialize an empty record, sized by `prepare`."""
        self.data = np.zeros((0, 5))

    def prepare(self, run: Run) -> None:
        """Make room for every frame of the run.

        Parameters
        ----------
        run : Run
            The run about to start.
        """
        self.data = np.zeros((run.n_frames, 5))

    def analyse(self, frame: Frame) -> None:
        """Record the number of clusters and their sizes in the current frame.

        Parameters
        ----------
        frame : Frame
            The current frame.
        """
        sizes = [cls.size for cls in frame.clusters.values()]
        if len(sizes) == 0:
            avg, min_size, max_size = 0, 0, 0
        else:
            avg = np.average(sizes)
            min_size = min(sizes)
            max_size = max(sizes)

        self.data[frame.index][0] = frame.time
        self.data[frame.index][1] = len(sizes)
        self.data[frame.index][2] = min_size
        self.data[frame.index][3] = avg
        self.data[frame.index][4] = max_size

    def finish(self, run: Run) -> None:
        """Log a summary of the cluster counts and sizes and write them to evo.txt.

        Parameters
        ----------
        run : Run
            The run that just ended.
        """
        n_clusters = self.data[:, 1]
        logger.info(
            f"Found {n_clusters.mean():.1f} cluster(s) per frame on average "
            f"({n_clusters.min():.0f}-{n_clusters.max():.0f}); the largest held "
            f"{self.data[:, 4].max():.0f} molecule(s)."
        )

        np.savetxt(
            run.output.path("evo.txt"),
            self.data,
            header="Time NClusters MinSize AvgSize MaxSize",
        )
