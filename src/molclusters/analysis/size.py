# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides `SizeEvolution`, the number and sizes of the clusters over time (evo.txt)."""

import numpy as np
from loguru import logger

from ..tracker import ClusterTracker


class SizeEvolution:
    """Records the number of clusters and their sizes, frame by frame.

    Attributes
    ----------
    data : np.ndarray
        One row per frame: time, number of clusters, and the minimum, average and
        maximum cluster size (all 0 in a frame without clusters).
    """

    __slots__ = ["data"]

    def __init__(self) -> None:
        """Initialize an empty record, sized by `prepare`."""
        self.data = np.zeros((0, 5))

    def prepare(self, tracker: ClusterTracker, n_frames: int) -> None:  # noqa: ARG002
        """Make room for `n_frames` frames.

        Parameters
        ----------
        tracker : ClusterTracker
            The tracker whose clusters will be analysed.
        n_frames : int
            The number of frames of the run.
        """
        self.data = np.zeros((n_frames, 5))

    def analyse(self, tracker: ClusterTracker, frame: int) -> None:
        """Update cluster size evolution information for a given frame.

        Parameters
        ----------
        tracker : ClusterTracker
            The tracker, up to date with the current frame.
        frame : int
            The frame index.
        """
        sizes = [cls.size for cls in tracker.clusters.values()]
        if len(sizes) == 0:
            avg, min_size, max_size = 0, 0, 0
        else:
            avg = np.average(sizes)
            min_size = min(sizes)
            max_size = max(sizes)

        time = tracker.uni.coord.time

        self.data[frame][0] = time
        self.data[frame][1] = len(tracker.clusters)
        self.data[frame][2] = min_size
        self.data[frame][3] = avg
        self.data[frame][4] = max_size

    def finish(self) -> list[str]:
        """Log a summary of the cluster counts and sizes and write them to evo.txt.

        Returns
        -------
        list[str]
            The names of the files written.
        """
        n_clusters = self.data[:, 1]
        logger.info(
            f"Found {n_clusters.mean():.1f} cluster(s) per frame on average "
            f"({n_clusters.min():.0f}-{n_clusters.max():.0f}); the largest held "
            f"{self.data[:, 4].max():.0f} molecule(s)."
        )

        np.savetxt(
            "evo.txt",
            self.data,
            header="Time NClusters MinSize AvgSize MaxSize",
        )
        return ["evo.txt"]
