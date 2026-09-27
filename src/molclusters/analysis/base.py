# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides `FrameAnalysis`, the interface of the analyses run on the tracked clusters."""

from typing import Protocol

from ..tracker import ClusterTracker


class FrameAnalysis(Protocol):
    """An analysis of the tracked clusters, run frame by frame.

    A run calls `prepare` once, then `analyse` for every frame (the first one
    included) once the tracker is up to date with it, and finally `finish`.
    """

    def prepare(self, tracker: ClusterTracker, n_frames: int) -> None:
        """Get ready to analyse `n_frames` frames.

        Parameters
        ----------
        tracker : ClusterTracker
            The tracker whose clusters will be analysed.
        n_frames : int
            The number of frames of the run.
        """
        ...

    def analyse(self, tracker: ClusterTracker, frame: int) -> None:
        """Analyse the clusters of the current frame.

        Parameters
        ----------
        tracker : ClusterTracker
            The tracker, up to date with the current frame.
        frame : int
            The index of the current frame in the run.
        """
        ...

    def finish(self) -> list[str]:
        """Write the results to the working directory.

        Returns
        -------
        list[str]
            The names of the files written, for the run's summary.
        """
        ...
