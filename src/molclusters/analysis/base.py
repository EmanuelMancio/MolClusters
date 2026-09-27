# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides `FrameAnalysis`, the base of the analyses run on the tracked clusters.

Classes:
--------
- FrameAnalysis: The base class of an analysis, run frame by frame.
- Run: What an analysis sees of the run as a whole.
- Frame: What an analysis sees of the current frame.
"""

from abc import ABC, abstractmethod
from collections.abc import Generator, Mapping, Sequence
from contextlib import contextmanager
from types import MappingProxyType

import MDAnalysis as mda

from ..cluster import Cluster
from ..config import MolClsConfig
from ..output import OutputFile, RunOutput
from ..tracker import ClusterTracker


class FrameAnalysis(ABC):
    """An analysis of the tracked clusters, run frame by frame.

    A run calls `prepare` once, then `analyse` for every frame (the first one
    included) once the clusters are up to date with it, and finally `finish`.
    The analyses of a run are called in order, so an analysis sees the results of
    the ones before it for the same frame.

    Only `analyse` must be implemented. An analysis keeps its results as its own
    attributes; one that writes files does so through `Run.output` and declares
    them in `outputs`.

    Attributes
    ----------
    outputs : tuple[OutputFile, ...]
        The files the analysis writes, to warn about results of an earlier run
        they would overwrite or be appended to, and to report them at the end.
        Usually a class attribute; an analysis whose files depend on its options
        sets it in ``__init__`` instead, since it's read before `prepare`.
    """

    outputs: tuple[OutputFile, ...] = ()

    def prepare(self, run: "Run") -> None:  # noqa: B027 (optional hook)
        """Get ready for a run. Does nothing unless overridden.

        Parameters
        ----------
        run : Run
            The run about to start.
        """

    @abstractmethod
    def analyse(self, frame: "Frame") -> None:
        """Analyse the clusters of the current frame.

        Parameters
        ----------
        frame : Frame
            The current frame.
        """

    def finish(self, run: "Run") -> None:  # noqa: B027 (optional hook)
        """Wrap up once every frame was analysed. Does nothing unless overridden.

        Parameters
        ----------
        run : Run
            The run that just ended.
        """


class Run:
    """What an analysis sees of the run as a whole.

    Attributes
    ----------
    universe : mda.Universe
        The Universe being analysed.
    config : MolClsConfig
        The analysis configuration.
    n_frames : int
        The number of frames of the run.
    output : RunOutput
        Where the analyses write their files.
    """

    __slots__ = ["universe", "config", "n_frames", "output", "_analyses", "_visible"]

    def __init__(
        self,
        universe: mda.Universe,
        config: MolClsConfig,
        n_frames: int,
        output: RunOutput,
        analyses: Sequence[FrameAnalysis],
    ) -> None:
        """Initialize the context of a run.

        Parameters
        ----------
        universe : mda.Universe
            The Universe being analysed.
        config : MolClsConfig
            The analysis configuration.
        n_frames : int
            The number of frames of the run.
        output : RunOutput
            Where the analyses write their files.
        analyses : Sequence[FrameAnalysis]
            The analyses of the run, in the order they're called.
        """
        self.universe = universe
        self.config = config
        self.n_frames = n_frames
        self.output = output
        self._analyses = analyses
        self._visible = len(analyses)

    def prepare_analyses(self) -> None:
        """Prepare every analysis of the run, in order (called by the runner)."""
        for i, analysis in enumerate(self._analyses):
            self._visible = i
            with _blame(analysis, "prepare"):
                analysis.prepare(self)
        self._visible = len(self._analyses)

    def analyse_frame(self, frame: "Frame") -> None:
        """Run every analysis of the run on `frame`, in order (called by the runner).

        Parameters
        ----------
        frame : Frame
            The current frame.
        """
        for analysis in self._analyses:
            with _blame(analysis, "analyse", frame.index):
                analysis.analyse(frame)

    def finish_analyses(self) -> None:
        """Finish every analysis of the run, in order (called by the runner)."""
        for analysis in self._analyses:
            with _blame(analysis, "finish"):
                analysis.finish(self)

    def analysis[T: FrameAnalysis](self, kind: type[T]) -> T | None:
        """Find the first analysis of the run of type `kind`, for its results.

        While the analyses are being prepared, only those before the one being
        prepared are found, since only their results for a frame are ready when
        it analyses the frame.

        Parameters
        ----------
        kind : type[T]
            The analysis class to look for (subclasses match too).

        Returns
        -------
        T | None
            The analysis, or None if the run has none (before this one).
        """
        for analysis in self._analyses[: self._visible]:
            if isinstance(analysis, kind):
                return analysis
        return None


class Frame:
    """What an analysis sees of the current frame.

    Attributes
    ----------
    index : int
        The index of the frame in the run (0 for its first frame).
    time : float
        The time of the frame, as the trajectory reports it.
    universe : mda.Universe
        The Universe, positioned at this frame.
    """

    __slots__ = ["index", "time", "universe", "_tracker"]

    def __init__(self, index: int, tracker: ClusterTracker) -> None:
        """Initialize the current frame of a run.

        Parameters
        ----------
        index : int
            The index of the frame in the run.
        tracker : ClusterTracker
            The tracker, up to date with the frame.
        """
        self.index = index
        self.time: float = tracker.uni.coord.time
        self.universe = tracker.uni
        self._tracker = tracker

    @property
    def clusters(self) -> Mapping[int, Cluster]:
        """The clusters of the frame, keyed by cluster id (read-only).

        Returns
        -------
        Mapping[int, Cluster]
            A read-only view of the clusters; the clusters must not be modified.
        """
        return MappingProxyType(self._tracker.clusters)

    def find(self, mol: int) -> int | bool:
        """Find the id of the cluster a molecule belongs to.

        Parameters
        ----------
        mol : int
            The molecule's resid.

        Returns
        -------
        int | bool
            The cluster id, or False if the molecule is in no cluster.
        """
        return self._tracker.find(mol)


@contextmanager
def _blame(
    analysis: FrameAnalysis, hook: str, frame: int | None = None
) -> Generator[None, None, None]:
    """Note which analysis raised an error, so it can be told from a bug of the run.

    Parameters
    ----------
    analysis : FrameAnalysis
        The analysis being called.
    hook : str
        The name of the method being called.
    frame : int | None
        The index of the frame being analysed, if any.

    Yields
    ------
    None
        Control, for the call to the analysis.
    """
    try:
        yield
    except Exception as err:
        where = "" if frame is None else f" on frame {frame}"
        err.add_note(f"Raised by {type(analysis).__name__}.{hook}(){where}")
        raise
