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

import time
from abc import ABC, abstractmethod
from collections.abc import Generator, Mapping, Sequence
from contextlib import contextmanager
from types import MappingProxyType

import MDAnalysis as mda
from loguru import logger

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
    them in `outputs`. The same analysis may be run more than once, so it should
    start its results over in `prepare`, not only in ``__init__``.

    An analysis logs with loguru's ``logger``, as the rest of the package does;
    the run tags whatever is logged during its hooks with its class name (the
    ``analysis`` key of the record's ``extra``), which the log shows. To keep the
    log readable over long trajectories:

    - `analyse` counts what's worth reporting (e.g. in a `collections.Counter`)
      instead of logging it, and `finish` logs a summary: one INFO line, or one
      warning saying how to fix the problem.
    - Anything logged per frame is DEBUG or TRACE, with loguru's own arguments
      (``logger.debug("cluster {}: ...", cls.id)``) rather than an f-string: a
      message filtered out by the log level is then never built.
    - A warning that could repeat is only logged once.

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
        """Get ready for a run, starting over any results. Does nothing unless overridden.

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
    durations : list[float]
        The time spent in each analysis' hooks so far, in seconds, in the order of
        the analyses.
    """

    __slots__ = [
        "universe",
        "config",
        "n_frames",
        "output",
        "durations",
        "_analyses",
        "_visible",
    ]

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
        self.durations = [0.0] * len(analyses)
        self._analyses = analyses
        self._visible = len(analyses)

    def prepare_analyses(self) -> None:
        """Prepare every analysis of the run, in order (called by the runner)."""
        for i, analysis in enumerate(self._analyses):
            self._visible = i
            with self._hook(i, "prepare"):
                analysis.prepare(self)
        self._visible = len(self._analyses)

    def analyse_frame(self, frame: "Frame") -> None:
        """Run every analysis of the run on `frame`, in order (called by the runner).

        Parameters
        ----------
        frame : Frame
            The current frame.
        """
        for i, analysis in enumerate(self._analyses):
            with self._hook(i, "analyse", frame.index):
                analysis.analyse(frame)

    def finish_analyses(self) -> None:
        """Finish every analysis of the run, in order (called by the runner)."""
        for i, analysis in enumerate(self._analyses):
            with self._hook(i, "finish"):
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

    @contextmanager
    def _hook(
        self, i: int, hook: str, frame: int | None = None
    ) -> Generator[None, None, None]:
        """Call a hook of the `i`-th analysis: tag its log, time it, blame its errors.

        What's logged during the call is tagged with the analysis' class name, and
        an error it raises is noted with it, so either can be told from the run's
        own; the time taken is added to `durations`.

        Parameters
        ----------
        i : int
            The index of the analysis being called.
        hook : str
            The name of the method being called.
        frame : int | None
            The index of the frame being analysed, if any.

        Yields
        ------
        None
            Control, for the call to the analysis.
        """
        name = type(self._analyses[i]).__name__
        start = time.perf_counter()
        try:
            with logger.contextualize(analysis=name):
                yield
        except Exception as err:
            where = "" if frame is None else f" on frame {frame}"
            err.add_note(f"Raised by {name}.{hook}(){where}")
            raise
        finally:
            self.durations[i] += time.perf_counter() - start


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
    output : RunOutput
        Where the run's analyses write their files (the same as `Run.output`),
        e.g. to append this frame to a file.
    """

    __slots__ = ["index", "time", "universe", "output", "_tracker"]

    def __init__(self, index: int, tracker: ClusterTracker, output: RunOutput) -> None:
        """Initialize the current frame of a run.

        Parameters
        ----------
        index : int
            The index of the frame in the run.
        tracker : ClusterTracker
            The tracker, up to date with the frame.
        output : RunOutput
            Where the run's analyses write their files.
        """
        self.index = index
        self.time: float = tracker.uni.coord.time
        self.universe = tracker.uni
        self.output = output
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

    def find(self, mol: int) -> int | None:
        """Find the id of the cluster a molecule belongs to.

        Parameters
        ----------
        mol : int
            The molecule's resid.

        Returns
        -------
        int | None
            The cluster id, or None if the molecule is in no cluster.
        """
        return self._tracker.find(mol)
