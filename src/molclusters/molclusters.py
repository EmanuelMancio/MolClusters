# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides the `MolClusters` class for analyzing molecular clusters in molecular dynamics simulations.

Classes:
--------
- MolClusters: Runs the analyses the config enables (and any given ones) on the
  clusters of every frame of a trajectory.

The clusters themselves, and their ids, are tracked by `tracker.ClusterTracker`;
the analyses are in the `analysis` package.

Dependencies:
-------------
- MDAnalysis: For molecular dynamics trajectory and structure analysis.
- NumPy: For numerical computations.
- tqdm: For progress tracking during analysis.
"""

import pathlib as path
import time
from collections.abc import Iterable

import MDAnalysis as mda
from loguru import logger
from tqdm import tqdm

from .analysis import Frame, FrameAnalysis, Run
from .analysis.builtins import build_builtins
from .config import MolClsConfig
from .conntable import check_resids
from .log import FILE_ONLY, format_duration
from .output import OutputFile, RunOutput
from .tracker import ClusterTracker


class MolClusters:
    """A class for analyzing molecular clusters in molecular dynamics simulations.

    It tracks the clusters frame by frame and runs the analyses on them: the
    built-ins the config enables (cluster sizes, solute-solvent, coordinates,
    nuclei, the JSON report; see `analysis.builtins`), then the ones given to it.
    A built-in runs when the config options it needs are set, unless the
    config's `analyses` turns it off.

    Attributes
    ----------
    uni : mda.Universe
        The MDAnalysis Universe object associated with the simulation.
    config : MolClsConfig
        The analysis configuration, with the defaults that depend on the
        topology filled in (the solvent).
    tracker : ClusterTracker | None
        Follows the clusters frame by frame, keeping their ids stable; created by
        `run` (None before), and left at the last frame.
    analyses : list[FrameAnalysis]
        The analyses run on the clusters of every frame: the built-ins the config
        enables, then the given ones (see `analysis` to find one).
    output : RunOutput
        Where the last run wrote its files.
    """

    __slots__ = [
        "uni",
        "config",
        "tracker",
        "analyses",
        "output",
    ]

    def __init__(
        self,
        universe: mda.Universe,
        config: MolClsConfig,
        analyses: Iterable[FrameAnalysis] = (),
    ) -> None:
        """Initialize the MolClusters object.

        The residues must be numbered 1 to N in topology order; otherwise a
        ValueError says so (see `conntable.check_resids`).

        Parameters
        ----------
        universe : mda.Universe
            The MDAnalysis Universe object associated with the simulation.
        config : MolClsConfig
            The analysis configuration.
        analyses : Iterable[FrameAnalysis]
            Analyses to run in addition to the ones the config enables, after them
            and in the order given.

        Raises
        ------
        TypeError
            If one of `analyses` isn't a `FrameAnalysis` instance.
        """
        check_resids(universe)
        extra = list(analyses)
        for analysis in extra:
            if not isinstance(analysis, FrameAnalysis):
                hint = (
                    " (pass an instance, not the class)"
                    if isinstance(analysis, type)
                    else ""
                )
                raise TypeError(
                    f"Analyses must be FrameAnalysis instances, got {analysis!r}{hint}."
                )

        if config.solvent is None and config.solute is not None:
            # on a copy: the caller's config may be reused with another topology
            solvent = sorted(set(universe.residues.resnames) - set(config.solute))
            config = config.model_copy(update={"solvent": solvent})
            logger.info(
                "No 'solvent' configured: using every non-solute residue name in the "
                f"topology: {config.solvent}"
            )

        self.uni = universe
        self.config = config

        self.__check_resnames()
        logger.info(f"Effective configuration:\n{config.describe()}")

        self.tracker: ClusterTracker | None = None  # created by run()

        self.analyses: list[FrameAnalysis] = [*build_builtins(config), *extra]
        if not self.analyses:
            logger.warning(
                "Every analysis is turned off ('analyses' in the config): the clusters "
                "will be tracked, but no results written."
            )

    def analysis[T: FrameAnalysis](self, kind: type[T]) -> T | None:
        """Find the first of the analyses of type `kind`, for its results.

        Parameters
        ----------
        kind : type[T]
            The analysis class to look for (subclasses match too).

        Returns
        -------
        T | None
            The analysis, or None if there's none (e.g. the config doesn't enable
            that built-in, or turns it off).
        """
        for analysis in self.analyses:
            if isinstance(analysis, kind):
                return analysis
        return None

    def __check_resnames(self) -> None:
        """Warn about residue names in the config that the topology doesn't have.

        A misspelled or missing name selects no atoms, so everything built on it
        would silently come out empty.
        """
        present = set(self.uni.residues.resnames)
        comps = self.config.ignore_composition or []
        fields = {
            "rules": self.config._rules.all_keys(),
            "solute": self.config.solute or [],
            "solvent": self.config.solvent or [],
            "nucleus": self.config.nucleus or [],
            "ignore_composition": {name for comp in comps for name in comp},
        }
        for field, names in fields.items():
            missing = sorted(set(names) - present)
            if missing:
                logger.warning(
                    f"'{field}' names residue(s) {missing} that are not in the "
                    f"topology (residue names found: {sorted(present)})."
                )

    def __check_previous_outputs(self, directory: path.Path) -> None:
        """Warn about output files left in the output directory by an earlier run.

        The files appended to (e.g. the per-frame .gro files) would end up mixing an
        earlier run's results with this one's; the other outputs are just overwritten.

        Parameters
        ----------
        directory : path.Path
            The directory this run writes to.
        """
        outputs = self.__declared_outputs()

        overwritten = [
            out.name
            for out in outputs
            if not out.append and (directory / out.name).exists()
        ]
        if overwritten:
            logger.info(f"Overwriting results of an earlier run: {overwritten}")

        appended = sorted(
            {
                file.relative_to(directory).as_posix()
                for out in outputs
                if out.append
                for file in out.existing(directory)
            }
        )
        if appended:
            shown = ", ".join(appended[:5]) + (", ..." if len(appended) > 5 else "")
            logger.warning(
                f"{len(appended)} file(s) from an earlier run are in {directory} "
                f"({shown}): this run appends to them, mixing both runs. "
                "Move or delete them first to keep the runs apart."
            )

    def __declared_outputs(self) -> list[OutputFile]:
        """List the files this run writes, in the order they're reported.

        Returns
        -------
        list[OutputFile]
            The analyses' outputs, in the order of the analyses.
        """
        return [out for analysis in self.analyses for out in analysis.outputs]

    def __log_rule_connections(self, n_frames: int) -> None:
        """Log how many connections each rule found, warning about unused rules.

        Parameters
        ----------
        n_frames : int
            The number of frames analysed.
        """
        rules = self.config._rules
        counts = self.tracker.conntab.rule_connections
        present = set(self.uni.residues.resnames)

        per_frame = ", ".join(
            f"{mi} - {mj} ({rules[mi, mj]}) {counts[mi, mj] / n_frames:.1f}"
            for mi, mj in sorted(rules)
        )
        logger.info(f"Connections per frame, by rule: {per_frame}")

        # rules naming residues missing from the topology were warned about already
        unused = [
            f"{mi} - {mj} ({rules[mi, mj]})"
            for mi, mj in sorted(rules)
            if counts[mi, mj] == 0 and {mi, mj} <= present
        ]
        if unused:
            logger.warning(
                f"Rule(s) {', '.join(unused)} never connected any molecules in "
                f"{n_frames} frame(s): check their cutoffs (distances in angstrom, "
                "angles in degrees)."
            )

    # TODO: break into single_step function to better use in MDRHConstant
    def run(self, output_dir: path.Path | str | None = None) -> None:
        """Track the clusters over the whole trajectory and run the analyses on them.

        Every run starts over from the first frame with a new tracker (cluster ids
        start from 1 again), and the analyses start over in `prepare`, so running
        again gives the same results.

        Parameters
        ----------
        output_dir : path.Path | str | None
            The directory to write the results to, created if missing; the current
            directory when None.
        """
        directory = (
            path.Path.cwd() if output_dir is None else path.Path(output_dir).absolute()
        )
        directory.mkdir(parents=True, exist_ok=True)

        n_frames = len(self.uni.trajectory)
        self.__check_previous_outputs(directory)
        logger.info(f"Tracking clusters over {n_frames} frame(s)")
        start = time.perf_counter()
        # the bar only shows on the terminal, so the log files get a line every 10%
        progress_step = max(1, n_frames // 10)

        self.uni.trajectory[0]
        self.tracker = ClusterTracker(self.uni, self.config)
        tracking = time.perf_counter() - start  # the analyses' time is kept by `run`

        self.output = RunOutput(directory)
        run = Run(self.uni, self.config, n_frames, self.output, self.analyses)

        # the appended files are buffered, so write what was already rendered even
        # if the run is interrupted, as appending it frame by frame used to
        with self.output:
            run.prepare_analyses()
            run.analyse_frame(Frame(0, self.tracker, self.output))

            with tqdm(total=n_frames, initial=1, mininterval=5, miniters=10) as pbar:
                for i, _ in enumerate(self.uni.trajectory[1:], start=1):
                    update_start = time.perf_counter()
                    self.tracker.update()
                    tracking += time.perf_counter() - update_start
                    run.analyse_frame(Frame(i, self.tracker, self.output))
                    pbar.update()

                    done = i + 1
                    if done % progress_step == 0 and done < n_frames:
                        logger.bind(**FILE_ONLY).info(
                            f"Frame {done}/{n_frames} ({done / n_frames:.0%}), "
                            f"{format_duration(time.perf_counter() - start)} elapsed"
                        )

        elapsed = time.perf_counter() - start
        logger.info(
            f"Tracked {n_frames} frame(s) in {format_duration(elapsed)} "
            f"({n_frames / max(elapsed, 1e-9):.1f} frames/s)"
        )
        self.__log_rule_connections(n_frames)

        run.finish_analyses()
        self.output.flush()  # in case an analysis appended to a file in `finish`

        durations = ", ".join(
            f"{type(analysis).__name__} {format_duration(seconds)}"
            for analysis, seconds in zip(self.analyses, run.durations, strict=True)
        )
        logger.info(
            f"Time spent tracking the clusters: {format_duration(tracking)}; "
            f"in each analysis: {durations}"
        )

        outputs = ", ".join(out.name for out in self.__declared_outputs())
        logger.info(f"Results written to {self.output.directory}: {outputs}")
