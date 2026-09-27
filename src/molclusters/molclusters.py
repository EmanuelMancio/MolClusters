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
import numpy as np
from loguru import logger
from tqdm import tqdm

from .analysis import (
    ClusterCoordinates,
    Frame,
    FrameAnalysis,
    JsonReport,
    Nucleus,
    Run,
    SizeEvolution,
    SoluteSolvent,
)
from .cluster import Cluster
from .config import MolClsConfig
from .log import FILE_ONLY, format_duration
from .output import OutputFile, RunOutput
from .tracker import ClusterTracker


class MolClusters:
    """A class for analyzing molecular clusters in molecular dynamics simulations.

    It tracks the clusters frame by frame and runs the analyses on them: the
    built-ins the config enables (cluster sizes, solute-solvent, coordinates,
    nuclei, the JSON report), then the ones given to it.

    Attributes
    ----------
    uni : mda.Universe
        The MDAnalysis Universe object associated with the simulation.
    config : dict
        The configuration dictionary containing analysis settings.
    tracker : ClusterTracker
        Follows the clusters frame by frame, keeping their ids stable.
    clusters : dict[int, Cluster]
        The tracker's clusters of the current frame, keyed by cluster ID.
    mol_clt : dict[int, int]
        The tracker's mapping of molecule IDs to their respective cluster IDs.
    analyses : list[FrameAnalysis]
        The analyses run on the clusters of every frame.
    size_evolution : SizeEvolution
        The number and sizes of the clusters over time (one of `analyses`).
    nucleus : Nucleus | None
        The nuclei inside the clusters, when the config enables them (one of
        `analyses`).
    clusters_size_evo : np.ndarray
        An array tracking the evolution of cluster sizes over time (see
        `size_evolution`).
    """

    __slots__ = [
        "uni",
        "config",
        "tracker",
        "analyses",
        "size_evolution",
        "radius_evolution",
        "nucleus",
        "output",
    ]

    def __init__(
        self,
        universe: mda.Universe,
        config: MolClsConfig,
        analyses: Iterable[FrameAnalysis] = (),
    ) -> None:
        """Initialize the MolClusters object.

        Parameters
        ----------
        universe : mda.Universe
            The MDAnalysis Universe object associated with the simulation.
        config : dict
            The configuration dictionary containing analysis settings.
        analyses : Iterable[FrameAnalysis]
            Analyses to run in addition to the ones the config enables, after them
            and in the order given.

        Raises
        ------
        TypeError
            If one of `analyses` isn't a `FrameAnalysis` instance.
        """
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

        self.uni = universe
        self.config = config

        if config.solvent is None and config.solute is not None:
            config.solvent = sorted(
                set(universe.residues.resnames) - set(config.solute)
            )
            logger.info(
                "No 'solvent' configured: using every non-solute residue name in the "
                f"topology: {config.solvent}"
            )

        self.__check_resnames()
        logger.info(f"Effective configuration:\n{config.describe()}")

        # TODO: move start to run
        self.tracker = ClusterTracker(self.uni, self.config)

        self.size_evolution = SizeEvolution()
        builtins: list[FrameAnalysis] = [self.size_evolution]
        if config.solute is not None:
            builtins += [
                SoluteSolvent(config.solute, config.solvent),
                ClusterCoordinates(config.solute, follow=config._follow_solute),
            ]
        self.nucleus = None if config.nucleus is None else Nucleus(config.nucleus)
        if self.nucleus is not None:
            builtins.append(self.nucleus)
        builtins.append(JsonReport())
        self.analyses: list[FrameAnalysis] = [*builtins, *extra]

        self.radius_evolution = {}

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

    def __check_previous_outputs(self) -> None:
        """Warn about output files left in the working directory by an earlier run.

        The files appended to (e.g. the per-frame .gro files) would end up mixing an
        earlier run's results with this one's; the other outputs are just overwritten.
        """
        cwd = path.Path.cwd()
        outputs = self.__declared_outputs()

        overwritten = [
            out.name for out in outputs if not out.append and (cwd / out.name).exists()
        ]
        if overwritten:
            logger.info(f"Overwriting results of an earlier run: {overwritten}")

        appended = sorted(
            file.name
            for file in cwd.iterdir()
            if any(out.append and out.matches(file.name) for out in outputs)
        )
        if appended:
            shown = ", ".join(appended[:5]) + (", ..." if len(appended) > 5 else "")
            logger.warning(
                f"{len(appended)} file(s) from an earlier run are in {cwd} "
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
    def run(self) -> None:
        """Run the molecular cluster analysis.

        This method performs cluster detection, solute-solvent analysis, nucleus analysis,
        and exports the results to files.
        """
        n_frames = len(self.uni.trajectory)
        self.__check_previous_outputs()
        logger.info(f"Tracking clusters over {n_frames} frame(s)")
        start = time.perf_counter()
        # the bar only shows on the terminal, so the log files get a line every 10%
        progress_step = max(1, n_frames // 10)

        self.output = RunOutput(path.Path.cwd())
        run = Run(self.uni, self.config, n_frames, self.output, self.analyses)

        # the appended files are buffered, so write what was already rendered even
        # if the run is interrupted, as appending it frame by frame used to
        with self.output:
            run.prepare_analyses()
            run.analyse_frame(Frame(0, self.tracker))

            with tqdm(total=n_frames, initial=1, mininterval=5, miniters=10) as pbar:
                for i, _ in enumerate(self.uni.trajectory[1:], start=1):
                    self.tracker.update()
                    run.analyse_frame(Frame(i, self.tracker))
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

        outputs = ", ".join(out.name for out in self.__declared_outputs())
        logger.info(f"Results written to {self.output.directory}: {outputs}")

    @property
    def clusters(self) -> dict[int, Cluster]:
        """The clusters of the current frame, keyed by cluster ID (see `tracker`)."""
        return self.tracker.clusters

    @property
    def mol_clt(self) -> dict[int, int]:
        """Molecule ID -> cluster ID for the current frame (see `tracker`)."""
        return self.tracker.mol_clt

    def find(self, mol: int) -> int | bool:
        """Find the cluster ID for a given molecule.

        Parameters
        ----------
        mol : int
            The molecule ID to search for.

        Returns
        -------
        int | bool
            The cluster ID if the molecule is found, or False if not found.
        """
        return self.tracker.find(mol)

    @property
    def clusters_size_evo(self) -> np.ndarray:
        """The number and sizes of the clusters, frame by frame (see `size_evolution`)."""
        return self.size_evolution.data

    def __print_clusters_index(self) -> None:
        """Write ndx file with cluster index."""
        cols = 15
        i = 1
        with path.Path("clusters_index.ndx").open("w+", encoding="utf-8") as ndx:
            for id, cluster in self.tracker.clusters.items():
                ndx.write(f"[ CLS-{id} ]\n")
                for mol in sorted(cluster):
                    for at in self.uni.residues[mol - 1].atoms:
                        if i % cols != 0:
                            ndx.write(f"{at.id + 1}\t")
                        else:
                            ndx.write(f"{at.id + 1}\n")

                        i += 1

                ndx.write("\n")

                # TODO: (low priority) Make the skipped lines work
