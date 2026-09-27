# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides `ClusterCoordinates`, the coordinates of the clusters holding solutes (.gro files)."""

import io
from collections import Counter
from collections.abc import Iterable
from pathlib import PurePosixPath

import MDAnalysis as mda
from loguru import logger
from MDAnalysis.lib.util import NamedStream

from ..cluster import Cluster
from ..output import OutputFile
from .base import Frame, FrameAnalysis, Run


class ClusterCoordinates(FrameAnalysis):
    """Writes the coordinates of every cluster holding a solute, frame by frame.

    Each such cluster, made whole, is appended to ``cls-n<size>.gro`` (the clusters
    of one size, pooled) and ``cls-id<id>.gro`` (one cluster's own trajectory). When
    following the solutes, a cluster holding exactly one solute is also appended to
    ``solute-<resid>.gro``. The files go in the ``folder`` of the output directory.

    Attributes
    ----------
    solutes : set[str]
        The residue names of the solutes.
    follow : bool
        Whether to write each solute's cluster to its own file.
    folder : str
        The folder of the output directory the files go in.
    solute_ids : set[int]
        The resids of the solute molecules, found by `prepare`.
    follow_skipped : Counter[int]
        Cluster id -> number of frames it wasn't followed in, for holding more than
        one solute.
    """

    def __init__(
        self,
        solutes: Iterable[str],
        *,
        follow: bool = False,
        folder: str = "coordinates",
    ) -> None:
        """Initialize the analysis.

        Parameters
        ----------
        solutes : Iterable[str]
            The residue names of the solutes.
        follow : bool
            Whether to write each solute's cluster to its own file.
        folder : str
            The folder of the output directory the files go in ("" for the
            output directory itself).
        """
        self.solutes = set(solutes)
        self.follow = follow
        self.folder = folder
        self.solute_ids: set[int] = set()
        self.follow_skipped: Counter[int] = Counter()

        self.outputs = (
            OutputFile(self._file("cls-n<size>.gro"), append=True),
            OutputFile(self._file("cls-id<id>.gro"), append=True),
        )
        if follow:
            self.outputs += (OutputFile(self._file("solute-<resid>.gro"), append=True),)

    def _file(self, name: str) -> str:
        """Place the file `name` in `folder`.

        Parameters
        ----------
        name : str
            The file name.

        Returns
        -------
        str
            The file's path relative to the output directory.
        """
        return str(PurePosixPath(self.folder, name))

    def prepare(self, run: Run) -> None:
        """Find the solute molecules.

        Parameters
        ----------
        run : Run
            The run about to start.
        """
        self.solute_ids = {
            int(resid)
            for name in sorted(self.solutes)
            for resid in run.universe.select_atoms(f"resname {name}").residues.resids
        }
        self.follow_skipped = Counter()

    def analyse(self, frame: Frame) -> None:
        """Append the clusters of the current frame that hold a solute.

        Parameters
        ----------
        frame : Frame
            The current frame.
        """
        for cls in frame.clusters.values():
            sol_ids = self.solute_ids.intersection(cls.resids)
            if not sol_ids:
                continue

            text = self._render(cls, frame.time)

            # pooled by size: an ensemble of what an N-mer looks like, across all
            # clusters that were ever that size, independent of cluster identity
            frame.output.append(self._file(f"cls-n{cls.size}.gro"), text)

            # pooled by identity: this specific cluster's own trajectory, tracked
            # across frames via the dominance algorithm regardless of size changes
            frame.output.append(self._file(f"cls-id{cls.id}.gro"), text)

            # TODO: change to support merges
            # FIXME: with changes in config this needs to be updated
            if self.follow:
                if len(sol_ids) > 1:
                    # TODO: make more feature-rich follow procedure
                    # counted, not logged: `finish` sums them up
                    self.follow_skipped[cls.id] += 1
                    continue

                (sol_id,) = sol_ids

                frame.output.append(self._file(f"solute-{sol_id}.gro"), text)

    @staticmethod
    def _render(cls: Cluster, time: float) -> str:
        """Render a cluster, made whole, as one frame of a .gro file.

        Parameters
        ----------
        cls : Cluster
            The cluster to render.
        time : float
            The time of the current frame, for the title line.

        Returns
        -------
        str
            The frame's text, title line included.
        """
        # the GRO writer only writes whole files, so render the frame in memory
        # (NamedStream keeps the buffer open when the writer closes it)
        buf = io.StringIO()
        with (
            mda.Writer(NamedStream(buf, "cluster.gro"), multiframe=False) as w,
            cls.whole() as atoms,
        ):
            w.write(atoms.sort())

        # replace the writer's fixed "Written by MDAnalysis" title line
        _, body = buf.getvalue().split("\n", 1)
        return f"Cluster-{cls.id} - Time = {time}\n{body}"

    def finish(self, run: Run) -> None:
        """Warn about the frames in which solutes weren't followed.

        Parameters
        ----------
        run : Run
            The run that just ended.
        """
        if self.follow_skipped:
            logger.warning(
                f"Solutes were not followed in {self.follow_skipped.total()} "
                f"frame(s) of {len(self.follow_skipped)} cluster(s) holding more than "
                "one solute: those frames are missing from the solute-<resid>.gro "
                "files. The cluster ids are logged at DEBUG level (--log-level DEBUG)."
            )
            logger.debug(f"Frames not followed, by cluster id: {self.follow_skipped}")
