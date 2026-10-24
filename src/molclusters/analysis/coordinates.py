# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Provides `ClusterCoordinates`, the coordinates of the clusters holding solutes (.gro files)."""

from collections import Counter
from collections.abc import Iterable
from itertools import chain
from pathlib import PurePosixPath

import numpy as np
from loguru import logger
from MDAnalysis import core, units
from MDAnalysis.exceptions import NoDataError

from ..cluster import Cluster
from ..output import OutputFile
from .base import Frame, FrameAnalysis, Run

# MDAnalysis' GRO writer's formats (`GROWriter.fmt`): an atom, its velocity, and
# the box, rectangular or not (the triclinic vectors' components, in the writer's
# order)
_GRO_ATOM = "{:>5d}{:<5.5s}{:>5.5s}{:>5d}{:8.3f}{:8.3f}{:8.3f}"
_GRO_VELOCITY = "{:8.4f}{:8.4f}{:8.4f}"
_GRO_BOX = "{:10.5f} {:9.5f} {:9.5f}\n"
_GRO_TRICLINIC_BOX = (
    "{0:10.5f} {4:9.5f} {8:9.5f} {1:9.5f} {2:9.5f} {3:9.5f} {5:9.5f} {6:9.5f} "
    "{7:9.5f}\n"
)
# .gro files are in nm and nm/ps, converted as the writer does
_NM = units.get_conversion_factor("length", "Angstrom", "nm")
_NM_PS = units.get_conversion_factor("speed", "Angstrom/ps", "nm/ps")
# the coordinates a .gro file can hold, in nm (`GROWriter.gro_coor_limits`)
_GRO_MIN, _GRO_MAX = -999.9995, 9999.9995


def gro_frame(atoms: core.groups.AtomGroup, title: str) -> str:
    """Render atoms as one frame of a .gro file, as MDAnalysis' GRO writer does.

    The writer only writes whole files, and formats every atom through calls of
    its own; this formats them all at once, several times faster, into the same
    text. The atoms are numbered from 1 in the order given, and their resids and
    numbers keep their last 5 digits, as the format has room for. Atoms without
    names are written as "X", as by the writer (which also warns about it).

    Parameters
    ----------
    atoms : core.groups.AtomGroup
        The atoms, at their current positions.
    title : str
        The frame's title line.

    Returns
    -------
    str
        The frame's text.

    Raises
    ------
    ValueError
        If a coordinate is out of the range a .gro file can hold.
    """
    n_atoms = len(atoms)
    positions = _NM * atoms.positions
    if not ((_GRO_MIN < positions).all() and (positions <= _GRO_MAX).all()):
        raise ValueError(
            f"GRO files must have coordinate values between {_GRO_MIN:.3f} and "
            f"{_GRO_MAX:.3f} nm, so {title!r} can't be written."
        )
    try:
        names = atoms.names
    except NoDataError:
        names = np.full(n_atoms, "X")
    columns = [
        atoms.resids % 100_000,
        atoms.resnames,
        names,
        np.arange(1, n_atoms + 1) % 100_000,
        *positions.T,
    ]
    line = _GRO_ATOM
    try:
        columns += [*(_NM_PS * atoms.velocities).T]
        line += _GRO_VELOCITY
    except NoDataError:
        pass
    # one format call for every atom: the columns, row by row
    rows = zip(*(column.tolist() for column in columns), strict=True)
    atom_lines = ((line + "\n") * n_atoms).format(*chain.from_iterable(rows))

    box = atoms.dimensions
    if box is None:
        footer = _GRO_BOX.format(0.0, 0.0, 0.0)
    elif np.allclose(box[3:], 90.0):
        footer = _GRO_BOX.format(*(_NM * box[:3]).tolist())
    else:
        vectors = _NM * atoms.universe.coord.triclinic_dimensions.flatten()
        footer = _GRO_TRICLINIC_BOX.format(*vectors.tolist())

    return f"{title}\n{n_atoms:5d}\n{atom_lines}{footer}"


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
            OutputFile(self._file("cls-n<size>.gro")),
            OutputFile(self._file("cls-id<id>.gro")),
        )
        if follow:
            self.outputs += (OutputFile(self._file("solute-<resid>.gro")),)

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
        with cls.whole() as atoms:
            return gro_frame(atoms.sort(), f"Cluster-{cls.id} - Time = {time}")

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
