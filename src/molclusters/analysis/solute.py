# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides `SoluteSolvent`, the clusters holding both solutes and solvents (solute_solvent.csv)."""

from collections.abc import Iterable, Iterator

import numpy as np
import pandas as pd

from ..cluster import Cluster
from ..output import OutputFile
from .base import Frame, FrameAnalysis, Run

COLUMNS = [
    "Time",
    "NCls",
    "NSolt",
    "NSolv",
    "Radius",
    "Density",
    "Charge",
    "Dipole",
    "Spher",
    "Shape",
]


class SoluteSolvent(FrameAnalysis):
    """Records the clusters that hold both solute and solvent molecules, frame by frame.

    Attributes
    ----------
    solutes : set[str]
        The residue names of the solutes.
    solvents : set[str]
        The residue names of the solvents.
    data : np.ndarray
        One row per frame (see `COLUMNS`): time, the number of such clusters, their
        average number of solutes and of solvents (0 without such clusters), and
        their average radius, density, charge, dipole moment, sphericity and shape
        parameter (NaN without such clusters).
    """

    outputs = (OutputFile("solute_solvent.csv"),)

    def __init__(self, solutes: Iterable[str], solvents: Iterable[str]) -> None:
        """Initialize an empty record, sized by `prepare`.

        Parameters
        ----------
        solutes : Iterable[str]
            The residue names of the solutes.
        solvents : Iterable[str]
            The residue names of the solvents.
        """
        self.solutes = set(solutes)
        self.solvents = set(solvents)
        self.data = np.zeros((0, len(COLUMNS)))

    def prepare(self, run: Run) -> None:
        """Make room for every frame of the run.

        Parameters
        ----------
        run : Run
            The run about to start.
        """
        self.data = np.zeros((run.n_frames, len(COLUMNS)))

    def _solute_solvent_clusters(self, frame: Frame) -> Iterator[Cluster]:
        """Generate the clusters of the frame holding both solutes and solvents.

        Parameters
        ----------
        frame : Frame
            The current frame.

        Yields
        ------
        Cluster
            Each cluster holding at least one solute and one solvent molecule.
        """
        for cls in frame.clusters.values():
            res = set(cls.resnames)
            if res & self.solutes and res & self.solvents:
                yield cls

    def analyse(self, frame: Frame) -> None:
        """Record the clusters holding both solutes and solvents in the current frame.

        Parameters
        ----------
        frame : Frame
            The current frame.
        """
        n_cls = 0
        n_solvents = []
        n_solutes = []
        radius = []
        dipole = []
        density = []
        sphericity = []
        shape = []
        charge = []

        for cls in self._solute_solvent_clusters(frame):
            n_cls += 1

            composition = cls.composition
            n_solvents.append(sum(composition.get(name, 0) for name in self.solvents))
            n_solutes.append(sum(composition.get(name, 0) for name in self.solutes))

            radius.append(cls.radius_of_gyration)
            dipole.append(cls.dipole_moment)
            density.append(cls.density)
            sphericity.append(cls.sphericity)
            shape.append(cls.shape_parameter)
            charge.append(cls.charge)

        row = self.data[frame.index]
        row[0] = frame.time
        row[1] = n_cls
        row[2] = 0 if n_cls == 0 else np.average(n_solutes)
        row[3] = 0 if n_cls == 0 else np.average(n_solvents)
        row[4] = np.nan if n_cls == 0 else np.average(radius)
        row[5] = np.nan if n_cls == 0 else np.average(density)
        row[6] = np.nan if n_cls == 0 else np.average(charge)
        row[7] = np.nan if n_cls == 0 else np.average(dipole)
        row[8] = np.nan if n_cls == 0 else np.average(sphericity)
        row[9] = np.nan if n_cls == 0 else np.average(shape)

    def finish(self, run: Run) -> None:
        """Write the record to solute_solvent.csv.

        Parameters
        ----------
        run : Run
            The run that just ended.
        """
        pd.DataFrame(self.data, columns=COLUMNS).to_csv(
            run.output.path("solute_solvent.csv"), index=False
        )
