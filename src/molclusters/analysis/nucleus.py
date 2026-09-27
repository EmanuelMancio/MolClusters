# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides `Nucleus`, the nuclei inside the clusters (nucleus_data.csv)."""

from collections.abc import Iterable

import networkx as nx
import numpy as np
import pandas as pd

from ..cluster import MolGroup
from ..output import OutputFile
from .base import Frame, FrameAnalysis, Run

COLUMNS = [
    "Time",
    "NNuc",
    "Size",
    "Radius",
    "Density",
    "Charge",
    "Dipole",
    "Spher",
    "Shape",
]


class Nucleus(FrameAnalysis):
    """Finds the nuclei inside the clusters, frame by frame.

    A nucleus is a connected group, inside a cluster, of molecules with one of the
    nucleus residue names.

    Attributes
    ----------
    resnames : set[str]
        The residue names of the molecules that make up nuclei.
    nuclei : dict[int, list[MolGroup]]
        Cluster id -> the cluster's nuclei in the current frame (an empty list for
        a cluster without nuclei), for the analyses after this one.
    data : np.ndarray
        One row per frame (see `COLUMNS`): time, the average number of nuclei of
        the clusters that have any (0 without nuclei), and the nuclei's average
        size, radius, density, charge, dipole moment, sphericity and shape
        parameter (NaN without nuclei).
    """

    outputs = (OutputFile("nucleus_data.csv"),)

    def __init__(self, resnames: Iterable[str]) -> None:
        """Initialize an empty record, sized by `prepare`.

        Parameters
        ----------
        resnames : Iterable[str]
            The residue names of the molecules that make up nuclei.
        """
        self.resnames = set(resnames)
        self.nuclei: dict[int, list[MolGroup]] = {}
        self.data = np.full((0, len(COLUMNS)), np.nan)

    def prepare(self, run: Run) -> None:
        """Make room for every frame of the run.

        Parameters
        ----------
        run : Run
            The run about to start.
        """
        self.data = np.full((run.n_frames, len(COLUMNS)), np.nan)

    def analyse(self, frame: Frame) -> None:
        """Find the nuclei of the current frame and record their averages.

        Parameters
        ----------
        frame : Frame
            The current frame.
        """
        self.nuclei = {}
        n_nucleus = []
        sizes = []
        radius = []
        dipole = []
        density = []
        sphericity = []
        shape = []
        charge = []
        for cid, cls in frame.clusters.items():
            possible_nucleus = [
                rid
                for rnm, rid in zip(cls.resnames, cls.resids, strict=True)
                if rnm in self.resnames
            ]

            self.nuclei[cid] = []
            subcomps = nx.induced_subgraph(cls.graph, possible_nucleus)
            n_nuc = 0
            for sg in nx.connected_components(subcomps):
                tp = MolGroup(frame.universe, list(sg))

                n_nuc += 1
                sizes.append(tp.size)
                radius.append(tp.radius)
                dipole.append(tp.dipole_moment)
                density.append(tp.density)
                sphericity.append(tp.sphericity)
                shape.append(tp.shape_parameter)
                charge.append(tp.charge)

                self.nuclei[cid].append(tp)

            if n_nuc != 0:
                n_nucleus.append(n_nuc)

        row = self.data[frame.index]
        row[0] = frame.time
        row[1] = 0 if len(n_nucleus) == 0 else np.average(n_nucleus)
        row[2] = np.nan if len(n_nucleus) == 0 else np.average(sizes)
        row[3] = np.nan if len(n_nucleus) == 0 else np.average(radius)
        row[4] = np.nan if len(n_nucleus) == 0 else np.average(density)
        row[5] = np.nan if len(n_nucleus) == 0 else np.average(charge)
        row[6] = np.nan if len(n_nucleus) == 0 else np.average(dipole)
        row[7] = np.nan if len(n_nucleus) == 0 else np.average(sphericity)
        row[8] = np.nan if len(n_nucleus) == 0 else np.average(shape)

    def finish(self, run: Run) -> None:
        """Write the record to nucleus_data.csv.

        Parameters
        ----------
        run : Run
            The run that just ended.
        """
        pd.DataFrame(self.data, columns=COLUMNS).to_csv(
            run.output.path("nucleus_data.csv"), index=False
        )
