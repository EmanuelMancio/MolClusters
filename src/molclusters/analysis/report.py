# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides `JsonReport`, every cluster of every frame and its properties (molclusters.json)."""

import json
import pathlib as path
from functools import reduce
from typing import Any

from ..cluster import Cluster, MolGroup
from ..output import OutputFile
from ..version import __version__
from .base import Frame, FrameAnalysis, Run
from .nucleus import Nucleus


# TODO: use orjson for better encoding options
class JsonReport(FrameAnalysis):
    """Records every cluster of every frame, with its properties and nuclei.

    The nuclei come from the `Nucleus` analysis, when it runs before this one.

    Attributes
    ----------
    data : dict
        The report: the software, input paths and configuration of the run, and
        under ``"MolClusters"`` one entry per frame with its clusters.
    """

    outputs = (OutputFile("molclusters.json"),)

    def __init__(self) -> None:
        """Initialize an empty report, started by `prepare`."""
        self.data: dict[str, Any] = {}
        self._nucleus: Nucleus | None = None

    def prepare(self, run: Run) -> None:
        """Start the report with the run's inputs and configuration.

        Parameters
        ----------
        run : Run
            The run about to start.
        """
        self._nucleus = run.analysis(Nucleus)
        self.data = {
            "Software": f"MolClusters {__version__}",
            "Trajectory": str(path.Path(run.universe.trajectory.filename).absolute()),
            "Topology": str(path.Path(run.universe.filename).absolute()),
            "Config": run.config.model_dump(),
            "MolClusters": [],
        }

    def analyse(self, frame: Frame) -> None:
        """Add the clusters of the current frame to the report.

        Parameters
        ----------
        frame : Frame
            The current frame.
        """
        data = {}

        data["Time"] = frame.time
        data["Frame"] = frame.universe.coord.frame
        data["NClusters"] = len(frame.clusters)

        molclusters_data = []
        for cid, cls in frame.clusters.items():
            cls_data = self.encode_cluster(cls)
            if self._nucleus is not None:
                cls_data["Nucleus"] = []
                if self._nucleus.nuclei.get(cid, False):
                    nuclei = reduce(lambda a, b: a + b, self._nucleus.nuclei[cid])
                    cls_data["NucleiDipole"] = nuclei.dipole_moment
                    for nuc in self._nucleus.nuclei[cid]:
                        cls_data["Nucleus"].append(JsonReport.encode_nucleus(nuc))

            molclusters_data.append(cls_data)

        data["Clusters"] = molclusters_data
        self.data["MolClusters"].append(data)

    def finish(self, run: Run) -> None:
        """Write the report to molclusters.json.

        Parameters
        ----------
        run : Run
            The run that just ended.
        """
        with run.output.path("molclusters.json").open("w+") as json_out:
            json.dump(self.data, json_out, indent=2)

    @staticmethod
    def encode_cluster(cls: Cluster) -> dict:
        """Encode a cluster into a dictionary.

        Parameters
        ----------
        cls : Cluster
            The cluster to encode.

        Returns
        -------
        dict
            The encoded cluster data.
        """
        data = {}
        data["ID"] = cls.id
        JsonReport.encode_properties(cls, data)
        return data

    @staticmethod
    def encode_nucleus(nuc: MolGroup) -> dict:
        """Encode a nucleus into a dictionary.

        Parameters
        ----------
        nuc : MolGroup
            The nucleus to encode.

        Returns
        -------
        dict
            The encoded nucleus data.
        """
        data = {}
        JsonReport.encode_properties(nuc, data)
        return data

    @staticmethod
    def encode_properties(obj: Cluster | MolGroup, data: dict) -> None:
        """Encode the properties of a cluster or nucleus.

        Parameters
        ----------
        obj : Cluster | MolGroup
            The object to encode.
        data : dict
            The dictionary to store the encoded properties.
        """
        data["Size"] = obj.size
        data["Composition"] = JsonReport.encode_composition(obj)

        if isinstance(obj, Cluster):
            data["Connections"] = JsonReport.encode_connections(obj)

        data["ResIDs"] = sorted([int(x) for x in obj.resids])
        data["Mass"] = obj.mass
        data["Volume"] = obj.volume
        data["Radius"] = obj.radius
        data["Diameter"] = obj.diameter
        data["Density"] = obj.density
        data["Charge"] = obj.charge
        data["Dipole Moment"] = obj.dipole_moment
        data["Sphericity"] = obj.sphericity
        data["Shape"] = obj.shape_parameter

    @staticmethod
    def encode_composition(obj: Cluster | MolGroup) -> list[dict]:
        """Encode the composition of a cluster or nucleus.

        Parameters
        ----------
        obj : Cluster | MolGroup
            The object to encode.

        Returns
        -------
        list[dict]
            A list of dictionaries representing the composition.
        """
        comp = {}
        for rnm, rid in zip(obj.resnames, map(int, obj.resids), strict=True):
            if rnm in comp:
                comp[rnm]["n"] += 1
                comp[rnm]["resids"].append(rid)
            else:
                comp[rnm] = {"resname": rnm, "n": 1, "resids": [rid]}

        return list(comp.values())

    @staticmethod
    def encode_connections(
        obj: Cluster,
    ) -> list[tuple[int, int, dict[str, Any]]]:
        """Encode the connections of a cluster.

        Parameters
        ----------
        obj : Cluster
            The object to encode.

        Returns
        -------
        list[tuple[int, int, dict[str, Any]]]
            A list of tuples representing the connections and their properties.
        """
        return [
            (
                int(edge[0]),
                int(edge[1]),
                # counts (n_hbonds) stay whole numbers
                {k: v if isinstance(v, int) else float(v) for k, v in edge[2].items()},
            )
            for edge in obj.graph.edges.data()
        ]
