# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides `JsonReport`, every cluster of every frame and its properties (molclusters.jsonl)."""

import operator
import pathlib as path
from functools import reduce
from typing import Any

from ..config import ReportCompression
from ..output import OutputFile
from ..report import (
    CONNECTION_COLUMNS,
    DECIMALS,
    FORMAT,
    UNITS,
    VERSION,
    cluster_record,
    dumps,
    group_record,
    report_name,
)
from ..version import __version__
from .base import Frame, FrameAnalysis, Run
from .nucleus import Nucleus


def _source(obj: object) -> str | None:
    """Give the absolute path of the file a Universe or trajectory was read from.

    Returns
    -------
    str | None
        The path, or None for one built in memory (e.g. by ``mda.Merge``).
    """
    filename = getattr(obj, "filename", None)
    return str(path.Path(filename).absolute()) if filename else None


class JsonReport(FrameAnalysis):
    """Records every cluster of every frame, with its properties and nuclei.

    The report is written frame by frame, as JSON Lines, compressed as asked; its
    format is described in `molclusters.report`, which also reads it back. The
    nuclei come from the `Nucleus` analysis, when it runs before this one.

    Attributes
    ----------
    compression : ReportCompression
        How the report is compressed: "zstd", "gzip" or "none".
    name : str
        The report's file name, e.g. ``molclusters.jsonl.zst``.
    """

    def __init__(self, compression: ReportCompression = "zstd") -> None:
        """Initialize the report, started by `prepare`.

        Parameters
        ----------
        compression : ReportCompression
            How to compress the report: "zstd", "gzip" or "none".
        """
        self.compression = compression
        self.name = report_name(compression)
        self.outputs = (OutputFile(self.name),)
        self._nucleus: Nucleus | None = None

    def prepare(self, run: Run) -> None:
        """Start the report with its header: the run's inputs and configuration.

        Parameters
        ----------
        run : Run
            The run about to start.
        """
        self._nucleus = run.analysis(Nucleus)
        header = {
            "format": FORMAT,
            "version": VERSION,
            "software": f"MolClusters {__version__}",
            "trajectory": _source(run.universe.trajectory),
            "topology": _source(run.universe),
            "n_frames": run.n_frames,
            "decimals": DECIMALS,
            "units": UNITS,
            "connection_columns": CONNECTION_COLUMNS,
            "config": run.config.model_dump(mode="json"),
        }
        run.output.append(self.name, dumps(header))

    def analyse(self, frame: Frame) -> None:
        """Add the clusters of the current frame to the report.

        Parameters
        ----------
        frame : Frame
            The current frame.
        """
        clusters = []
        for cls in frame.clusters.values():
            record = cluster_record(cls)
            if self._nucleus is not None:
                record["Nucleus"] = self._nuclei(cls.id)
            clusters.append(record)

        record = {
            "frame": int(frame.universe.coord.frame),
            "time": round(float(frame.time), DECIMALS),
            "clusters": clusters,
        }
        frame.output.append(self.name, dumps(record))

    def _nuclei(self, cluster_id: int) -> dict[str, Any]:
        """Encode the nuclei of a cluster.

        Parameters
        ----------
        cluster_id : int
            The cluster's id.

        Returns
        -------
        dict[str, Any]
            ``nuclei``, one `group_record` per nucleus, and
            ``combined_dipole_moment``, the dipole moment of all of them together
            (None without nuclei).
        """
        nuclei = self._nucleus.nuclei.get(cluster_id, [])
        combined = (
            round(float(reduce(operator.add, nuclei).dipole_moment), DECIMALS)
            if nuclei
            else None
        )
        return {
            "nuclei": [group_record(nucleus) for nucleus in nuclei],
            "combined_dipole_moment": combined,
        }
