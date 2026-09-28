# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides `JsonReport`, every cluster of every frame and its properties (molclusters.jsonl)."""

import pathlib as path
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
    report_name,
    rounded,
)
from ..version import __version__
from .base import Frame, FrameAnalysis, Run


def _source(obj: object) -> str | None:
    """Give the absolute path of the file a Universe or trajectory was read from.

    Returns
    -------
    str | None
        The path, or None for one built in memory (e.g. by ``mda.Merge``).
    """
    filename = getattr(obj, "filename", None)
    return str(path.Path(filename).absolute()) if filename else None


def _overrides(analysis: FrameAnalysis, hook: str) -> bool:
    """Tell whether an analysis overrides one of `FrameAnalysis`'s report hooks.

    Returns
    -------
    bool
        True if `hook` is the analysis' own, not the base class' (which adds
        nothing).
    """
    return getattr(type(analysis), hook) is not getattr(FrameAnalysis, hook)


def _contribution(
    analysis: FrameAnalysis, hook: str, *args: object
) -> dict[str, Any] | None:
    """Get what an analysis adds to a record, rounded as the report's floats are.

    Parameters
    ----------
    analysis : FrameAnalysis
        The analysis.
    hook : str
        Its hook to call: "report_frame" or "report_cluster".
    *args : object
        The hook's arguments.

    Returns
    -------
    dict[str, Any] | None
        The fields, or None if the analysis adds nothing.
    """
    try:
        fields = getattr(analysis, hook)(*args)
        return None if fields is None else rounded(fields)
    except Exception as err:
        # the run notes JsonReport's hook, not the one it called
        err.add_note(f"In {type(analysis).__name__}.{hook}(), for the report")
        raise


class JsonReport(FrameAnalysis):
    """Records every cluster of every frame, with its properties.

    The report is written frame by frame, as JSON Lines, compressed as asked; its
    format is described in `molclusters.report`, which also reads it back. Every
    other analysis of the run can add fields to it (see
    `FrameAnalysis.report_frame` and `report_cluster`), so it runs after them all.

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
        self._frame_contributors: list[tuple[str, FrameAnalysis]] = []
        self._cluster_contributors: list[tuple[str, FrameAnalysis]] = []

    def prepare(self, run: Run) -> None:
        """Start the report with its header: the run's inputs and configuration.

        Parameters
        ----------
        run : Run
            The run about to start.

        Raises
        ------
        ValueError
            If two analyses of the same class name would add to the report.
        """
        contributors = [
            analysis
            for analysis in run.analyses
            if _overrides(analysis, "report_frame")
            or _overrides(analysis, "report_cluster")
        ]
        names = [type(analysis).__name__ for analysis in contributors]
        repeated = sorted({name for name in names if names.count(name) > 1})
        if repeated:
            # their fields would go under the same key
            raise ValueError(
                f"More than one {' and '.join(repeated)} analysis would add to the "
                "report: run only one of each."
            )

        self._frame_contributors = [
            (name, analysis)
            for name, analysis in zip(names, contributors, strict=True)
            if _overrides(analysis, "report_frame")
        ]
        self._cluster_contributors = [
            (name, analysis)
            for name, analysis in zip(names, contributors, strict=True)
            if _overrides(analysis, "report_cluster")
        ]

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
            # this one runs last, so the run's analyses are those before it, and it
            "analyses": [type(a).__name__ for a in (*run.analyses, self)],
            "contributors": names,
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
            for name, analysis in self._cluster_contributors:
                fields = _contribution(analysis, "report_cluster", frame, cls)
                if fields is not None:
                    record[name] = fields
            clusters.append(record)

        record = {
            "frame": int(frame.universe.coord.frame),
            "time": round(float(frame.time), DECIMALS),
            "clusters": clusters,
        }
        for name, analysis in self._frame_contributors:
            fields = _contribution(analysis, "report_frame", frame)
            if fields is not None:
                record[name] = fields
        frame.output.append(self.name, dumps(record))
