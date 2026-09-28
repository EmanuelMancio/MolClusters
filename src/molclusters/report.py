# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides the format of the JSON report (molclusters.jsonl), and its reader.

The report, written by `analysis.JsonReport`, is JSON Lines: one JSON object per
line, compressed with zstd (``molclusters.jsonl.zst``, the default), gzip
(``.jsonl.gz``) or not at all (``.jsonl``), as the config's `report_compression`
says. It's written frame by frame, so a run never holds it in memory, and an
interrupted run leaves every frame it analysed readable.

The first line is the header::

    {"format": "molclusters-report", "version": 2, "software": "MolClusters x.y.z",
     "trajectory": <absolute path, or null for an in-memory trajectory>,
     "topology": <absolute path, or null>, "n_frames": <frames in the run>,
     "decimals": 6, "units": {<quantity>: <unit>, ...},
     "connection_columns": ["i", "j", "distance", "angle", "n_hbonds"],
     "contributors": [<the analyses that add fields, see below>],
     "config": <the effective configuration>}

Then one line per frame::

    {"frame": <trajectory frame number>, "time": <ps>, "clusters": [<cluster>, ...]}

and each cluster is::

    {
        "id": 3,
        "size": 3,
        "resids": [107, 169, 171],
        "composition": {"MAL": 3},
        "mass": 402.260990,
        "volume": 834.555760,
        "radius": 5.840575,
        "diameter": 11.681149,
        "density": 0.800390,
        "charge": 0.0,
        "dipole_moment": 15.696198,
        "sphericity": 0.738690,
        "shape_parameter": 0.138796,
        "connections": [[169, 107, 2.684737, 171.434846, 1], ...],
    }

with each connection a row of ``connection_columns``: the two resids, then the
connection's distance (and, for H-bonds, the D-H-A angle and the number of H-bonds
between the pair; null for "cm" rules).

Other analyses add fields of their own (see `FrameAnalysis.report_frame` and
`report_cluster`), under their class name in a frame's or a cluster's record. The
`Nucleus` analysis adds to each cluster ``"Nucleus": {"nuclei": [<group>, ...],
"combined_dipole_moment": <dipole of all its nuclei together, or null>}``, a group
being a cluster without its ``id`` and ``connections``.

Floats are rounded to `DECIMALS` decimal places, what other analyses add included,
and NaN is written as null (a bare NaN isn't valid JSON). Units are those of `UNITS`.

Classes:
--------
- Report: A report file, read frame by frame, or as tables.

Functions:
----------
- read_report: Open a report file.
- report_name: The report's file name for a compression.
- dumps: Encode a record as one line of the report.
- rounded: Round every float in a value, as the report's are.
- group_record: Encode a group of molecules and its properties.
- cluster_record: Encode a cluster, its properties and connections.
"""

import gzip
import io
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from pathlib import Path
from typing import Any, BinaryIO

import numpy as np
import orjson
import pandas as pd
import zstandard
from loguru import logger

from .cluster import Cluster, MolGroup
from .config import ReportCompression

FORMAT = "molclusters-report"
VERSION = 2  # 1 was molclusters.json, a single JSON document

# decimal places of every float in the report: well below what the properties
# are good for (e.g. 1e-6 angstrom), yet well under the 17 digits of a float
DECIMALS = 6

UNITS = {
    "time": "ps",
    "mass": "amu",
    "volume": "angstrom^3",
    "radius": "angstrom",
    "diameter": "angstrom",
    "density": "g/cm^3",
    "charge": "e",
    "dipole_moment": "D",
    "distance": "angstrom",
    "angle": "degrees",
}

# the columns of a connection row; the ones after the resids are edge attributes
CONNECTION_COLUMNS = ("i", "j", "distance", "angle", "n_hbonds")

_SUFFIXES: dict[str, str] = {"zstd": ".zst", "gzip": ".gz", "none": ""}

_OPTIONS = orjson.OPT_APPEND_NEWLINE | orjson.OPT_SERIALIZE_NUMPY

_GZIP_MAGIC = b"\x1f\x8b"
_ZSTD_MAGIC = b"\x28\xb5\x2f\xfd"


def report_name(compression: ReportCompression) -> str:
    """Give the report's file name for a compression.

    Parameters
    ----------
    compression : ReportCompression
        "zstd", "gzip" or "none".

    Returns
    -------
    str
        ``molclusters.jsonl`` with the compression's suffix, if any.
    """
    return f"molclusters.jsonl{_SUFFIXES[compression]}"


def dumps(record: dict[str, Any]) -> bytes:
    """Encode a record as one line of the report.

    Parameters
    ----------
    record : dict[str, Any]
        The record: JSON types, and numpy scalars or arrays.

    Returns
    -------
    bytes
        The record's JSON, newline included.
    """
    return orjson.dumps(record, option=_OPTIONS)


def _round(value: float) -> float:
    """Round a float to `DECIMALS` decimal places.

    Returns
    -------
    float
        The rounded value, a Python float (NaN stays NaN, written as null).
    """
    return round(float(value), DECIMALS)


def _key(key: object) -> str:
    """Give a mapping's key as a JSON object's: a string.

    Returns
    -------
    str
        The key, an integer one (e.g. a resid) written in decimal, as the json
        module does.

    Raises
    ------
    TypeError
        If the key is neither a string nor an integer.
    """
    match key:
        case str():
            return key
        case int() | np.integer() if not isinstance(key, bool):
            return str(int(key))
    raise TypeError(f"A {type(key).__name__} key can't be written to the report.")


def rounded(value: Any) -> Any:  # noqa: ANN401 (any JSON-able value)
    """Round every float in a value to `DECIMALS` decimal places, however nested.

    For what analyses add to the report, whose floats are rounded as the report's
    own are.

    Parameters
    ----------
    value : Any
        JSON types (mappings, sequences, numbers, strings, None), numpy scalars or
        arrays.

    Returns
    -------
    Any
        The value, mappings as dicts, other sequences as lists, float arrays
        rounded.

    Raises
    ------
    TypeError
        If the value holds something that can't be written as JSON.
    """
    match value:
        case float() | np.floating():
            return _round(value)
        case None | bool() | int() | str() | np.integer() | np.bool_():
            return value
        case np.ndarray() if value.dtype.kind == "f":
            return np.round(value, DECIMALS)
        case np.ndarray() if value.dtype.kind in "biu":
            return np.ascontiguousarray(value)  # orjson takes contiguous arrays
        case Mapping():
            return {_key(key): rounded(item) for key, item in value.items()}
        case list() | tuple() | np.ndarray():
            return [rounded(item) for item in value]
    raise TypeError(f"A {type(value).__name__} can't be written to the report.")


def group_record(group: MolGroup) -> dict[str, Any]:
    """Encode a group of molecules and its properties.

    Parameters
    ----------
    group : MolGroup
        The group, e.g. a nucleus.

    Returns
    -------
    dict[str, Any]
        Its size, resids, composition and properties (see the module's docstring).
    """
    return {
        "size": group.size,
        "resids": sorted(map(int, group.resids)),
        "composition": {str(name): n for name, n in group.composition.items()},
        "mass": _round(group.mass),
        "volume": _round(group.volume),
        "radius": _round(group.radius),
        "diameter": _round(group.diameter),
        "density": _round(group.density),
        "charge": _round(group.charge),
        "dipole_moment": _round(group.dipole_moment),
        "sphericity": _round(group.sphericity),
        "shape_parameter": _round(group.shape_parameter),
    }


def _attribute(value: object) -> int | float | None:
    """Encode an edge attribute: counts stay whole numbers, the rest are rounded.

    Returns
    -------
    int | float | None
        The value, or None if the edge doesn't have the attribute.
    """
    if value is None or isinstance(value, int):
        return value
    if hasattr(value, "dtype") and value.dtype.kind in "iu":
        return int(value)
    return _round(value)


def connection_rows(cluster: Cluster) -> list[list[int | float | None]]:
    """Encode the connections of a cluster, as rows of `CONNECTION_COLUMNS`.

    Parameters
    ----------
    cluster : Cluster
        The cluster.

    Returns
    -------
    list[list[int | float | None]]
        One row per connection, in the order of the cluster's graph.
    """
    attributes = CONNECTION_COLUMNS[2:]
    return [
        [int(i), int(j), *(_attribute(edge.get(name)) for name in attributes)]
        for i, j, edge in cluster.graph.edges.data()
    ]


def cluster_record(cluster: Cluster) -> dict[str, Any]:
    """Encode a cluster, its properties and connections.

    Parameters
    ----------
    cluster : Cluster
        The cluster.

    Returns
    -------
    dict[str, Any]
        Its id, then `group_record`, then its connections.
    """
    return {
        "id": cluster.id,
        **group_record(cluster),
        "connections": connection_rows(cluster),
    }


def read_report(path: Path | str) -> "Report":
    """Open a report file.

    Parameters
    ----------
    path : Path | str
        The report, compressed or not.

    Returns
    -------
    Report
        The report, its header read.
    """
    return Report(path)


class Report:
    """A report file, read frame by frame, or as tables.

    The file is read again by each call, frame by frame, so even a report too
    large for memory can be gone through. Its compression is told from its
    first bytes, not its name.

    Attributes
    ----------
    path : Path
        The report file.
    header : dict[str, Any]
        The header: the software, inputs, units and configuration of the run.
    """

    __slots__ = ["path", "header"]

    def __init__(self, path: Path | str) -> None:
        """Open a report, reading its header.

        Parameters
        ----------
        path : Path | str
            The report, compressed or not.

        Raises
        ------
        ValueError
            If the file isn't a MolClusters report this version can read.
        """
        self.path = Path(path)
        with self.__lines() as lines:
            first = next(lines, b"")

        try:
            header = orjson.loads(first)
        except orjson.JSONDecodeError:
            header = None
        if not isinstance(header, dict) or header.get("format") != FORMAT:
            raise ValueError(
                f"{self.path} is not a MolClusters report. The molclusters.json of "
                "MolClusters 0.6 and earlier is a single JSON document: read it "
                "with the json module."
            )
        if header.get("version", 0) > VERSION:
            raise ValueError(
                f"{self.path} is a version {header['version']} report, and this "
                f"MolClusters reads up to version {VERSION}: update MolClusters."
            )
        self.header: dict[str, Any] = header

    @contextmanager
    def __lines(self) -> Iterator[Iterator[bytes]]:
        """Open the file, decompressed, as lines.

        Yields
        ------
        Iterator[bytes]
            The file's lines, newlines included.
        """
        with self.path.open("rb") as raw:
            magic = raw.read(4)
            raw.seek(0)
            stream: BinaryIO
            if magic.startswith(_GZIP_MAGIC):
                stream = gzip.GzipFile(fileobj=raw)
            elif magic == _ZSTD_MAGIC:
                reader = zstandard.ZstdDecompressor().stream_reader(
                    raw, read_across_frames=True, closefd=False
                )
                stream = io.BufferedReader(reader)
            else:
                stream = raw
            with stream:
                yield iter(stream)

    def frames(self) -> Iterator[dict[str, Any]]:
        """Read the frames, one at a time.

        A report that ends early (the run was interrupted, or killed) gives the
        frames it has, with a warning.

        Yields
        ------
        dict[str, Any]
            Each frame's record (see the module's docstring).
        """
        n_read = 0
        complete = True
        with self.__lines() as lines:
            next(lines)  # the header
            try:
                for line in lines:
                    if not line.endswith(b"\n"):
                        complete = False  # cut off while being written
                        break
                    yield orjson.loads(line)
                    n_read += 1
            except EOFError:
                complete = False  # a gzip stream that was never finished

        expected = self.header.get("n_frames")
        if not complete or (expected is not None and n_read < expected):
            of = "" if expected is None else f" of the {expected} of the run"
            logger.warning(
                f"{self.path} ends early, after {n_read} frame(s){of}: the run "
                "was interrupted before writing the rest."
            )

    def __iter__(self) -> Iterator[dict[str, Any]]:
        """Read the frames, one at a time (see `frames`).

        Returns
        -------
        Iterator[dict[str, Any]]
            Each frame's record.
        """
        return self.frames()

    def clusters(self) -> pd.DataFrame:
        """Tabulate the clusters: one row per cluster per frame.

        Returns
        -------
        pd.DataFrame
            ``frame`` and ``time``, then the cluster's fields but its connections
            (see `connections`); the composition gets a ``composition.<resname>``
            column per residue name (NaN when absent), and the fields another
            analysis added a ``<Analysis>.<field>`` column each.
        """
        rows = [
            {"frame": frame["frame"], "time": frame["time"], **cluster}
            for frame in self.frames()
            for cluster in frame["clusters"]
        ]
        for row in rows:
            del row["connections"]
        return pd.json_normalize(rows, max_level=1)

    def connections(self) -> pd.DataFrame:
        """Tabulate the connections: one row per connection per frame.

        Returns
        -------
        pd.DataFrame
            ``frame``, ``time`` and ``cluster`` (its id), then the header's
            ``connection_columns``; ``n_hbonds`` is a nullable integer (``<NA>``
            for "cm" rules).
        """
        columns = ["frame", "time", "cluster", *self.header["connection_columns"]]
        rows = [
            [frame["frame"], frame["time"], cluster["id"], *row]
            for frame in self.frames()
            for cluster in frame["clusters"]
            for row in cluster["connections"]
        ]
        table = pd.DataFrame(rows, columns=columns)
        if "n_hbonds" in table:
            table["n_hbonds"] = table["n_hbonds"].astype("Int64")
        return table
