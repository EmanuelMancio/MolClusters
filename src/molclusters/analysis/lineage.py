# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides `Lineage`, where the clusters came from, what became of them and when.

It writes cluster_events.csv, the clusters' formations, splits, merges and
dissolutions frame by frame, and cluster_lifetimes.csv, each cluster's birth,
end and lifetime.
"""

import math
from collections import Counter
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from loguru import logger

from ..output import OutputFile
from .base import Frame, FrameAnalysis, Run

EVENTS = "cluster_events.csv"
LIFETIMES = "cluster_lifetimes.csv"

EVENT_COLUMNS = ["Frame", "Time", "Event", "Cluster", "Other", "NMols"]

LIFETIME_COLUMNS = [
    "Id",
    "BirthTime",
    "DeathTime",
    "Lifetime",
    "NFrames",
    "BornAtStart",
    "AliveAtEnd",
    "Origin",
    "Parents",
    "Fate",
    "MergedInto",
    "BirthSize",
    "MaxSize",
    "LastSize",
]

# what each event row says, for the log
_DESCRIPTIONS = {
    "formation": "cluster {cluster} formed of {n} free molecule(s)",
    "split": "cluster {cluster} split off cluster {other} with {n} molecule(s)",
    "merge": "cluster {other} merged into cluster {cluster} with {n} molecule(s)",
    "dissolution": "cluster {cluster} ({n} molecule(s)) dissolved",
}


@dataclass(slots=True)
class _Life:
    """What `Lineage` knows of a cluster: its birth, its size so far and its end."""

    birth_time: float
    born_at_start: bool
    origin: str
    parents: tuple[int, ...]
    birth_size: int
    max_size: int
    last_size: int
    n_frames: int = 1
    death_time: float = math.nan
    fate: str = "alive"
    merged_into: int | None = None


class Lineage(FrameAnalysis):
    """Records where each cluster came from, what became of it, and when.

    Events are read from each frame's `Frame.transition` and stamped with the
    first frame that shows them, so a cluster's `DeathTime` is the time of the
    first frame it's gone from, and its lifetime, ``DeathTime - BirthTime``, a
    whole number of frame spacings (a cluster seen in a single frame lived one
    spacing).

    cluster_events.csv has one row per event (see `EVENT_COLUMNS`); ``Other`` is
    the other cluster involved, if any, and ``NMols`` the molecules it concerns:

    - ``formation``: `Cluster` is new, made only of free molecules (``NMols``, its
      size);
    - ``split``: `Cluster` is new, and took ``NMols`` molecules from cluster
      ``Other``. A cluster made of pieces of several clusters has a row for each,
      and one breaking into several pieces a row for each new piece;
    - ``merge``: cluster ``Other`` ended, merged into `Cluster` with ``NMols`` of
      its molecules; a merge of several clusters has a row for each one absorbed;
    - ``dissolution``: `Cluster` ended without merging (``NMols``, its last size).

    The clusters of the run's first frame have no event: whatever formed them
    happened before the run.

    cluster_lifetimes.csv has one row per cluster (see `LIFETIME_COLUMNS`):
    ``Origin`` is "initial" (there at the first frame), "formation" or "split",
    ``Parents`` the clusters it took molecules from (";"-separated, most
    molecules first), ``Fate`` "merge", "dissolution" or "alive" (at the last
    frame), and ``MergedInto`` the cluster that absorbed it. The lifetimes of
    clusters born at the start or alive at the end are censored: they lived at
    least that long. One alive at the end has no ``DeathTime``, and its
    ``Lifetime`` is up to the last frame.

    Each frame's `Frame.transition` also goes in the report (see
    `report_frame`), for the flows of molecules between clusters.

    Attributes
    ----------
    lives : dict[int, _Life]
        Cluster id -> what is known of the cluster so far, for every cluster seen.
    """

    __slots__ = ["lives", "_time"]

    outputs = (OutputFile(EVENTS), OutputFile(LIFETIMES))

    def __init__(self) -> None:
        """Initialize an empty record, started over by `prepare`."""
        self.lives: dict[int, _Life] = {}
        self._time = math.nan

    def prepare(self, run: Run) -> None:
        """Start the record over, and cluster_events.csv with its header.

        Parameters
        ----------
        run : Run
            The run about to start.
        """
        self.lives = {}
        self._time = math.nan
        run.output.append(EVENTS, ",".join(EVENT_COLUMNS) + "\n")

    def analyse(self, frame: Frame) -> None:
        """Record the current frame's events, and the clusters' sizes.

        Parameters
        ----------
        frame : Frame
            The current frame.
        """
        transition = frame.transition
        time = float(frame.time)
        self._time = time
        rows: list[tuple[str, int, int | None, int]] = []

        for cid in sorted(transition.born):
            size = frame.clusters[cid].size
            parents = sorted(
                ((prev, n) for prev, n in transition.sources(cid).items() if prev),
                key=lambda source: (-source[1], source[0]),
            )
            if frame.index == 0:
                origin = "initial"
            elif parents:
                origin = "split"
                rows.extend(("split", cid, prev, n) for prev, n in parents)
            else:
                origin = "formation"
                rows.append(("formation", cid, None, size))
            self.lives[cid] = _Life(
                birth_time=time,
                born_at_start=frame.index == 0,
                origin=origin,
                parents=tuple(prev for prev, _ in parents),
                birth_size=size,
                max_size=size,
                last_size=size,
            )

        for lost, into in sorted(transition.merged.items()):
            rows.append(("merge", into, lost, transition.flows[lost, into]))
            self._end(lost, time, "merge", into)
        for cid in sorted(transition.dissolved):
            rows.append(("dissolution", cid, None, self.lives[cid].last_size))
            self._end(cid, time, "dissolution")

        for cid, cls in frame.clusters.items():
            if cid not in transition.born:
                life = self.lives[cid]
                life.n_frames += 1
                life.last_size = cls.size
                life.max_size = max(life.max_size, life.last_size)

        number = int(frame.universe.coord.frame)
        for event, cid, other, n in rows:
            logger.debug(
                "Frame {}: " + _DESCRIPTIONS[event],
                number,
                cluster=cid,
                other=other,
                n=n,
            )
            other_text = "" if other is None else other
            frame.output.append(
                EVENTS, f"{number},{time!r},{event},{cid},{other_text},{n}\n"
            )

    def report_frame(self, frame: Frame) -> dict[str, Any]:
        """Give how the frame's clusters came from the previous frame's, for the report.

        Parameters
        ----------
        frame : Frame
            The current frame.

        Returns
        -------
        dict[str, Any]
            The frame's `Frame.transition`, sorted: ``flows``, rows of (previous
            id, current id, molecules), id 0 for no cluster; ``born``; ``merged``,
            rows of (absorbed id, id it merged into); and ``dissolved``.
        """
        transition = frame.transition
        return {
            "flows": [
                [prev, cur, n] for (prev, cur), n in sorted(transition.flows.items())
            ],
            "born": sorted(transition.born),
            "merged": [list(pair) for pair in sorted(transition.merged.items())],
            "dissolved": sorted(transition.dissolved),
        }

    def _end(
        self, cid: int, time: float, fate: str, merged_into: int | None = None
    ) -> None:
        """Record the end of a cluster.

        Parameters
        ----------
        cid : int
            The cluster's id.
        time : float
            The time of the first frame without it.
        fate : str
            "merge" or "dissolution".
        merged_into : int | None
            The cluster that absorbed it, for a merge.
        """
        life = self.lives[cid]
        life.death_time = time
        life.fate = fate
        life.merged_into = merged_into

    def lifetimes(self) -> pd.DataFrame:
        """Tabulate every cluster seen so far, as cluster_lifetimes.csv has them.

        Returns
        -------
        pd.DataFrame
            One row per cluster, by id (see `LIFETIME_COLUMNS`).
        """
        rows = [
            [
                cid,
                life.birth_time,
                life.death_time,
                (self._time if life.fate == "alive" else life.death_time)
                - life.birth_time,
                life.n_frames,
                life.born_at_start,
                life.fate == "alive",
                life.origin,
                ";".join(map(str, life.parents)),
                life.fate,
                life.merged_into,
                life.birth_size,
                life.max_size,
                life.last_size,
            ]
            for cid, life in sorted(self.lives.items())
        ]
        table = pd.DataFrame(rows, columns=LIFETIME_COLUMNS)
        table["MergedInto"] = table["MergedInto"].astype("Int64")
        return table

    def finish(self, run: Run) -> None:
        """Log a summary of the events and lifetimes, and write cluster_lifetimes.csv.

        Parameters
        ----------
        run : Run
            The run that just ended.
        """
        table = self.lifetimes()
        table.to_csv(run.output.path(LIFETIMES), index=False)

        origins = Counter(table["Origin"])
        fates = Counter(table["Fate"])
        message = (
            f"Cluster lineage: {origins['formation']} cluster(s) formed of free "
            f"molecules and {origins['split']} split off others; "
            f"{fates['merge']} merged into another and {fates['dissolution']} "
            "dissolved."
        )
        # lifetimes neither censored at the start nor at the end
        whole = table[~table["BornAtStart"] & ~table["AliveAtEnd"]]
        if len(whole):
            single = int((whole["NFrames"] == 1).sum())
            message += (
                f" The {len(whole)} born and ended within the run lived "
                f"{np.median(whole['Lifetime']):.4g} ps (median); {single} of them "
                "for a single frame."
            )
        logger.info(message)
