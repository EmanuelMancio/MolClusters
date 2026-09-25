# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides the `MolClusters` class for analyzing molecular clusters in molecular dynamics simulations.

Classes:
--------
- MolClusters: A class for managing and analyzing molecular clusters, including solute-solvent interactions,
  cluster evolution, and nucleus analysis.
- MolClustersData: A helper class for encoding and storing cluster data for output.

The clusters themselves, and their ids, are tracked by `tracker.ClusterTracker`.

Dependencies:
-------------
- MDAnalysis: For molecular dynamics trajectory and structure analysis.
- NetworkX: For graph-based operations on molecular clusters.
- NumPy: For numerical computations.
- Pandas: For data manipulation and exporting results.
- tqdm: For progress tracking during analysis.
"""

import io
import json
import pathlib as path
import re
import time
from collections import Counter, defaultdict
from functools import reduce
from typing import Any, Generator

import MDAnalysis as mda
import networkx as nx
import numpy as np
import pandas as pd
from loguru import logger
from MDAnalysis.lib.util import NamedStream
from tqdm import tqdm

from .cluster import Cluster, MDAResidueGroupAnalyzer
from .config import MolClsConfig
from .log import FILE_ONLY, format_duration
from .tracker import ClusterTracker
from .version import __version__

# TODO: create analysis class to declutter MolClusters


class _AppendBuffer:
    """Collects text to append to files, writing each file in one go when flushed.

    Opening a file costs milliseconds on Windows (antivirus scanning, on any file
    size), which dominated runs that appended every frame to its .gro files. Used
    as a context manager, it flushes on exit.

    Attributes
    ----------
    max_chars : int
        The buffered size, in characters, that triggers a flush.
    """

    __slots__ = ["max_chars", "_pending", "_size"]

    def __init__(self, max_chars: int = 64 * 2**20) -> None:
        self.max_chars = max_chars
        self._pending: defaultdict[str, list[str]] = defaultdict(list)
        self._size = 0

    def append(self, filename: str, text: str) -> None:
        """Queue `text` to be appended to `filename`, flushing if the buffer is full."""
        self._pending[filename].append(text)
        self._size += len(text)
        if self._size >= self.max_chars:
            self.flush()

    def flush(self) -> None:
        """Append the queued text to its files, opening each file once."""
        for filename, chunks in self._pending.items():
            with path.Path(filename).open("a+") as out:
                out.writelines(chunks)
        self._pending.clear()
        self._size = 0

    def __enter__(self) -> "_AppendBuffer":
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.flush()


class MolClusters:
    """A class for analyzing molecular clusters in molecular dynamics simulations.

    This class manages molecular clusters, performs solute-solvent analysis, tracks cluster evolution,
    and performs nucleus analysis.

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
    clusters_size_evo : np.ndarray
        An array tracking the evolution of cluster sizes over time.
    solutes : list[int]
        A list of solute molecule IDs.
    solvents : list[str]
        A list of solvent residue names.
    solute_data : np.ndarray
        An array storing solute-solvent analysis results.
    """

    __slots__ = [
        "uni",
        "config",
        "tracker",
        "clusters_size_evo",
        "radius_evolution",
        "solutes",
        "solvents",
        "solute_resnames",
        "solvent_resnames",
        "solute_data",
        "nucleus_data",
        "nucleus_holder",
        "data_holder",
        "follow_skipped",
        "gro_out",
    ]

    def __init__(self, universe: mda.Universe, config: MolClsConfig) -> None:
        """Initialize the MolClusters object.

        Parameters
        ----------
        universe : mda.Universe
            The MDAnalysis Universe object associated with the simulation.
        config : dict
            The configuration dictionary containing analysis settings.
        """
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

        self.clusters_size_evo = np.zeros((len(self.uni.trajectory), 5))
        self.radius_evolution = {}
        self.gro_out = _AppendBuffer()

        self.__get_clusters_info(0)
        if self.config.nucleus is not None:
            self.nucleus_data = np.empty(
                (len(self.uni.trajectory), 9)
            )  # Value 9 accounts for time column and 8 property columns
            self.nucleus_data.fill(np.nan)
            self.__nucleus_analysis(0)

        self.data_holder = MolClustersData(self)
        self.data_holder.parse_frame()

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

        The per-frame .gro files are appended to, so an earlier run's frames would
        end up mixed with this one's; the other outputs are just overwritten.
        """
        cwd = path.Path.cwd()
        names = ["evo.txt", "molclusters.json"]
        patterns = []
        if self.config.solute is not None:
            names.append("solute_solvent.csv")
            patterns += [r"cls-n\d+\.gro", r"cls-id\d+\.gro"]
            if self.config._follow_solute:
                patterns.append(r"solute-\d+\.gro")
        if self.config.nucleus is not None:
            names.append("nucleus_data.csv")

        overwritten = [name for name in names if (cwd / name).exists()]
        if overwritten:
            logger.info(f"Overwriting results of an earlier run: {overwritten}")

        appended = sorted(
            file.name
            for file in cwd.glob("*.gro")
            if any(re.fullmatch(pattern, file.name) for pattern in patterns)
        )
        if appended:
            shown = ", ".join(appended[:5]) + (", ..." if len(appended) > 5 else "")
            logger.warning(
                f"{len(appended)} .gro file(s) from an earlier run are in {cwd} "
                f"({shown}): this run appends its frames to them, mixing both runs. "
                "Move or delete them first to keep the runs apart."
            )

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

    def __nucleus_analysis(self, frame: int) -> None:
        """Perform nucleus analysis for a given frame.

        Parameters
        ----------
        frame : int
            index of the current frame in the trajectory.
        """
        self.nucleus_holder: dict[int, list[MDAResidueGroupAnalyzer]] = {}
        n_nucleus = []
        sizes = []
        radius = []
        dipole = []
        density = []
        sphericity = []
        shape = []
        charge = []
        for cid, cls in self.tracker.clusters.items():
            possible_nucleus = []
            for rnm, rid in zip(cls.resnames, cls.resids, strict=True):
                if rnm in self.config.nucleus:
                    possible_nucleus.append(rid)

            self.nucleus_holder[cid] = []
            subcomps = nx.induced_subgraph(cls.cluster, possible_nucleus)
            n_nuc = 0
            for sg in nx.connected_components(subcomps):
                tp = MDAResidueGroupAnalyzer(self.uni, list(sg))

                n_nuc += 1
                sizes.append(tp.size)
                radius.append(tp.radius_of_gyration)
                dipole.append(tp.dipole_moment)
                density.append(tp.density)
                sphericity.append(tp.sphericity)
                shape.append(tp.shape_parameter)
                charge.append(tp.charge)

                self.nucleus_holder[cid].append(tp)

            if n_nuc != 0:
                n_nucleus.append(n_nuc)

        self.nucleus_data[frame][0] = self.uni.coord.time
        self.nucleus_data[frame][1] = (
            0 if len(n_nucleus) == 0 else np.average(n_nucleus)
        )
        self.nucleus_data[frame][2] = (
            np.nan if len(n_nucleus) == 0 else np.average(sizes)
        )
        self.nucleus_data[frame][3] = (
            np.nan if len(n_nucleus) == 0 else np.average(radius)
        )
        self.nucleus_data[frame][4] = (
            np.nan if len(n_nucleus) == 0 else np.average(density)
        )
        self.nucleus_data[frame][5] = (
            np.nan if len(n_nucleus) == 0 else np.average(charge)
        )
        self.nucleus_data[frame][6] = (
            np.nan if len(n_nucleus) == 0 else np.average(dipole)
        )
        self.nucleus_data[frame][7] = (
            np.nan if len(n_nucleus) == 0 else np.average(sphericity)
        )
        self.nucleus_data[frame][8] = (
            np.nan if len(n_nucleus) == 0 else np.average(shape)
        )

    def __start_solute_solvent(self) -> None:
        """Initialize solute-solvent analysis."""
        self.solutes: list[int] = [
            id
            for sel in self.config.solute
            for id in self.tracker.sels[sel].residues.resids
        ]
        self.solvents = self.config.solvent

        self.solute_resnames = set(self.config.solute)
        self.solvent_resnames = set(self.config.solvent)

        self.solute_data = np.zeros(
            (len(self.uni.trajectory), 10)
        )  # Value 8 accounts for time column and 7 property columns

    def __solute_solvent_clusters(self) -> Generator[Cluster, None, None]:
        """Generate clusters that contain both solute and solvent residues.

        Yields
        ------
        Generator[Cluster]
            A generator that yields clusters containing both solute and solvent residues.
        """
        for cls in self.tracker.clusters.values():
            res = set(cls.resnames)
            if (
                len(res.intersection(self.solute_resnames)) > 0
                and len(res.intersection(self.solvent_resnames)) > 0
            ):
                yield cls

    def __solute_solvent_analysis(self, frame: int) -> None:
        """Perform solute-solvent analysis for a given frame.

        Parameters
        ----------
        frame : int
            index of the current frame in the trajectory.
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

        for cls in self.__solute_solvent_clusters():
            n_cls += 1

            n_solv = 0
            n_solt = 0
            mol_pop = Counter(cls.resnames)
            for solvent in self.solvent_resnames:
                n_solv += mol_pop.get(solvent, 0)

            for solute in self.solute_resnames:
                n_solt += mol_pop.get(solute, 0)

            n_solvents.append(n_solv)
            n_solutes.append(n_solt)

            radius.append(cls.radius_of_gyration)
            dipole.append(cls.dipole_moment)
            density.append(cls.density)
            sphericity.append(cls.sphericity)
            shape.append(cls.shape_parameter)
            charge.append(cls.charge)

        self.solute_data[frame][0] = self.uni.coord.time
        self.solute_data[frame][1] = n_cls
        self.solute_data[frame][2] = 0 if n_cls == 0 else np.average(n_solutes)
        self.solute_data[frame][3] = 0 if n_cls == 0 else np.average(n_solvents)
        self.solute_data[frame][4] = np.nan if n_cls == 0 else np.average(radius)
        self.solute_data[frame][5] = np.nan if n_cls == 0 else np.average(density)
        self.solute_data[frame][6] = np.nan if n_cls == 0 else np.average(charge)
        self.solute_data[frame][7] = np.nan if n_cls == 0 else np.average(dipole)
        self.solute_data[frame][8] = np.nan if n_cls == 0 else np.average(sphericity)
        self.solute_data[frame][9] = np.nan if n_cls == 0 else np.average(shape)

    def __write_coordinates(self) -> None:
        """Write the coordinates of clusters to files."""
        solutes = set(self.solutes)
        for cls in self.tracker.clusters.values():
            sol_ids = solutes.intersection(cls.resids)
            if not sol_ids:
                continue

            # the GRO writer only writes whole files, so render the frame in memory
            # (NamedStream keeps the buffer open when the writer closes it) and
            # append it to each output below
            buf = io.StringIO()
            with (
                mda.Writer(NamedStream(buf, "cluster.gro"), multiframe=False) as w,
                cls.whole() as atoms,
            ):
                w.write(atoms.sort())

            # replace the writer's fixed "Written by MDAnalysis" title line
            _, body = buf.getvalue().split("\n", 1)
            frame = f"Cluster-{cls.id} - Time = {self.uni.coord.time}\n{body}"

            # pooled by size: an ensemble of what an N-mer looks like, across all
            # clusters that were ever that size, independent of cluster identity
            self.gro_out.append(f"cls-n{cls.size}.gro", frame)

            # pooled by identity: this specific cluster's own trajectory, tracked
            # across frames via the dominance algorithm regardless of size changes
            self.gro_out.append(f"cls-id{cls.id}.gro", frame)

            # TODO: change to support merges
            # FIXME: with changes in config this needs to be updated
            if self.config._follow_solute:
                if len(sol_ids) > 1:
                    # TODO: make more feature-rich follow procedure
                    logger.debug(
                        f"Cluster {cls.id}: more than one solute, will not follow"
                    )
                    self.follow_skipped[cls.id] += 1
                    continue

                (sol_id,) = sol_ids

                self.gro_out.append(f"solute-{sol_id}.gro", frame)

    # TODO: break into single_step function to better use in MDRHConstant
    def run(self) -> None:
        """Run the molecular cluster analysis.

        This method performs cluster detection, solute-solvent analysis, nucleus analysis,
        and exports the results to files.
        """
        n_frames = len(self.uni.trajectory)
        self.__check_previous_outputs()
        logger.info(f"Tracking clusters over {n_frames} frame(s)")
        self.follow_skipped: Counter[int] = Counter()
        start = time.perf_counter()
        # the bar only shows on the terminal, so the log files get a line every 10%
        progress_step = max(1, n_frames // 10)

        if self.config.solute is not None:
            self.__start_solute_solvent()
            self.__solute_solvent_analysis(0)

        # the .gro frames are buffered, so write the ones already rendered even if the
        # run is interrupted, as appending them frame by frame used to
        with (
            tqdm(total=n_frames, initial=1, mininterval=5, miniters=10) as pbar,
            self.gro_out,
        ):
            for i, _ in enumerate(self.uni.trajectory[1:], start=1):
                self.tracker.update()
                self.__get_clusters_info(i)
                if self.config.solute is not None:
                    self.__solute_solvent_analysis(i)
                    self.__write_coordinates()

                if self.config.nucleus is not None:
                    self.__nucleus_analysis(i)

                self.data_holder.parse_frame()
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

        if self.follow_skipped:
            logger.warning(
                f"Solutes were not followed in {self.follow_skipped.total()} "
                f"frame(s) of {len(self.follow_skipped)} cluster(s) holding more than "
                "one solute: those frames are missing from the solute-<resid>.gro "
                "files. The cluster ids are logged at DEBUG level (--log-level DEBUG)."
            )
            logger.debug(f"Frames not followed, by cluster id: {self.follow_skipped}")

        n_clusters = self.clusters_size_evo[:, 1]
        logger.info(
            f"Found {n_clusters.mean():.1f} cluster(s) per frame on average "
            f"({n_clusters.min():.0f}-{n_clusters.max():.0f}); the largest held "
            f"{self.clusters_size_evo[:, 4].max():.0f} molecule(s)."
        )

        outputs = ["evo.txt", "molclusters.json"]

        np.savetxt(
            "evo.txt",
            self.clusters_size_evo,
            header="Time NClusters MinSize AvgSize MaxSize",
        )
        if self.config.solute is not None:
            self.solute_data = pd.DataFrame(
                self.solute_data,
                columns=[
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
                ],
            )
            self.solute_data.to_csv("solute_solvent.csv", index=False)
            outputs += ["solute_solvent.csv", "cls-n<size>.gro", "cls-id<id>.gro"]
            if self.config._follow_solute:
                outputs.append("solute-<resid>.gro")

        if self.config.nucleus is not None:
            self.nucleus_data = pd.DataFrame(
                self.nucleus_data,
                columns=[
                    "Time",
                    "NNuc",
                    "Size",
                    "Radius",
                    "Density",
                    "Charge",
                    "Dipole",
                    "Spher",
                    "Shape",
                ],
            )
            self.nucleus_data.to_csv("nucleus_data.csv", index=False)
            outputs.append("nucleus_data.csv")

        with path.Path("molclusters.json").open("w+") as json_out:
            json.dump(self.data_holder.data, json_out, indent=2)

        logger.info(f"Results written to {path.Path.cwd()}: {', '.join(outputs)}")

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

    def __get_clusters_info(self, k: int) -> None:
        """Update cluster size evolution information for a given frame.

        Parameters
        ----------
        k : int
            The frame index.
        """
        sizes = [cls.size for cls in self.tracker.clusters.values()]
        if len(sizes) == 0:
            avg, min_size, max_size = 0, 0, 0
        else:
            avg = np.average(sizes)
            min_size = min(sizes)
            max_size = max(sizes)

        time = self.uni.coord.time

        self.clusters_size_evo[k][0] = time
        self.clusters_size_evo[k][1] = len(self.tracker.clusters)
        self.clusters_size_evo[k][2] = min_size
        self.clusters_size_evo[k][3] = avg
        self.clusters_size_evo[k][4] = max_size

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


# TODO: use orjson for better encoding options
class MolClustersData:
    """A helper class for encoding and storing molecular cluster data.

    Attributes
    ----------
    molcls : MolClusters
        The parent MolClusters object.
    data : dict
        A dictionary for storing encoded cluster data.
    """

    def __init__(self, molclusters: MolClusters) -> None:
        """Initialize the MolClustersData object.

        Parameters
        ----------
        molclusters : MolClusters
            The parent MolClusters object.
        """
        self.molcls = molclusters

        conf = self.molcls.config.model_dump()

        self.data = {
            "Software": f"MolClusters {__version__}",
            "Trajectory": str(
                path.Path(self.molcls.uni.trajectory.filename).absolute()
            ),
            "Topology": str(path.Path(self.molcls.uni.filename).absolute()),
            "Config": conf,
            "MolClusters": [],
        }

    def parse_frame(self) -> None:
        """Parse the current frame and encode cluster data."""
        data = {}

        data["Time"] = self.molcls.uni.coord.time
        data["Frame"] = self.molcls.uni.coord.frame
        data["NClusters"] = len(self.molcls.tracker.clusters)

        molclusters_data = []
        for cid, cls in self.molcls.tracker.clusters.items():
            cls_data = self.encode_cluster(cls)
            if self.molcls.config.nucleus is not None:
                cls_data["Nucleus"] = []
                if self.molcls.nucleus_holder.get(cid, False):
                    nuclei = reduce(lambda a, b: a + b, self.molcls.nucleus_holder[cid])
                    cls_data["NucleiDipole"] = nuclei.dipole_moment
                    for nuc in self.molcls.nucleus_holder[cid]:
                        cls_data["Nucleus"].append(MolClustersData.encode_nucleus(nuc))

            molclusters_data.append(cls_data)

        data["Clusters"] = molclusters_data
        self.data["MolClusters"].append(data)

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
        MolClustersData.encode_properties(cls, data)
        return data

    @staticmethod
    def encode_nucleus(nuc: MDAResidueGroupAnalyzer) -> dict:
        """Encode a nucleus into a dictionary.

        Parameters
        ----------
        nuc : MDAResidueGroupAnalyzer
            The nucleus to encode.

        Returns
        -------
        dict
            The encoded nucleus data.
        """
        data = {}
        MolClustersData.encode_properties(nuc, data)
        return data

    @staticmethod
    def encode_properties(obj: Cluster | MDAResidueGroupAnalyzer, data: dict) -> None:
        """Encode the properties of a cluster or nucleus.

        Parameters
        ----------
        obj : Cluster | MDAResidueGroupAnalyzer
            The object to encode.
        data : dict
            The dictionary to store the encoded properties.
        """
        data["Size"] = obj.size
        data["Composition"] = MolClustersData.encode_composition(obj)

        if isinstance(obj, Cluster):
            data["Connections"] = MolClustersData.encode_connections(obj)

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
    def encode_composition(obj: Cluster | MDAResidueGroupAnalyzer) -> list[dict]:
        """Encode the composition of a cluster or nucleus.

        Parameters
        ----------
        obj : Cluster | MDAResidueGroupAnalyzer
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
            (int(edge[0]), int(edge[1]), {k: float(v) for k, v in edge[2].items()})
            for edge in obj.cluster.edges.data()
        ]
