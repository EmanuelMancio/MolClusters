# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides the `MolClusters` class for analyzing molecular clusters in molecular dynamics simulations.

Classes:
--------
- MolClusters: A class for managing and analyzing molecular clusters, including solute-solvent interactions,
  cluster evolution, and nucleus analysis.
- MolClustersData: A helper class for encoding and storing cluster data for output.

Dependencies:
-------------
- MDAnalysis: For molecular dynamics trajectory and structure analysis.
- NetworkX: For graph-based operations on molecular clusters.
- NumPy: For numerical computations.
- Pandas: For data manipulation and exporting results.
- tqdm: For progress tracking during analysis.
"""

import json
import pathlib as path
import tempfile
from collections import Counter, defaultdict
from functools import reduce
from typing import Any, Generator

import MDAnalysis as mda
import networkx as nx
import numpy as np
import pandas as pd
from loguru import logger
from MDAnalysis import core
from tqdm import tqdm

from .cluster import Cluster, MDAResidueGroupAnalyzer
from .config import MolClsConfig
from .conntable import ConnectionTable
from .version import __version__

# TODO: create analysis class to declutter MolClusters


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
    sels : dict[str, core.groups.AtomGroup]
        Atom groups for each residue type.
    conntab : ConnectionTable
        The connectivity table for molecular clusters.
    clusters : dict[int, Cluster]
        A dictionary of detected clusters, keyed by cluster ID.
    mol_clt : dict[int, int]
        A mapping of molecule IDs to their respective cluster IDs.
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
        "sels",
        "conntab",
        "clusters",
        "mol_clt",
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
        self.sels: dict[str, core.groups.AtomGroup] = {
            res: self.uni.select_atoms(f"resname {res}")
            for res in config._rules.all_keys()
        }

        self.conntab = ConnectionTable(self.uni, self.config._rules, self.sels)
        self.clusters: dict[int, Cluster] = {}
        self.mol_clt: dict[int, int] = {}

        self.clusters_size_evo = np.zeros((len(self.uni.trajectory), 5))
        self.radius_evolution = {}

        # TODO: move start to run
        self.__start_clusters()
        self.__get_clusters_info(0)
        if self.config.nucleus is not None:
            self.nucleus_data = np.empty(
                (len(self.uni.trajectory), 9)
            )  # Value 9 accounts for time column and 8 property columns
            self.nucleus_data.fill(np.nan)
            self.__nucleus_analysis(0)

        self.data_holder = MolClustersData(self)
        self.data_holder.parse_frame()

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
        for cid, cls in self.clusters.items():
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
            id for sel in self.config.solute for id in self.sels[sel].residues.resids
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
        for cls in self.clusters.values():
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

    def __start_clusters(self) -> None:
        """Initialize clusters at the beginning of the analysis."""
        for subconn in self.conntab.subconntables():
            if self.config.is_ignored_composition(subconn.resnames):
                continue

            cls_id = self.__create_new_cluster(subconn)

            for mol in subconn:
                self.mol_clt[mol] = cls_id

    # TODO: make a better name for this function
    def __gen_origin_cluster_counter(
        self, subconn: ConnectionTable._SubConnTable
    ) -> Counter:
        """Count where the molecules of a new connected group came from.

        Each molecule of `subconn` (a connected group of the current frame) is
        looked up in `mol_clt`, which still holds the *previous* frame's assignment,
        so the result says how many of the connected group's molecules each
        previous cluster contributed. Molecules that were free in the previous frame
        are counted under id ``0`` (safe as a sentinel because cluster ids start
        at 1).

        For example, a connected group made of 4 molecules of cluster 7, 2 of
        cluster 9 and 1 free molecule gives ``Counter({7: 4, 9: 2, 0: 1})``.

        Parameters
        ----------
        subconn : ConnectionTable._SubConnTable
            A connected group of the current frame.

        Returns
        -------
        Counter
            Previous-frame cluster id (``0`` for free molecules) -> number of the
            connected group's molecules that came from it.
        """
        mols_origin_clusters = {mol: self.mol_clt.get(mol, 0) for mol in subconn}
        return Counter(mols_origin_clusters.values())

    @staticmethod
    def __best_groups(
        conn_info: list[tuple[ConnectionTable._SubConnTable, Counter]],
    ) -> dict[int, int]:
        """Find the connected group that best continues each previous cluster.

        A previous cluster's best group is the one holding the most of its
        molecules. On a tie, the smallest group wins, i.e. the purest one (the
        largest share of its molecules came from the cluster), so a fragment made
        only of the cluster's own molecules beats one mixed with newcomers. If the
        groups are also the same size, the first one in `conn_info` wins.

        A contribution of a single molecule is ignored: one molecule joining other
        molecules doesn't carry the cluster's identity. A cluster with no
        contribution of two or more molecules has no best group and dies.

        Every group is looked at before any choice is made, so the result doesn't
        depend on the order the groups are processed in.

        Parameters
        ----------
        conn_info : list[tuple[ConnectionTable._SubConnTable, Counter]]
            Every connected group of the current frame, largest first, each with its
            origin counter (see `__gen_origin_cluster_counter`).

        Returns
        -------
        dict[int, int]
            Previous-frame cluster id -> index in `conn_info` of its best group.
        """
        best: dict[int, tuple[tuple[int, int], int]] = {}

        for i, (subconn, origin_clusters) in enumerate(conn_info):
            for cls_id, n in origin_clusters.items():
                if not cls_id or n == 1:  # free molecules, or a single molecule
                    continue

                key = (n, -len(subconn))
                if cls_id not in best or key > best[cls_id][0]:
                    best[cls_id] = (key, i)

        return {cls_id: i for cls_id, (_, i) in best.items()}

    def __create_new_cluster(self, subconn: ConnectionTable._SubConnTable) -> int:
        """Create a new cluster from a subconnection table.

        Parameters
        ----------
        subconn : ConnectionTable._SubConnTable
            The subconnection table from which to create the new cluster.

        Returns
        -------
        int
            The ID of the newly created cluster.
        """
        cluster = Cluster(self.uni, subconn)
        id = cluster.id
        self.clusters[id] = cluster

        return id

    @staticmethod
    def __get_older_cluster(candidates: list[int], origin_clusters: Counter) -> int:
        """Pick which candidate cluster keeps its id in a merge.

        The candidate that contributed the most molecules wins, so a large cluster
        absorbing a small one keeps its id even if it is younger. If several
        candidates tie for the most molecules, the oldest wins. Cluster ids are
        handed out in increasing order, so the oldest is the lowest id. The order of
        `candidates` doesn't matter.

        Parameters
        ----------
        candidates : list[int]
            Non-empty list of the previous clusters whose best group (see
            `__best_groups`) is this connected group.
        origin_clusters : Counter
            Origin counter of the connected group.

        Returns
        -------
        int
            Id of the cluster that continues as the connected group.
        """
        return max(candidates, key=lambda c: (origin_clusters[c], -c))

    def __update_clusters(self) -> None:
        """Reconcile the current frame's connected groups with the previous clusters.

        This is the dominance algorithm, which keeps cluster ids stable while
        clusters grow, shrink, split and merge.

        A *connected group* (``subconn`` in the code, a
        `ConnectionTable._SubConnTable`) is a set of molecules linked to each
        other, directly or through other molecules, by the config's rules in the
        current frame, and to nothing outside the set. Free (unlinked) molecules
        form no connected group. A connected group has no id yet; a *cluster* is
        what it becomes once this method gives it one, new or inherited from the
        previous frame.

        The connection table is rebuilt for the current frame, and the previous
        clusters that the molecules of each connected group (ignored compositions
        excluded) came from are counted (see `__gen_origin_cluster_counter`). Ids
        are then assigned in two steps, each of which sees the whole frame:

        1. Each previous cluster picks its best group: the one holding the most of
           its molecules, the purest on a tie (see `__best_groups`).
        2. Each connected group chosen by one or more previous clusters continues
           the one that contributed the most molecules, the oldest on a tie (see
           `__get_older_cluster`). A group chosen by nobody is a new cluster.

        The usual events follow from these two rules:

        - formation: a group of free molecules is chosen by nobody -> new;
        - growth or shrinking: a group chosen by one cluster continues it;
        - split: the cluster continues in its best fragment and the other fragments,
          chosen by nobody, are new;
        - merge: several clusters choose the same group, one continues and the
          others die, even if they also left a remnant elsewhere;
        - dissolution: a cluster that contributed at most one molecule to every
          group has no best group and dies.

        So a mixed dimer (one molecule from each origin) is always new, and each id
        continues in at most one group. Previous clusters that no connected group
        continues are dropped, along with the `mol_clt` entries of molecules that
        are now free.
        """
        self.conntab.update()

        modified_mols = set()
        modified_clusters = set()

        conn_info = [
            (sub, self.__gen_origin_cluster_counter(sub))
            for sub in self.conntab.subconntables()
            if not self.config.is_ignored_composition(sub.resnames)
        ]

        candidates = defaultdict(list)
        for cls_id, i in self.__best_groups(conn_info).items():
            candidates[i].append(cls_id)

        for i, (subconn, origin_clusters) in enumerate(conn_info):
            if candidates[i]:
                id = self.__get_older_cluster(candidates[i], origin_clusters)
                self.clusters[id].update_from_conntable(subconn)
            else:
                id = self.__create_new_cluster(subconn)

            for mol in subconn:
                self.mol_clt[mol] = id

            modified_clusters.add(id)
            modified_mols.update(subconn)

        for mol in set(self.mol_clt.keys()).difference(modified_mols):
            self.mol_clt.pop(mol)

        # TODO: deal with clusters that weren't modified. Needs to consider that some clusters merged (for log filing)  # noqa: E501
        for cls in set(self.clusters.keys()).difference(modified_clusters):
            self.clusters.pop(cls)

    def __write_coordinates(self) -> None:
        """Write the coordinates of clusters to files."""
        for cls in self.clusters.values():
            if not set(self.config.solute).intersection(set(cls.resnames)):
                continue

            with tempfile.TemporaryDirectory() as tmp_dir:
                tmp_path = path.Path(tmp_dir) / "cluster.gro"

                with (
                    mda.Writer(str(tmp_path), multiframe=False) as w,
                    cls.whole() as atoms,
                ):
                    w.write(atoms.sort())

                with tmp_path.open() as tmp:
                    dt = tmp.readlines()

            dt[0] = f"Cluster-{cls.id} - Time = {self.uni.coord.time}\n"

            # pooled by size: an ensemble of what an N-mer looks like, across all
            # clusters that were ever that size, independent of cluster identity
            with path.Path(f"cls-n{cls.size}.gro").open("a+") as out:
                out.write("".join(dt))

            # pooled by identity: this specific cluster's own trajectory, tracked
            # across frames via the dominance algorithm regardless of size changes
            with path.Path(f"cls-id{cls.id}.gro").open("a+") as out:
                out.write("".join(dt))

            # TODO: change to support merges
            # FIXME: with changes in config this needs to be updated
            if self.config._follow_solute:
                sol_id: set[int] = set(cls.resids).intersection(set(self.solutes))
                if len(sol_id) == 0:
                    logger.error(f"Cluster {cls.id}: expected a solute, found none")
                    continue
                elif len(sol_id) > 1:
                    # TODO: make more feature-rich follow procedure
                    logger.warning(
                        f"Cluster {cls.id}: more than one solute, will not follow"
                    )
                    continue

                sol_id = sol_id.pop()

                with path.Path(f"solute-{sol_id}.gro").open("a+") as out:
                    out.write("".join(dt))

    # TODO: break into single_step function to better use in MDRHConstant
    def run(self) -> None:
        """Run the molecular cluster analysis.

        This method performs cluster detection, solute-solvent analysis, nucleus analysis,
        and exports the results to files.
        """
        if self.config.solute is not None:
            self.__start_solute_solvent()
            self.__solute_solvent_analysis(0)

        with tqdm(
            total=len(self.uni.trajectory[1:]), initial=1, mininterval=5, miniters=10
        ) as pbar:
            for i, _ in enumerate(self.uni.trajectory[1:], start=1):
                self.__update_clusters()
                self.__get_clusters_info(i)
                if self.config.solute is not None:
                    self.__solute_solvent_analysis(i)
                    self.__write_coordinates()

                if self.config.nucleus is not None:
                    self.__nucleus_analysis(i)

                self.data_holder.parse_frame()
                pbar.update()

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

        with path.Path("molclusters.json").open("w+") as json_out:
            json.dump(self.data_holder.data, json_out, indent=2)

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
        return self.mol_clt.get(mol, False)

    def __get_clusters_info(self, k: int) -> None:
        """Update cluster size evolution information for a given frame.

        Parameters
        ----------
        k : int
            The frame index.
        """
        sizes = [cls.size for cls in self.clusters.values()]
        if len(sizes) == 0:
            avg, min_size, max_size = 0, 0, 0
        else:
            avg = np.average(sizes)
            min_size = min(sizes)
            max_size = max(sizes)

        time = self.uni.coord.time

        self.clusters_size_evo[k][0] = time
        self.clusters_size_evo[k][1] = len(self.clusters)
        self.clusters_size_evo[k][2] = min_size
        self.clusters_size_evo[k][3] = avg
        self.clusters_size_evo[k][4] = max_size

    def __print_clusters_index(self) -> None:
        """Write ndx file with cluster index."""
        cols = 15
        i = 1
        with path.Path("clusters_index.ndx").open("w+", encoding="utf-8") as ndx:
            for id, cluster in self.clusters.items():
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
        data["NClusters"] = len(self.molcls.clusters)

        molclusters_data = []
        for cid, cls in self.molcls.clusters.items():
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
