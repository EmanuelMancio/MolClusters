#!python

"""Cluster Analyzer Script."""

from collections import Counter
from typing import Dict, Type, Union

import json
import MDAnalysis as mda
import networkx as nx
import numpy as np
import pandas as pd
import pathlib as path

from MDAnalysis import core
from tqdm import tqdm

from .cluster import Cluster, MDAAtomGroupAnalyzer
from .conntable import ConnTable
from . import __version__

# TODO: Create analysis class to declutter MolClusters


class MolClusters:
    __slots__ = [
        "uni",
        "config",
        "sels",
        "conntab",
        "clusters",
        "mol_clt",
        "clusters_size_evo",
        "solutes",
        "solvents",
        "solute_data",
        "__dict__",
    ]

    def __init__(
        self,
        universe: Type[mda.Universe],
        config,
    ) -> None:
        self.uni = universe
        self.config = config
        self.sels: Dict[str, Type[core.groups.AtomGroup]] = {
            res: self.uni.select_atoms(f"resname {res}") for res in config["rules"]
        }

        self.conntab = ConnTable(self.uni, self.config["rules"], self.sels)
        self.clusters: Dict[int, Type[Cluster]] = {}
        self.mol_clt: Dict[int, int] = {}

        self.clusters_size_evo = np.zeros((len(self.uni.trajectory), 5))
        self.radius_evolution = {}

        # TODO: move start to run function
        self.__start_clusters()
        self.__get_clusters_info(0)
        if self.config.get("solute", False):
            self.nucleus_data = np.empty(
                (len(self.uni.trajectory), 9)
            )  # Value 8 accounts for time column and 7 property columns
            self.nucleus_data.fill(np.nan)
            self.__nucleus_analysis(0)

        self.data_holder = MolClustersData(self)
        self.data_holder.parse_frame()

    def __nucleus_analysis(self, frame):
        self.nucleus_holder = {}
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
            for rnm, rid in zip(cls.resnames, cls.resids):
                if rnm in self.config["nucleus"]:
                    possible_nucleus.append(rid)

            self.nucleus_holder[cid] = []
            subcomps = nx.induced_subgraph(cls.cluster, possible_nucleus)
            n_nuc = 0
            for sg in nx.connected_components(subcomps):
                tp = MDAAtomGroupAnalyzer(self.uni, list(sg))

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
        self.nucleus_data[frame][1] = np.average(n_nucleus)
        self.nucleus_data[frame][2] = np.average(sizes)
        self.nucleus_data[frame][3] = np.average(radius)
        self.nucleus_data[frame][4] = np.average(density)
        self.nucleus_data[frame][5] = np.average(charge)
        self.nucleus_data[frame][6] = np.average(dipole)
        self.nucleus_data[frame][7] = np.average(sphericity)
        self.nucleus_data[frame][8] = np.average(shape)

    def __start_solute_solvent(self):
        self.solutes = [
            id for sel in self.config["solute"] for id in self.sels[sel].residues.resids
        ]
        self.solvents = self.config[
            "solvent"
        ]  # TODO: solvents should be automatically discovered but file takes precedent

        self.solute_resnames = set(self.config["solute"])
        self.solvent_resnames = set(self.config["solvent"])

        self.solute_data = np.zeros(
            (len(self.uni.trajectory), 10)
        )  # Value 8 accounts for time column and 7 property columns

    def __solute_solvent_clusters(self):
        for cls in self.clusters.values():
            res = set(cls.resnames)
            if (
                len(res.intersection(self.solute_resnames)) > 0
                and len(res.intersection(self.solvent_resnames)) > 0
            ):
                yield cls

    def __solute_solvent_analysis(self, frame):
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

    def __start_clusters(self):
        for subconn in self.conntab.subconntables():
            cls_id = self.__create_new_cluster(subconn)

            for mol in subconn:
                self.mol_clt[mol] = cls_id

    def __gen_origin_cluster_counter(self, subconn):
        mols_origin_clusters = {mol: self.mol_clt.get(mol, 0) for mol in subconn}
        return Counter(mols_origin_clusters.values())

    def __check_dominance(self, cls_id, n, conn_info, conn_skip):
        if self.clusters[cls_id].size == n:
            return True

        for subconn, count in conn_info[conn_skip + 1 :]:
            if len(subconn) == 2:
                return False
            if count[cls_id] > n:
                return False

        return True

    def __create_new_cluster(self, subconn):
        cluster = Cluster(self.uni, subconn)
        id = cluster.id
        self.clusters[id] = cluster

        return id

    def __construct_dominance(self, origin_clusters, modified_clusters, conn_info, i):
        dominances = {True: [], False: []}

        for cls_id, n in origin_clusters.most_common():
            if not cls_id:
                continue

            if n == 1:
                dominances[False].append(cls_id)
                continue

            if cls_id not in modified_clusters:
                dom = self.__check_dominance(cls_id, n, conn_info, i)
            else:
                dom = False

            dominances[dom].append(cls_id)

        return dominances

    def __get_older_cluster(self, dominances, origin_clusters):
        id = dominances[True][0]

        # this loop makes sure that in case the cluster that comes from merges of same
        # size agglomerates take the oldest one
        for j in dominances[True][1:]:
            if origin_clusters[id] == origin_clusters[j] and j < id:
                id = j
            else:
                break

        return id

    def __update_clusters(self):
        self.conntab.update()

        modified_mols = set()
        modified_clusters = set()
        merged_clusters = set()  # needs this because of dominance devolution

        conn_info = [
            (sub, self.__gen_origin_cluster_counter(sub)) for sub in self.conntab.subconntables()
        ]

        for i, (subconn, origin_clusters) in enumerate(conn_info):
            if len(origin_clusters) == 1:
                id = list(origin_clusters)[0]
                if id == 0:  # cluster formation
                    id = self.__create_new_cluster(subconn)
                elif id in modified_clusters:  # cluster separation
                    id = self.__create_new_cluster(subconn)
                else:
                    self.clusters[id].update_from_conntable(subconn)
            elif len(subconn) == 2:  # dimer is always new
                id = self.__create_new_cluster(subconn)
            else:
                dominances = self.__construct_dominance(
                    origin_clusters, modified_clusters, conn_info, i
                )

                if not dominances[True]:
                    id = self.__create_new_cluster(subconn)
                else:
                    id = self.__get_older_cluster(dominances, origin_clusters)

                    self.clusters[id].update_from_conntable(subconn)
                    merged_clusters.update(set(dominances[True]) - {id})

            for mol in subconn:
                self.mol_clt[mol] = id

            modified_clusters.add(id)
            modified_mols.update(subconn)

        for mol in set(self.mol_clt.keys()).difference(modified_mols):
            self.mol_clt.pop(mol)

        # TODO: deal with clusters that weren't modified. Needs to consider that some clusters merged (for log filing)  # noqa: E501
        for cls in set(self.clusters.keys()).difference(modified_clusters):
            self.clusters.pop(cls)

    def __write_coordinates(self):
        for cls in self.clusters.values():
            if not set(self.config["solute"]).intersection(set(cls.resnames)):
                continue

            with mda.Writer("tmp.gro", multiframe=False) as w:
                w.write(cls.ag.atoms.sort())

            with open(f"cls-n{cls.size}.gro", "a+") as out:
                with open("tmp.gro", "r") as tmp:
                    dt = tmp.readlines()

                dt[0] = f"Cluster-{cls.id} - Time = {self.uni.coord.time}\n"

                out.write("".join(dt))

    def run(self):
        if self.config.get("solute", False):
            self.__start_solute_solvent()
            self.__solute_solvent_analysis(0)

        with tqdm(
            total=len(self.uni.trajectory[1:]), initial=1, mininterval=5, miniters=10
        ) as pbar:
            for i, _ in enumerate(self.uni.trajectory[1:], start=1):
                self.__update_clusters()
                self.__get_clusters_info(i)
                if self.config.get("solute", False):
                    self.__solute_solvent_analysis(i)
                    self.__write_coordinates()

                if self.config.get("solute", False):
                    self.__nucleus_analysis(i)

                self.data_holder.parse_frame()
                pbar.update()

        # self.__print_clusters_index()
        np.savetxt("evo.txt", self.clusters_size_evo)
        if self.config.get("solute", False):
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

        if self.config.get("nucleus", False):
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

        with open("molclusters.json", "w+") as json_out:
            json.dump(self.data_holder.data, json_out, indent=4)

    def find(self, mol: int) -> Union[int, bool]:
        return self.mol_clt.get(mol, False)

    def __get_clusters_info(self, k):
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

    def __print_clusters_index(self):
        cols = 15
        i = 1
        with open("clusters_index.ndx", "w+", encoding="utf-8") as ndx:
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


class MolClustersData:
    def __init__(self, molclusters: MolClusters):
        self.molcls = molclusters
        self.data = {
            "Software": f"MolClusters {__version__}",
            "Trajectory": str(path.Path(self.molcls.uni.trajectory.filename).absolute()),
            "Topology": str(path.Path(self.molcls.uni.filename).absolute()),
            "Config": self.molcls.config,
            "MolClusters": [],
        }

    def parse_frame(self):
        data = {}

        data["Time"] = self.molcls.uni.coord.time
        data["Frame"] = self.molcls.uni.coord.frame
        data["NClusters"] = len(self.molcls.clusters)

        molclusters_data = []
        for cid, cls in self.molcls.clusters.items():
            cls_data = self.encode_cluster(cls)
            if self.molcls.config.get("nucleus", False):
                cls_data["Nucleus"] = []
                if self.molcls.nucleus_holder.get(cid, False):
                    for nuc in self.molcls.nucleus_holder[cid]:
                        cls_data["Nucleus"].append(MolClustersData.encode_nucleus(nuc))

            molclusters_data.append(cls_data)

        data["Clusters"] = molclusters_data
        self.data["MolClusters"].append(data)

    @staticmethod
    def encode_cluster(cls: Cluster):
        data = {}
        data["ID"] = cls.id
        MolClustersData.encode_properties(cls, data)
        return data

    @staticmethod
    def encode_nucleus(nuc: MDAAtomGroupAnalyzer):
        data = {}
        MolClustersData.encode_properties(nuc, data)
        return data

    @staticmethod
    def encode_properties(obj: Union[Cluster, MDAAtomGroupAnalyzer], data: dict):
        data["Size"] = obj.size
        data["Composition"] = MolClustersData.encode_composition(obj)
        data["ResIDs"] = sorted([int(x) for x in obj.resids])
        data["Mass"] = obj.mass
        data["Volume"] = obj.volume
        data["Density"] = obj.density
        data["Charge"] = obj.charge
        data["Dipole Moment"] = obj.dipole_moment
        data["Sphericity"] = obj.sphericity
        data["Shape"] = obj.shape_parameter

    @staticmethod
    def encode_composition(obj: Union[Cluster, MDAAtomGroupAnalyzer]):
        comp = {}
        for rnm, rid in zip(obj.resnames, map(int, obj.resids)):
            if rnm in comp:
                comp[rnm]["n"] += 1
                comp[rnm]["resids"].append(rid)
            else:
                comp[rnm] = {"resname": rnm, "n": 1, "resids": [rid]}

        return list(comp.values())
