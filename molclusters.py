#!python

"""Cluster Analyzer Script."""

import argparse as arg
import copy
import warnings
from collections import Counter
from typing import Dict, Iterable, Iterator, List, Optional, Tuple, Type, Union

# with warnings.catch_warnings(record=False,module='MDAnalysis',action="ignore"):
import MDAnalysis as mda

import networkx as nx
import numpy as np
import yaml
from MDAnalysis import core
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
# from MDAnalysis.transformations import center_in_box
from tqdm import tqdm


class Cluster:
    __cls_id = 1

    __slots__ = ["uni", "initial_time", "cluster", "_cm", "_id", "_ag","__centered_time"]

    def __init__(
        self,
        universe: Type[mda.Universe],
        subconntab: Optional[Type["ConnTable.SubConnTable"]] = None,
    ) -> None:
        self.uni = universe
        self.initial_time: float = universe.trajectory.time
        self.cluster = None

        if subconntab:
            self.cluster: Type[nx.Graph] = nx.Graph(subconntab.graph)
            self._cm = (
                subconntab._cm
            )  # TODO: Change to keep the sum of center of masses of molecules
        else:
            self.cluster = nx.Graph()
            self._cm = np.empty(3)

        self._ag = core.groups.ResidueGroup(np.array(self.cluster) - 1, self.uni)

        self._id = Cluster.__cls_id
        Cluster.__cls_id += 1
        self.__centered_time = -np.inf

    @classmethod
    def _from_graph(cls, uni: Type[mda.Universe], graph: Type[nx.Graph]) -> "Cluster":
        tmp_cls = cls(uni, None)
        tmp_cls.cluster = graph
        tmp_cls.__recalculate_cm()
        return tmp_cls

    @property
    def ag(self):
        return self._ag

    @property
    def id(self) -> int:
        return self._id

    @property
    def cm(self) -> Type[np.ndarray]:
        return self._cm

    @cm.setter
    def cm(self, value: Type[np.ndarray]):
        self._cm = value

    def __recalculate_cm(self):
        self._cm = self._ag.center_of_mass()

    def __update_ag(self) -> None:
        self._ag = core.groups.ResidueGroup(np.array(self.cluster) - 1, self.uni)

    def add_mol(
        self,
        ref_mol: int,
        mol: int,
        resname: str,
        dist: float,
    ) -> None:
        if ref_mol not in self:
            raise ValueError(f"ref_mol {ref_mol} not in the cluster")

        if mol in self:
            raise ValueError(f"mol {mol} already in the cluster, use add_con instead")

        self.cluster.add_node(mol, name=resname)
        self.add_con(ref_mol, mol, dist)
        self._ag += core.groups.ResidueGroup([mol - 1], self.uni)
        self.__recalculate_cm()

    def add_con(
        self,
        moli: int,
        molj: int,
        dist: float,
    ) -> None:
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster, use add_mol instead")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster, use add_mol instead")

        if dist == 0.0:
            raise ValueError("dist is zero. Check your trajectory")

        self.cluster.add_edge(moli, molj, distance=dist, weight=np.exp(1 / dist))

    def remove_mol(self, mol: int) -> None:
        if mol not in self:
            raise ValueError(f"mol {mol} not in the cluster")

        self.cluster.remove_node(mol)
        self._ag -= core.groups.ResidueGroup([mol - 1], self.uni)
        self.__recalculate_cm()

    def remove_con(self, moli: int, molj: int) -> None:
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster")

        self.cluster.remove_edge(moli, molj)

    def remove_cons(self, ref_mol: int, cons: Iterable) -> None:
        for con in cons:
            self.remove_con(ref_mol, con)

    def get_age(self, time: float) -> float:
        return time - self.initial_time

    def merge(self, other: "Cluster") -> None:
        self.cluster = nx.compose(self.cluster, other._graph)
        self.__recalculate_cm()

    def separate(self) -> List[Type[nx.Graph]]:
        sub_clusters: List[Type[nx.Graph]] = [
            self.cluster.subgraph(c).copy()
            for c in sorted(
                nx.connected_components(self.cluster), key=len, reverse=True
            )
        ]

        self.cluster = sub_clusters[0]
        self.__update_ag()
        self.__recalculate_cm()

        return [Cluster._from_graph(self.uni, nx.Graph(g)) for g in sub_clusters[1:]]

    @property
    def _graph(self) -> Type[nx.Graph]:
        return self.cluster

    @property
    def size(self) -> int:
        return len(self)

    def update_from_conntable(self, conn: Type["ConnTable.SubConnTable"]):
        # removed_mols = set(self).difference(conn)
        # removed_cons = self.cluster.edges - conn.graph.edges
        # new_cons = conn.graph.edges - self.cluster.edges

        self.cluster = nx.Graph(conn.graph)
        self._cm = conn.cm
        self.__update_ag()

    def get_dist(self, moli: int, molj: int) -> float:
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster")

        return self.cluster[moli][molj]["distance"]

    def set_dist(self, moli, molj, dist):
        if dist == 0.0:
            raise ValueError("dist is zero")

        self.cluster[moli][molj]["weight"] = np.exp(1 / dist)
        self.cluster[moli][molj]["distance"] = dist

    def get_weight(self, moli: int, molj: int) -> float:
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster")

        return self.cluster[moli][molj]["weight"]

    def __make_whole(self):
        #! This should NOT be used before printing
        if self.__centered_time != self.uni.trajectory.time:
            boxcenter = np.sum(self.uni.trajectory.ts.triclinic_dimensions, axis=0) / 2
            # self._ag.atoms.unwrap(compound="residues",reference="cog",inplace=True)

            ref_mol_cm = self._ag[:1].center_of_mass(unwrap=True)
            vector = boxcenter - ref_mol_cm
            self._ag.atoms.positions += vector
            self._ag.atoms.unwrap(compound="residues",reference="cog",inplace=True)
            # center_in_box(self._ag,point=ref_mol_cm)(self.uni.trajectory.ts)
            # center_in_box(self._ag)(self.uni.trajectory.ts)
            self.__centered_time = self.uni.trajectory.time

    @property
    def resnames(self):
        return self._ag.resnames

    @property
    def mass(self):
        return self._ag.total_mass()

    @property
    def sphericity(self):
        self.__make_whole()
        return 1 - self._ag.asphericity()

    @property
    def dipole_moment(self):
        self.__make_whole()
        return self._ag.atoms.dipole_moment()

    @property
    def dipole(self):
        self.__make_whole()
        return self._ag.atoms.dipole_vector()

    @property
    def shape_parameter(self):
        self.__make_whole()
        return self._ag.shape_parameter()

    @property
    def bsphere(self):
        self.__make_whole()
        return self._ag.bsphere()

    @property
    def radius_of_gyration(self):
        self.__make_whole()
        return self._ag.radius_of_gyration()

    @property
    def volume(self):
        r = self.radius_of_gyration
        return 4 * np.pi * r**3 / 3  # angstrom^3

    @property
    def density(self):
        return (self.mass / self.volume) * 0.602214076  # g/cm^3

    @property
    def charge(self):
        return self._ag.total_charge()
    def __contains__(self, item: int) -> bool:
        return item in self.cluster

    def __iter__(self) -> Iterator:
        return iter(self.cluster.copy())

    def __getitem__(self, key: int):
        return self.cluster[key]

    def __eq__(self, other) -> bool:
        if not isinstance(other, Cluster):
            return False

        return nx.utils.graphs_equal(self.cluster, other.cluster)

    def __len__(self) -> int:
        return len(self.cluster)

    def __str__(self):
        return self.cluster.edges.data().__str__()


class ConnTable:
    __slots__ = ["uni", "clst_args", "sels", "conntab", "cms", "hbs"]

    class SubConnTable:
        __slots__ = ["conntab", "_graph", "_cm"]

        def __init__(
            self, subgraph: Type[nx.Graph], conntable: Type["ConnTable"]
        ) -> None:
            self.conntab = conntable
            self._graph = subgraph

            self._cm = self.__selection().center_of_mass()

        def __selection(self) -> Type[core.groups.ResidueGroup]:
            return core.groups.ResidueGroup(
                [m - 1 for m in self.graph], self.conntab.uni
            )

        def __getitem__(self, key):
            return self.conntab[key]

        def __len__(self):
            return len(self._graph)

        @property
        def cm(self) -> Type[np.ndarray]:
            return self._cm

        @property
        def graph(self):
            return self._graph

        def __iter__(self):
            return iter(self._graph)

    def __init__(
        self,
        universe: Type[mda.Universe],
        cluster_args: Dict[str, Dict[str, Tuple[str, float]]],
        selections: Dict[str, Type[core.groups.AtomGroup]],
    ) -> None:
        self.uni = universe
        self.clst_args = cluster_args
        self.sels = selections

        self.cms = {}
        self.hbs: Dict[str, Dict[str, Type[HydrogenBondAnalysis]]] = {}
        self.__get_hbonds()
        self.update()

    def __get_mass_centers(self):
        for res in self.clst_args:
            self.cms[res] = self.sels[res].center_of_mass(compound="residues")

    def __get_hbonds(self):
        for resi in self.clst_args:
            for resj in self.clst_args[resi]:
                if self.clst_args[resi][resj] == "cm":
                    continue

                if (
                    (resi not in self.hbs)
                    or (resj not in self.hbs[resi])
                    or (resj not in self.hbs)
                    or (resi not in self.hbs[resj])
                ):
                    hb = HydrogenBondAnalysis(
                        self.uni,
                        between=[f"resname {resi}", f"resname {resj}"],
                        d_a_cutoff=self.clst_args[resi][resj][1]["d"],
                        d_h_a_angle_cutoff=self.clst_args[resi][resj][1]["a"],
                        update_selections=False,
                    )

                    hb._prepare()

                    if resi not in self.hbs:
                        self.hbs[resi] = {resj: [hb, 0]}
                    else:
                        self.hbs[resi][resj] = [hb, 0]

                    if resj not in self.hbs:
                        self.hbs[resj] = {resi: [hb, 0]}
                    else:
                        self.hbs[resj][resi] = [hb, 0]

    def __get_pair_and_distances(
        self,
        resi: str,
        resj: str,
        box: Type[np.ndarray],
    ) -> Tuple[Type[np.ndarray], Type[np.ndarray]]:
        pairs: Type[np.ndarray]
        distances: Type[np.ndarray]

        if self.clst_args[resi][resj][0] == "cm":
            cm1: Type[np.ndarray] = self.cms[resi]
            cutoff = self.clst_args[resi][resj][1]
            if resi != resj:
                cm2: Type[np.ndarray] = self.cms[resj]
                pairs, distances = mda.lib.distances.capped_distance(
                    cm1, cm2, cutoff, box=box
                )
            else:
                pairs, distances = mda.lib.distances.self_capped_distance(
                    cm1, cutoff, box=box
                )

            for k, (moli, molj) in enumerate(pairs):
                pairs[k, 0] = self.sels[resi].residues[moli].resid
                pairs[k, 1] = self.sels[resj].residues[molj].resid
        else:
            hb = self.hbs[resi][resj][0]
            hb._ts = self.uni.trajectory.ts

            # suppress warnings when there are no HBonds
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                hb._single_frame()

            frame_id = self.hbs[resi][resj][1]
            res = (np.asarray(hb.results.hbonds).T)[frame_id:, -4:-1]
            distances = res[:, -2]
            pairs = np.empty((0, 2), int)
            for h_ati, a_ati, _ in res:
                h_ati, a_ati = int(h_ati), int(a_ati)
                moli = self.uni.atoms[h_ati].resid
                molj = self.uni.atoms[a_ati].resid
                pairs = np.append(pairs, [[moli, molj]], axis=0)

            self.hbs[resi][resj][1] += len(res)

        return pairs, distances

    def __construct_table(self) -> Type[nx.Graph]:
        self.conntab = nx.Graph()

        analyzed = dict.fromkeys(self.clst_args.keys())
        for k in analyzed:
            analyzed[k] = []

        for resi in self.clst_args:
            for resj in self.clst_args[resi]:
                if resj not in analyzed[resi]:
                    pairs: Type[np.ndarray]
                    distances: Type[np.ndarray]

                    pairs, distances = self.__get_pair_and_distances(
                        resi, resj, self.uni.dimensions
                    )

                    for k, (ri, rj) in enumerate(pairs):
                        self.conntab.add_edge(ri, rj, d=distances[k])

                    analyzed[resi].append(resj)
                    analyzed[resj].append(resi)

    def update(self) -> None:
        self.__get_mass_centers()
        self.__construct_table()

    def __getitem__(self, key: Union[Tuple[int, int], int]) -> Union[float, List[int]]:
        if isinstance(key, tuple):
            if len(key) > 2:
                raise KeyError(
                    "ConnTable accepts only one or two parameters to access data!"
                )

            if key[0] not in self:
                raise IndexError(f"{key[0]} not in ConnTable")

            if key[1] not in self:
                raise IndexError(f"{key[1]} not in ConnTable")

            return self.conntab[key[0]][key[1]]["d"]

        if key not in self:
            raise IndexError(f"{key} not in ConnTable")

        return self.mols_connected_to(key)

    def __contains__(self, item: int) -> bool:
        return item in self.conntab

    def __graph_from_mol(self, mol):
        desc = nx.node_connected_component(self.conntab, mol)
        return self.conntab.subgraph(desc)

    def __iter__(self):
        return iter(self.conntab)

    def connections_from(self, mol):
        return list(self.conntab.edges(mol))

    def connection_tree_from(self, mol):
        g = self.__graph_from_mol(mol)
        return list(g.edges)

    def mols_connected_to(self, mol):
        return list(self.conntab[mol])

    def mols_connected_tree_to(self, mol):
        g = self.__graph_from_mol(mol)
        return list(g.nodes)

    def __subgraphs(self):
        for c in nx.connected_components(self.conntab):
            yield self.conntab.subgraph(c)

    def subconntables(self):
        for s in sorted(
            self.__subgraphs(), key=lambda x: len(x), reverse=True
        ):  # sorts to guarantee that in case of separation the biggest cluster keeps
            # the id
            yield self.SubConnTable(s, self)


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
        self.solute_radius = None
        self.solute_dipole = None

        self.__start_clusters()
        self.__get_clusters_info(0)

    def __start_solute_solvent(self):
        self.solutes = [
            id for sel in self.config["solute"] for id in self.sels[sel].residues.resids
        ]
        self.solvents = self.config[
            "solvent"
        ]  # TODO: solvents should be automatically discovered but file takes precedent

        self.solute_data = np.zeros((len(self.uni.trajectory), len(self.solutes)))
        self.solute_radius = np.zeros((len(self.uni.trajectory), len(self.solutes)))
        self.solute_dipole = np.zeros((len(self.uni.trajectory), len(self.solutes)))
        self.solute_density = np.zeros((len(self.uni.trajectory), len(self.solutes)))
        self.solute_sphericity = np.zeros((len(self.uni.trajectory), len(self.solutes)))
        self.solute_shape = np.zeros((len(self.uni.trajectory), len(self.solutes)))
        self.solute_charge = np.zeros((len(self.uni.trajectory), len(self.solutes)))
        self.solute_n_cluster = np.zeros(len(self.uni.trajectory))
        # self.solute_hb = {}
        # self.solvent_hb = {}

    def __solute_solvent_analysis(self, frame):
        all_solute_clusters = []
        for i, solute_id in enumerate(self.solutes):
            clst_id = self.find(solute_id)

            n_solvents = 0
            radius = 0
            dipole = np.nan
            density = 0
            sphericity = np.nan
            shape = np.nan
            charge = np.nan
            if clst_id:
                mol_pop = Counter(self.clusters[clst_id].resnames)
                for solvent in self.solvents:
                    n_solvents += mol_pop.get(solvent, 0)

                radius = self.clusters[clst_id].radius_of_gyration
                dipole = self.clusters[clst_id].dipole_moment
                density = self.clusters[clst_id].density
                sphericity = self.clusters[clst_id].sphericity
                shape = self.clusters[clst_id].shape_parameter
                charge = self.clusters[clst_id].charge

                all_solute_clusters.append(clst_id)

                # solvent_ids = set(self.clusters[clst_id])
                # solvent_ids.remove(solute_id)
                # solvent_selection = ""
                # for id in solvent_ids:
                #     solvent_selection += f"resid {id} or"

                # solvent_selection = solvent_selection[:-2]

                # # self.solute_hb[frame] = {}
                # # self.solvent_hb[frame] = {}

                # # TODO: fix for when the aggregate has more than one solute molecules
                # hb_solute = HydrogenBondAnalysis(
                #     self.uni,
                #     between=[f"resid {solute_id}", solvent_selection],
                #     d_a_cutoff=3.5,
                #     update_selections=False,
                # )
                # hb_solvent = HydrogenBondAnalysis(
                #     self.uni,
                #     between=[solvent_selection, solvent_selection],
                #     d_a_cutoff=3.5,
                #     update_selections=False,
                # )

                # hb_solute._prepare()
                # hb_solvent._prepare()

                # hb_solute._ts = self.uni.trajectory.ts
                # hb_solvent._ts = self.uni.trajectory.ts

                # with warnings.catch_warnings():
                #     warnings.simplefilter("ignore")
                #     hb_solute._single_frame()
                #     hb_solvent._single_frame()

                # hb_solute._conclude()
                # hb_solvent._conclude()

                # # self.solute_hb[frame][clst_id] = hb_solute.results
                # # self.solvent_hb[frame][clst_id] = hb_solvent.results

                # if len(hb_solute.results.hbonds) > 0:
                #     with open(f"{solute_id}_solute_hb.txt","a") as f:
                #         np.savetxt(f,hb_solute.results.hbonds)

                # if len(hb_solvent.results.hbonds) > 0:
                #     with open(f"{solute_id}_solvent_hb.txt","a") as f:
                #         np.savetxt(f,hb_solvent.results.hbonds)

            self.solute_data[frame][i] = n_solvents
            self.solute_radius[frame][i] = radius
            self.solute_dipole[frame][i] = dipole
            self.solute_density[frame][i] = density
            self.solute_sphericity[frame][i] = sphericity
            self.solute_shape[frame][i] = shape
            self.solute_charge[frame][i] = charge

        self.solute_n_cluster[frame] = len(set(all_solute_clusters))

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
            (sub, self.__gen_origin_cluster_counter(sub))
            for sub in self.conntab.subconntables()
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

            with mda.Writer("tmp.gro",multiframe=False) as w:
                w.write(cls.ag.atoms.sort())

            with open(f"cls-n{cls.size}.gro","a+") as out:
                with open("tmp.gro",'r') as tmp:
                    dt = tmp.readlines()

                dt[0] = f"Cluster-{cls.id} - Time = {self.uni.trajectory.time}\n"

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

                pbar.update()

        self.__print_clusters_index()
        np.savetxt("evo.txt", self.clusters_size_evo)
        np.savetxt("solute_solvent.txt", self.solute_data)
        np.savetxt("solute_radius.txt", self.solute_radius)
        np.savetxt("solute_dipole.txt", self.solute_dipole)
        np.savetxt("solute_density.txt", self.solute_density)
        np.savetxt("solute_sphericity.txt", self.solute_sphericity)
        np.savetxt("solute_shape.txt", self.solute_shape)
        np.savetxt("solute_charge.txt", self.solute_charge)
        np.savetxt("solute_n_clusters.txt", self.solute_n_cluster)

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

        for cls in self.clusters.values():
            if cls.id in self.radius_evolution:
                self.radius_evolution[cls.id].append((time, cls.radius_of_gyration))
            else:
                self.radius_evolution[cls.id] = [(time, cls.radius_of_gyration)]

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


def parse_input_file(
    in_file: Type[arg.FileType],
) -> Dict[str, Dict[str, Tuple[float, str]]]:
    config = yaml.safe_load(in_file)

    for mi, val in copy.deepcopy(config["rules"]).items():
        for mj, rule in val.items():
            rule = rule.split()
            op = rule[0].lower()
            if op == "cm":
                dist = float(rule[1])
                if mj in config["rules"]:
                    config["rules"][mj][mi] = (op, dist)
                else:
                    config["rules"][mj] = {mi: (op, dist)}
                config["rules"][mi][mj] = (op, dist)
            else:
                dist = float(rule[rule.index("d") + 1]) if "d" in rule else 3.5
                ang = float(rule[rule.index("a") + 1]) if "a" in rule else 150.0

                if mj in config["rules"]:
                    config["rules"][mj][mi] = (op, {"d": dist, "a": ang})
                else:
                    config["rules"][mj] = {mi: (op, {"d": dist, "a": ang})}
                config["rules"][mi][mj] = (op, {"d": dist, "a": ang})

    # TODO: Implement input correctness analysis

    return config


if __name__ == "__main__":
    parser = arg.ArgumentParser()

    parser.add_argument("traj", type=str, help="Trajectory File")
    parser.add_argument("top", type=str, help="Topology file")
    parser.add_argument(
        "inp", type=arg.FileType("r"), help="Input file with analysis settings"
    )

    args = parser.parse_args()

    cls_args = parse_input_file(args.inp)

    uni = mda.Universe(
        args.top, args.traj, in_memory=True
    )  # TODO: add in_memory_step as option on cmdline

    molclusters = MolClusters(uni, cls_args)
    molclusters.run()
