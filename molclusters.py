#!python

"""Cluster Analyzer Script."""

import argparse as arg
from collections import Counter
from typing import Dict, Iterable, Iterator, List, Tuple, Type, Union

import MDAnalysis as mda
import networkx as nx
import numpy as np
from MDAnalysis import core


class Cluster:
    __cls_id = 1

    __slots__ = ["uni", "initial_time", "cluster", "_cm", "_id"]

    def __init__(
        self,
        universe: Type[mda.Universe],
        subconntab: Type["ConnTable.SubConnTable"],
    ) -> None:
        self.uni = universe
        self.initial_time: float = universe.trajectory.time

        if subconntab:
            self.cluster: Type[nx.Graph] = nx.Graph(subconntab.graph)
            self._cm = (
                subconntab._cm
            )  # TODO: Change to keep the sum of center of masses of molecules

        self._id = Cluster.__cls_id
        Cluster.__cls_id += 1

    @classmethod
    def _from_graph(cls, uni: Type[mda.Universe], graph: Type[nx.Graph]) -> "Cluster":
        tmp_cls = cls(uni, None)
        tmp_cls.cluster = graph
        tmp_cls.__recalculate_cm()
        return tmp_cls

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
        self._cm = self.to_atomgroup().center_of_mass()

    def to_atomgroup(self) -> Type[core.groups.ResidueGroup]:
        # TODO: modify to have a selection as variable and to not recalculate it every time
        return core.groups.ResidueGroup([m - 1 for m in self.cluster], self.uni)

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
    __slots__ = ["uni", "clst_args", "sels", "conntab", "cms"]

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
        cluster_args: Dict[str, Dict[str, Tuple[float, str]]],
        selections: Dict[str, Type[core.groups.AtomGroup]],
    ) -> None:
        self.uni = universe
        self.clst_args = cluster_args
        self.sels = selections

        self.cms = {}
        self.__get_mass_centers()
        self.__construct_table()

    def __get_mass_centers(self):
        for res in self.clst_args:
            self.cms[res] = self.sels[res].center_of_mass(compound="residues")

    def __get_pair_and_distances(
        self,
        resi: str,
        resj: str,
        box: Type[np.ndarray],
    ) -> Tuple[Type[np.ndarray], Type[np.ndarray]]:
        cm1: Type[np.ndarray] = self.cms[resi]
        cutoff = self.clst_args[resi][resj][0]

        pairs: Type[np.ndarray]
        distances: Type[np.ndarray]
        if resi != resj:
            cm2: Type[np.ndarray] = self.cms[resj]
            pairs, distances = mda.lib.distances.capped_distance(
                cm1, cm2, cutoff, box=box
            )
        else:
            pairs, distances = mda.lib.distances.self_capped_distance(
                cm1, cutoff, box=box
            )

        return pairs, distances

    def __construct_table(self) -> Type[nx.Graph]:
        self.conntab = nx.Graph()
        for resi in self.clst_args:
            for resj in self.clst_args[resi]:
                pairs: Type[np.ndarray]
                distances: Type[np.ndarray]

                pairs, distances = self.__get_pair_and_distances(
                    resi, resj, self.uni.dimensions
                )

                for k, (i, j) in enumerate(pairs):
                    ri = self.sels[resi].residues[i].resid
                    rj = self.sels[resj].residues[j].resid

                    self.conntab.add_edge(ri, rj, d=distances[k])

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
        ):  # sorts to guarantee that in case of separation the biggest cluster keeps the id
            yield self.SubConnTable(s, self)


class MolClusters:
    __slots__ = ["uni", "args", "sels", "conntab", "clusters", "mol_clt","clusters_size_evo"]

    def __init__(
        self,
        universe: Type[mda.Universe],
        cluster_args: Dict[str, Dict[str, Tuple[float, str]]],
    ) -> None:
        self.uni = universe
        self.args = cluster_args
        self.sels: Dict[str, Type[core.groups.AtomGroup]] = {
            res: self.uni.select_atoms(f"resname {res}") for res in cluster_args
        }

        self.conntab = ConnTable(self.uni, self.args, self.sels)
        self.clusters: Dict[int, Type[Cluster]] = {}
        self.mol_clt: Dict[int, int] = {}

        self.clusters_size_evo = np.zeros((len(uni.trajectory), 5))

        self.__start_clusters()
        self.__get_clusters_info(0)

    def __start_clusters(self):
        for subconn in self.conntab.subconntables():
            cls_id = self.__create_new_cluster(subconn)

            for mol in subconn:
                self.mol_clt[mol] = cls_id

    def __gen_origin_cluster_counter(self, subconn):
        mols_origin_clusters = {
            mol: self.mol_clt.get(mol, 0) for mol in subconn
        }
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

    def __construct_dominance(self,origin_clusters,modified_clusters,conn_info,i):
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

        # this loop makes sure that in case the cluster that comes from merges of same size agglomerates take the oldest one
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
        merged_clusters = set() # needs this because of dominance devolution

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
                dominances = self.__construct_dominance(origin_clusters,modified_clusters,conn_info,i)

                if not dominances[True]:
                    id = self.__create_new_cluster(subconn)
                else:
                    id = self.__get_older_cluster(dominances,origin_clusters)

                    self.clusters[id].update_from_conntable(subconn)
                    merged_clusters.update(set(dominances[True])-{id})

            for mol in subconn:
                self.mol_clt[mol] = id

            modified_clusters.add(id)
            modified_mols.update(subconn)

        for mol in set(self.mol_clt.keys()).difference(modified_mols):
            self.mol_clt.pop(mol)

        # TODO: deal with clusters that weren't modified. Needs to consider that some clusters merged
        for cls in set(self.clusters.keys()).difference(modified_clusters):
            self.clusters.pop(cls)

    def run(self):
        for i, _ in enumerate(self.uni.trajectory[1:],start=1):
            self.__update_clusters()
            self.__get_clusters_info(i)
            if i == 8:
                self.__print_clusters_index()

        np.savetxt("evo.txt",self.clusters_size_evo)

    def find(self, mol) -> Union[int, bool]:
        return self.mol_clt.get(mol, False)

    def __get_clusters_info(self, k):
        sizes = [cls.size for cls in self.clusters.values()]
        avg = np.average(sizes)
        min_size = min(sizes, default=0)
        max_size = max(sizes, default=0)
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


def parse_input_file(
    in_file: Type[arg.FileType],
) -> Dict[str, Dict[str, Tuple[float, str]]]:
    cls_args: Dict[str, Dict[str, Tuple[float, str]]] = {}

    for line in in_file:
        line_elements = line.split()
        try:
            if line_elements[0] in cls_args:
                cls_args[line_elements[0]][line_elements[1]] = (
                    float(line_elements[2]),
                    "cm",
                )
            else:
                cls_args[line_elements[0]] = {
                    line_elements[1]: (float(line_elements[2]), "cm")
                }

            if line_elements[1] in cls_args:
                cls_args[line_elements[1]][line_elements[0]] = (
                    float(line_elements[2]),
                    "cm",
                )
            else:
                cls_args[line_elements[1]] = {
                    line_elements[0]: (float(line_elements[2]), "cm")
                }
        except IndexError:
            pass

    # TODO: Implement input correctness analysis

    return cls_args


if __name__ == "__main__":
    parser = arg.ArgumentParser()

    parser.add_argument("traj", type=str, help="Trajectory File")
    parser.add_argument("top", type=str, help="Topology file")
    parser.add_argument(
        "inp", type=arg.FileType("r"), help="Input file with distances information"
    )

    args = parser.parse_args()

    cls_args = parse_input_file(args.inp)

    uni = mda.Universe(
        args.top, args.traj, in_memory=True
    )  # TODO: add in_memory_step as option on cmdline

    molclusters = MolClusters(uni,cls_args)
    molclusters.run()
