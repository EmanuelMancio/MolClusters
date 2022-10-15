#!python

"""Cluster Analyzer Script."""

import argparse as arg
from typing import Dict, Iterable, Iterator, List, Tuple, Type, Union

import MDAnalysis as mda
import networkx as nx
import numpy as np
from MDAnalysis import core


class Cluster:
    __cls_id = 1

    def __init__(
        self,
        moli: int,
        molj: int,
        resi: str,
        resj: str,
        dist: float,
        initial_time: float,
    ) -> None:
        self._cluster: Type[nx.Graph] = nx.Graph()
        self._cluster.add_node(moli, name=resi)
        self._cluster.add_node(molj, name=resj)
        self.add_con(moli, molj, dist)
        self._initial_time: float = initial_time

        self._id = Cluster.__cls_id
        Cluster.__cls_id += 1

    @classmethod
    def from_graph(cls, graph: Type[nx.Graph], time: float) -> "Cluster":
        tmp_cls = cls(0, 0, "", "", 0.0, time)
        tmp_cls._cluster = graph
        return tmp_cls

    @property
    def id(self) -> int:
        return self._id

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

        self._cluster.add_node(mol, name=resname)
        self.add_con(ref_mol, mol, dist)

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
            self._cluster.add_edge(moli, molj, distance=dist, weight=np.inf)
        else:
            self._cluster.add_edge(moli, molj, distance=dist, weight=np.exp(1 / dist))

    def remove_mol(self, mol: int) -> None:
        if mol not in self:
            raise ValueError(f"mol {mol} not in the cluster")

        self._cluster.remove_node(mol)

    def remove_con(self, moli: int, molj: int) -> None:
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster")

        self._cluster.remove_edge(moli, molj)

    def remove_cons(self, ref_mol: int, cons: Iterable) -> None:
        for con in cons:
            self.remove_con(ref_mol, con)

    def get_lifetime(self, time: float) -> float:
        return time - self._initial_time

    def merge(self, other: "Cluster") -> None:
        self._cluster = nx.compose(self._cluster, other._graph)

    def separate(self) -> List[Type[nx.Graph]]:
        sub_clusters: List[Type[nx.Graph]] = [
            self._cluster.subgraph(c).copy()
            for c in sorted(
                nx.connected_components(self._cluster), key=len, reverse=True
            )
        ]

        self._cluster = sub_clusters[0]

        return sub_clusters[1:]

    @property
    def _graph(self) -> Type[nx.Graph]:
        return self._cluster

    @property
    def size(self) -> int:
        return len(self)

    def get_dist(self, moli: int, molj: int) -> float:
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster")

        return self._cluster[moli][molj]["distance"]

    def set_dist(self, moli, molj, dist):
        if dist == 0.0:
            self._cluster[moli][molj]["weight"] = np.inf
        else:
            self._cluster[moli][molj]["weight"] = np.exp(1 / dist)

        self._cluster[moli][molj]["distance"] = dist

    def get_weight(self, moli: int, molj: int) -> float:
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster")

        return self._cluster[moli][molj]["weight"]

    def __contains__(self, item: int) -> bool:
        return item in self._cluster

    def __iter__(self) -> Iterator:
        return iter(self._cluster.copy())

    def __getitem__(self, key: int):
        return self._cluster[key]

    def __eq__(self, other) -> bool:
        if not isinstance(other, Cluster):
            return False

        return nx.utils.graphs_equal(self._cluster, other._cluster)

    def __len__(self) -> int:
        return len(self._cluster)

    def __str__(self):
        return self._cluster.edges.data().__str__()


class ConnTable:
    def __init__(self, universe, cluster_args, selections) -> None:
        self._uni = universe
        self._clst_args = cluster_args
        self._sels = selections
        self._conntab = self.__construct_table()

    def __get_pair_and_distances(
        seli: Type[core.groups.AtomGroup],
        selj: Type[core.groups.AtomGroup],
        cutoff: float,
        box: Type[np.ndarray],
    ) -> Tuple[Type[np.ndarray], Type[np.ndarray]]:
        cm1: Type[np.ndarray] = seli.center_of_mass(compound="residues")

        pairs: Type[np.ndarray]
        distances: Type[np.ndarray]
        if seli != selj:
            cm2: Type[np.ndarray] = selj.center_of_mass(compound="residues")
            pairs, distances = mda.lib.distances.capped_distance(
                cm1, cm2, cutoff, box=box
            )
        else:
            pairs, distances = mda.lib.distances.self_capped_distance(
                cm1, cutoff, box=box
            )

        return pairs, distances

    def __construct_table(self) -> Type[nx.Graph]:
        conn_tab = {}
        for resi in self._clst_args:
            seli = self._sels[resi]

            for resj in self._clst_args[resi]:
                selj = self._sels[resj]

                pairs: Type[np.ndarray]
                distances: Type[np.ndarray]

                pairs, distances = self.__get_pair_and_distances(
                    seli, selj, self._clst_args[resi][resj][0], self._uni.dimensions
                )

                for k, [i, j] in enumerate(pairs):
                    ri = seli.residues[i].resid
                    rj = selj.residues[j].resid

                    if ri in conn_tab:
                        conn_tab[ri][rj] = {"d": distances[k]}
                    else:
                        conn_tab[ri] = {rj: {"d": distances[k]}}

                    if rj in conn_tab:
                        conn_tab[rj][ri] = {"d": distances[k]}
                    else:
                        conn_tab[rj] = {ri: {"d": distances[k]}}

        return nx.Graph(conn_tab)

    def update(self) -> None:
        self._conntab = self.__construct_table()

    def __getitem__(self, key: Union[Tuple[int, int], int]):
        if isinstance(key, tuple):
            if len(key) > 2:
                raise KeyError("ConnTable only accpets two parameter: ConnTable[i,j]")

            return self._conntab[key[0]][key[1]]["d"]
        elif isinstance(key, int):
            return self.mols_connected_to(key)

        raise KeyError(f"{type(key)} not accepted as key for ConnTable")

    def __contains__(self, item: int) -> bool:
        return item in self._conntab

    def __graph_from_mol(self, mol):
        desc = nx.node_connected_component(self._conntab, mol)
        return self._conntab.subgraph(desc)

    def connections_from(self, mol):
        return list(self._conntab.edges(mol))

    def all_connections_from(self, mol):
        g = self.__graph_from_mol(mol)
        return list(g.edges)

    def mols_connected_to(self, mol):
        return list(self._conntab[mol])

    def all_mols_connected_to(self, mol):
        g = self.__graph_from_mol(mol)
        return list(g.nodes)

    def _subgraphs(self):
        for c in nx.connected_components(self._cluster):
            yield self._conntab.subgraph(c)


def init_clusters(universe, conn_tab: Type[ConnTable]):
    clusters = {}
    r_cluster = {}
    for subgraph in conn_tab._subgraphs():
        cons = list(subgraph.edges)
        moli, molj = cons[0]
        resi = universe.residues[moli - 1].resname
        resj = universe.residues[molj - 1].resname
        cls = Cluster(moli, molj, resi, resj, conn_tab[moli, molj], uni.coord.time)
        id = cls.id

        r_cluster[moli] = id
        r_cluster[molj] = id

        for moli, molj in cons[1:]:
            if moli in cls and molj not in cls:
                resj = universe.residues[molj - 1].resname
                cls.add_mol(moli, molj, resj, conn_tab[moli, molj])
                r_cluster[molj] = id
            elif moli not in cls and molj in cls:
                resi = universe.residues[moli - 1].resname
                cls.add_mol(molj, moli, resi, conn_tab[moli, molj])
                r_cluster[moli] = id
            else:
                cls.add_con(moli, molj, conn_tab[moli, molj])

        clusters[cls.id] = cls

    return clusters, r_cluster


def correct_clusters_index(clusters: List[Type[Cluster]], r_cluster, from_pos: int):
    # TODO: Really bad code, should change this in the future
    Cluster._Cluster__cls_id -= 1
    if len(clusters) == from_pos:
        return
    for cls in clusters[from_pos:]:
        cls._id -= 1
        i = cls.id
        for mol in cls:
            r_cluster[mol] = i


def get_clusters_info(clusters_size_evo, clusters, uni: Type[mda.Universe], k):
    sizes = [cls.size for cls in clusters]
    avg = np.average(sizes)
    min_size = min(sizes)
    max_size = max(sizes)
    time = uni.coord.time

    clusters_size_evo[k][0] = time
    clusters_size_evo[k][1] = len(clusters)
    clusters_size_evo[k][2] = min_size
    clusters_size_evo[k][3] = avg
    clusters_size_evo[k][4] = max_size


def merge_clusters(
    clusters: Dict[int, Type[Cluster]],
    conn_tab: Type[ConnTable],
    i,
    j,
    r_cluster,
    connections,
):
    if i == j:
        return

    if clusters[i].size >= clusters[j].size:
        clusters[i].merge(clusters[j])
        for mol in clusters[j]:
            r_cluster[mol] = i

        for mi, mj in connections:
            clusters[i].add_con(mi, mj, conn_tab[mi, mj])

        clusters[j] = clusters[i].id
        return j
    else:
        clusters[j].merge(clusters[i])
        for mol in clusters[i]:
            r_cluster[mol] = j

        for mi, mj in connections:
            clusters[j].add_con(mi, mj, conn_tab[mi, mj])

        clusters[i] = clusters[j].id
        return i


def get_cluster_index(clusters, start):
    if not isinstance(clusters[start], int):
        return start, []
    else:
        id, clst = get_cluster_index(clusters, clusters[start])
        clusters[start] = id
        clst.append(start)
        return id, clst


def analyze_trajectory(
    universe: Type[mda.Universe],
    cluster_args: Dict[str, Dict[str, Tuple[float, str]]],
) -> None:
    """Analyses the formation of clusters in the trajectory

    Parameters
    ----------
    universe : Type[mda.Universe]
        mdanalysis Universe object containing the simulation data
    cluster_residues : List[str]
        List of the residues to consider in the analysis
    cutoff : float
        cm-cm distance to be considered in the analyses
    """

    selections: Dict[str, Type[core.groups.AtomGroup]] = {
        res: universe.select_atoms(f"resname {res}") for res in cluster_args
    }

    conn_tab = ConnTable(universe, cluster_args, selections)

    clusters, r_cluster = init_clusters(universe, ConnTable)

    # print_clusters_index(universe, clusters)

    clusters_size_evo = np.zeros((len(uni.trajectory), 5))

    get_clusters_info(clusters_size_evo, clusters.values(), uni, 0)

    for i, conf in enumerate(universe.trajectory[1:], start=1):
        conn_tab.update()

        # order of operations:
        # 1° Remove connections
        # 2° Add connections
        # 3° Add molecules
        # 4° Cluster Separation
        # 5° Remove molecules
        # 6° Cluster Merge
        # 7° Cluster Formation

        merge = {}

        for id, clst in clusters.copy().items():
            clst_will_merge = False
            for mol in clst:
                mol_con = set(clst[mol])

                if mol in conn_tab:
                    new_con = set(conn_tab[mol])
                else:
                    clst.remove_mol(mol)
                    r_cluster.pop(mol)
                    continue

                con_to_add = new_con.difference(mol_con)
                con_to_rem = mol_con.difference(new_con)
                con_to_edit = mol_con.intersection(new_con)

                clst.remove_cons(mol, con_to_rem)  # remove connections

                for con in con_to_add:
                    if con in r_cluster:
                        # identify clusters to merge
                        if r_cluster[con] != id:
                            clst_will_merge = True
                            to_merge = [id, r_cluster[con]]
                            to_merge.sort()
                            to_merge = tuple(to_merge)
                            if to_merge in merge:
                                merge[to_merge].append((mol, con))
                            else:
                                merge[to_merge] = [(mol, con)]

                        # add connections between molecules of the cluster
                        else:
                            clst.add_con(mol, con, conn_tab[mol, con])
                    else:
                        # add free molecules to cluster
                        resname = uni.residues[con - 1].resname
                        clst.add_mol(mol, con, resname, conn_tab[mol, con])
                        r_cluster[con] = id

                # change distances to new values
                for con in con_to_edit:
                    clst.set_dist(mol, con, conn_tab[mol, con])

            if clst.size == 1:
                if clst_will_merge:
                    pass
                else:
                    for mol in clst:
                        r_cluster.pop(mol)
                    clusters.pop(id)
            elif clst.size == 0:
                clusters.pop(id)

        merged_clusters = []

        # merge clusters
        for to_merge in merge:
            m, n = to_merge
            m, cls = get_cluster_index(clusters, m)
            n, cls = get_cluster_index(clusters, n)

            merged = merge_clusters(
                clusters, conn_tab, m, n, r_cluster, merge[to_merge]
            )
            if merged:
                merged_clusters.append(merged)

        for cls in merged_clusters:
            clusters.pop(cls)

        for clst in clusters.copy():
            for graph in clusters[clst].separate():
                # remove molecules
                if len(graph) == 1:
                    r_cluster.pop(list(graph)[0])
                else:
                    # generate new cluster from separation
                    tmp_cls = Cluster.from_graph(graph, conf.time)
                    clusters[tmp_cls.id] = tmp_cls

                    for mol in tmp_cls:
                        r_cluster[mol] = tmp_cls.id

        # Identify new clusters
        for mol in conn_tab:
            if mol not in r_cluster:
                cons = conn_tab.all_connections_from(mol)
                mi, mj = cons[0]
                new_cluster = Cluster(
                    mi,
                    mj,
                    universe.residues[mi - 1].resname,
                    universe.residues[mj - 1].resname,
                    conn_tab[mi,mj],
                    conf.time,
                )

                r_cluster[mi] = new_cluster.id
                r_cluster[mj] = new_cluster.id

                for mi, mj in cons[1:]:
                    if mi in new_cluster and mj in new_cluster:
                        new_cluster.add_con(mi, mj, conn_tab[mi,mj])
                    elif mi in new_cluster:
                        new_cluster.add_mol(
                            mi,
                            mj,
                            universe.residues[mj - 1].resname,
                            conn_tab[mi, mj],
                        )
                        r_cluster[mj] = new_cluster.id
                    elif mj in new_cluster:
                        new_cluster.add_mol(
                            mj,
                            mi,
                            universe.residues[mi - 1].resname,
                            conn_tab[mi,mj],
                        )
                        r_cluster[mi] = new_cluster.id

                clusters[new_cluster.id] = new_cluster

        get_clusters_info(clusters_size_evo, clusters.values(), uni, i)

    np.savetxt("evo.txt", clusters_size_evo)


def print_clusters_index(uni, clusters):
    cols = 15
    i = 1
    with open("clusters_index.ndx", "w+", encoding="utf-8") as ndx:
        for cluster in clusters:
            ndx.write(f"[ CLS-{cluster} ]\n")
            for mol in sorted(clusters[cluster]):
                for at in uni.residues[mol - 1].atoms:
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

    # args = parser.parse_args()

    args = parser.parse_args(
        "data//met-mal//met-traj.pdb data//met-mal//met-mal.tpr data//met-mal//cls.in".split()
    )

    cls_args = parse_input_file(args.inp)

    uni = mda.Universe(
        args.top, args.traj, in_memory=True
    )  # TODO: add in_memory_step as option on cmdline

    analyze_trajectory(uni, cls_args)
