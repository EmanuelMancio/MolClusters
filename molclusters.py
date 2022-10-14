#!python

"""Cluster Analyzer Script."""

import argparse as arg
from typing import Dict, Iterable, Iterator, List, Tuple, Type

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
        self._cluster = nx.compose(self._cluster, other.graph)

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
    def graph(self) -> Type[nx.Graph]:
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


def get_pair_and_distances(
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
        pairs, distances = mda.lib.distances.capped_distance(cm1, cm2, cutoff, box=box)
    else:
        pairs, distances = mda.lib.distances.self_capped_distance(cm1, cutoff, box=box)

    return pairs, distances


def init_clusters(selections, cluster_args, box, initial_time):
    clusters = []
    conn_tab = {}
    r_cluster = {}
    for resi in cluster_args:
        seli = selections[resi]
        for resj in cluster_args[resi]:
            selj = selections[resj]

            pairs, distances = get_pair_and_distances(
                seli, selj, cluster_args[resi][resj][0], box
            )

            clusters_from_selection(
                clusters,
                conn_tab,
                r_cluster,
                pairs,
                distances,
                seli,
                selj,
                initial_time,
            )

    clusters = {cl.id: cl for cl in clusters}

    return clusters, conn_tab, r_cluster


def clusters_from_selection(
    clusters: List[Type[Cluster]],
    conn_tab,
    r_cluster,
    pairs: Type[np.ndarray],
    distances: Type[np.ndarray],
    seli: Type[core.groups.AtomGroup],
    selj: Type[core.groups.AtomGroup],
    time: float,
) -> None:
    for k, [i, j] in enumerate(pairs):
        resi = seli.residues[i].resid
        resj = selj.residues[j].resid
        ni = selj.residues[i].resname
        nj = selj.residues[j].resname
        try:
            ri_clusters = [resi in cls for cls in clusters]

            rj_clusters = [resj in cls for cls in clusters]
        except AttributeError:
            clusters.append(Cluster(resi, resj, ni, nj, distances[k], time))
            continue

        ri_in_clusters = any(ri_clusters)
        rj_in_clusters = any(rj_clusters)

        if not ri_in_clusters and not rj_in_clusters:
            clusters.append(Cluster(resi, resj, ni, nj, distances[k], time))
            cls_id = clusters[-1].id
        elif not ri_in_clusters and rj_in_clusters:
            ind = rj_clusters.index(True)
            clusters[ind].add_mol(resj, resi, ni, distances[k])
            cls_id = clusters[ind].id
        elif ri_in_clusters and not rj_in_clusters:
            ind = ri_clusters.index(True)
            clusters[ind].add_mol(resi, resj, nj, distances[k])
            cls_id = clusters[ind].id
        else:
            ri_cluster = r_cluster[resi] - 1
            rj_cluster = r_cluster[resj] - 1

            # Both molecules in the same cluster: just add connection
            if ri_cluster == rj_cluster:
                clusters[rj_cluster].add_con(resi, resj, distances[k])
                cls_id = ri_cluster + 1

            # Molecules in different clusters: needs to consider order of cluster
            # in the clusters list so that the id of clusters and r_cluster are correctly changed
            elif ri_cluster < rj_cluster:
                clusters[ri_cluster].merge(clusters[rj_cluster])
                cls_id = ri_cluster + 1
                for mol in clusters[rj_cluster]:
                    r_cluster[mol] = cls_id

                clusters.pop(rj_cluster)
                clusters[ri_cluster].add_con(resi, resj, distances[k])
                correct_clusters_index(clusters, r_cluster, rj_cluster)
            else:
                clusters[rj_cluster].merge(clusters[ri_cluster])
                cls_id = rj_cluster + 1
                for mol in clusters[ri_cluster]:
                    r_cluster[mol] = cls_id

                clusters.pop(ri_cluster)
                clusters[rj_cluster].add_con(resi, resj, distances[k])
                correct_clusters_index(clusters, r_cluster, ri_cluster)

        r_cluster[resi] = cls_id
        r_cluster[resj] = cls_id

        if resi in conn_tab:
            conn_tab[resi][resj] = distances[k]
        else:
            conn_tab[resi] = {resj: distances[k]}
        if resj in conn_tab:
            conn_tab[resj][resi] = distances[k]
        else:
            conn_tab[resj] = {resi: distances[k]}


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


def connection_table(
    selections: Dict[str, Type[core.groups.AtomGroup]],
    cluster_args: Dict[str, Dict[str, Tuple[float, str]]],
    box: Type[np.ndarray],
) -> None:
    conn_tab = {}
    for resi in cluster_args:
        seli = selections[resi]

        for resj in cluster_args[resi]:
            selj = selections[resj]

            pairs: Type[np.ndarray]
            distances: Type[np.ndarray]

            pairs, distances = get_pair_and_distances(
                seli, selj, cluster_args[resi][resj][0], box
            )

            for k, [i, j] in enumerate(pairs):
                ri = seli.residues[i].resid
                rj = selj.residues[j].resid

                if ri in conn_tab:
                    conn_tab[ri][rj] = distances[k]
                else:
                    conn_tab[ri] = {rj: distances[k]}

                if rj in conn_tab:
                    conn_tab[rj][ri] = distances[k]
                else:
                    conn_tab[rj] = {ri: distances[k]}

    return conn_tab


def get_clusters_info(clusters_size_evo, clusters, uni: Type[mda.Universe], k):
    sizes = [cls.size for cls in clusters]
    avg = np.average(sizes)
    min_size = min(sizes)
    max_size = max(sizes)
    time = uni.coord.time

    clusters_size_evo[k][0] = time
    clusters_size_evo[k][4] = len(clusters)
    clusters_size_evo[k][2] = min_size
    clusters_size_evo[k][3] = avg
    clusters_size_evo[k][4] = max_size


def merge_clusters(clusters, i, j, r_cluster):
    if i == j:
        return
    elif clusters[i].size >= clusters[j].size:
        clusters[i].merge(clusters[j])
        for mol in clusters[j]:
            r_cluster[mol] = i

        clusters[j] = clusters[i].id
        return j
    else:
        clusters[j].merge(clusters[i])
        for mol in clusters[i]:
            r_cluster[mol] = j

        clusters[i] = clusters[j].id
        return i


def get_cluster_index(clusters, start):
    if type(clusters[start]) != int:
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

    clusters, conn_tab, r_cluster = init_clusters(
        selections, cluster_args, uni.dimensions, uni.trajectory.time
    )

    print_clusters_index(universe, clusters)

    clusters_size_evo = np.zeros((len(uni.trajectory), 5))

    get_clusters_info(clusters_size_evo, clusters.values(), uni, 0)

    for i, conf in enumerate(universe.trajectory[1:], start=1):
        conn_tab = connection_table(selections, cluster_args, uni.dimensions)

        # order of operations:
        # 1° Remove connections
        # 2° Add connections
        # 3° Add molecules
        # 4° Cluster Separation
        # 5° Remove molecules
        # 6° Cluster Merge
        # 7° Cluster Formation

        merge = set()
        graphs_from_sep = []

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
                            merge.add(tuple(to_merge))

                        # add connections between molecules of the cluster
                        else:
                            clst.add_con(mol, con, conn_tab[mol][con])
                    else:
                        # add free molecules to cluster
                        resname = uni.residues[con - 1].resname
                        clst.add_mol(mol, con, resname, conn_tab[mol][con])
                        r_cluster[con] = id

                # change distances to new values
                for con in con_to_edit:
                    clst.set_dist(mol, con, conn_tab[mol][con])

            if clst.size > 1:
                graphs_from_sep += clst.separate()
            elif clst.size == 1:
                if clst_will_merge:
                    pass
                else:
                    for mol in clst:
                        r_cluster.pop(mol)
                    clusters.pop(id)
            else:
                clusters.pop(id)

        for graph in graphs_from_sep:
            # remove molecules
            if len(graph) == 1:
                r_cluster.pop(list(graph)[0])
            else:
                # generate new cluster from separation
                tmp_cls = Cluster.from_graph(graph, conf.time)
                clusters[tmp_cls.id] = tmp_cls

                for mol in tmp_cls:
                    r_cluster[mol] = tmp_cls.id

        merged_clusters = []

        # merge clusters
        for m, n in merge:
            m, cls = get_cluster_index(clusters, m)
            n, cls = get_cluster_index(clusters, n)

            merged = merge_clusters(clusters, m, n, r_cluster)
            if merged:
                merged_clusters.append(merged)

        for cls in merged_clusters:
            clusters.pop(cls)

        # Identify new clusters
        for mol in conn_tab:
            if mol not in r_cluster:
                mols_to_clus = list(conn_tab[mol])
                new_cluster = Cluster(
                    mol,
                    mols_to_clus[0],
                    universe.residues[mol].resname,
                    universe.residues[mols_to_clus[0]].resname,
                    conn_tab[mol][mols_to_clus[0]],
                    conf.time,
                )

                r_cluster[mol] = new_cluster.id

                for m in mols_to_clus[1:]:
                    new_cluster.add_mol(
                        mol, m, universe.residues[mol].resname, conn_tab[mol][m]
                    )
                    r_cluster[m] = new_cluster.id

                clusters[new_cluster.id] = new_cluster

        get_clusters_info(clusters_size_evo, clusters.values(), uni, i)

    np.savetxt(clusters_size_evo)


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

    args = parser.parse_args()

    cls_args = parse_input_file(args.inp)

    uni = mda.Universe(args.top, args.traj, in_memory_step=1000) # TODO: add in_memory_step as option on cmdline

    analyze_trajectory(uni, cls_args)
