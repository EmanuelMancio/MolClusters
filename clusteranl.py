#!python

""" Cluster Analyzer Script """

import argparse as arg
from typing import Iterator, Type, List

import MDAnalysis as mda
import networkx as nx
import numpy as np


class Cluster:
    def __init__(self, moli: int, moli_name: str, molj: int, molj_name: str, initial_time: float) -> None:
        self._cluster: Type[nx.Graph] = nx.Graph()
        self._cluster.add_node(moli,name=moli_name)
        self._cluster.add_node(molj,name=molj_name)
        self._cluster.add_edge(moli, molj)
        self.initial_time: float = initial_time

    def add_mol(self, ref_mol: int, molid: str, mol_name: str) -> None:
        if ref_mol not in self:
            raise ValueError(f"ref_mol {ref_mol} not in the cluster")

        if molid in self:
            raise ValueError(f"mol {molid} already in the cluster, use add_con instead")

        self._cluster.add_node(molid,name=mol_name)
        self._cluster.add_edge(ref_mol, molid)

    def add_con(self, moli:int, molj:int) -> None:
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster, use add_mol instead")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster, use add_mol instead")

        self._cluster.add_edge(moli, molj)

    def remove_mol(self, molid):
        if molid not in self:
            raise ValueError(f"mol {molid} not in the cluster")

        self._cluster.remove_node(molid)

    def remove_con(self, moli, molj):
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster")

        self._cluster.remove_edge(moli, molj)

    def __contains__(self,item):
        return item in self._cluster

    def __iter__(self) -> Iterator:
        return iter(self._cluster)

    def __getitem__(self,key: int):
        return self._cluster[key]

    def __eq__(self, other: "Cluster"):
        if not isinstance(other,Cluster):
            return NotImplemented

        this_molids = list(self._cluster)
        other_molids = list(other)

        return this_molids.sort() == other_molids.sort()

    def __len__(self) -> int:
        return len(self._cluster)

def parse_pairs(pairs: Type[np.ndarray], sel1, sel2, time: float) -> List[Type[Cluster]]:
    # TODO: Implements parse_pairs
    return []


def analyze_trajectory(
    universe: Type[mda.Universe], cluster_residues: List[str], cutoff: float
) -> None:
    """Analyses the formation of clusters in the trajectory

    Args:
        universe (Type[mda.Universe]): MDAnalysis Universe object containing the simulation data
        cluster_residues (List[str]): List of the residues to consider in the analysis
        cutoff (float): cm-cm distance to be considered in the analyses
    """

    clusters : List[Type[Cluster]]= []

    sel1 = universe.select_atoms(f"resname {cluster_residues[0]}")
    sel2 = universe.select_atoms(f"resname {cluster_residues[1]}")

    for conf in universe.trajectory:
        cm1 = sel1.center_of_mass()

        # ? Should I consider only packed clusters or segments are acceptable?
        if cluster_residues[0] != cluster_residues[1]:
            cm2 = sel2.center_of_mass()
            pairs: np.ndarray = mda.lib.distances.capped_distance(
                cm1, cm2, cutoff, box=universe.dimensions
            )
        else:
            pairs: np.ndarray = mda.lib.distances.self_capped_distance(
                cm1, cutoff, box=universe.dimensions
            )

        clusters = parse_pairs(pairs, sel1, sel2, conf.time)

        # TODO: Implement cluster evolution analysis



if __name__ == "__main__":
    parser = arg.ArgumentParser()

    parser.add_argument("traj", type=str, help="Trajectory File")
    parser.add_argument("top", type=str, help="Topology file")
    parser.add_argument("cutoff", type=float)
    parser.add_argument("residues", type=str, nargs=2)

    args = parser.parse_args()

    clst_residues = args.residues

    uni = mda.Universe(args.top, args.traj)
