#!python

""" Cluster Analyzer Script """

import argparse as arg
from typing import Iterator, Type, List

import MDAnalysis as mda
from MDAnalysis import core
import networkx as nx
import numpy as np


class Cluster:
    def __init__(
        self,
        moli: Type[core.groups.Residue],
        molj: Type[core.groups.Residue],
        initial_time: float,
        dist: float,
    ) -> None:
        self._cluster: Type[nx.Graph] = nx.Graph()
        self._cluster.add_node(moli.resid, name=moli.resname)
        self._cluster.add_node(molj.resid, name=molj.resname)
        self._cluster.add_edge(moli.resid, molj.resid, weight=dist, distance=dist)
        self.initial_time: float = initial_time

    def add_mol(
        self, ref_mol: Type[core.groups.Residue], mol: Type[core.groups.Residue]
    ) -> None:
        if ref_mol.resid not in self:
            raise ValueError(f"ref_mol {ref_mol.resid} not in the cluster")

        if mol.resid in self:
            raise ValueError(
                f"mol {mol.resid} already in the cluster, use add_con instead"
            )

        self._cluster.add_node(mol.resid, name=mol.resname)
        self._cluster.add_edge(ref_mol.resid, mol.resid)

    def add_con(self, moli: int, molj: int, dist: float) -> None:
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster, use add_mol instead")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster, use add_mol instead")

        self._cluster.add_edge(moli, molj, weight=dist, distance=dist)

    def remove_mol(self, mol: Type[core.groups.Residue]) -> None:
        if mol.resid not in self:
            raise ValueError(f"mol {mol.resid} not in the cluster")

        self._cluster.remove_node(mol.resid)

    def remove_con(
        self, moli: Type[core.groups.Residue], molj: Type[core.groups.Residue]
    ) -> None:
        if moli.resid not in self:
            raise ValueError(f"mol {moli.resid} not in the cluster")
        if molj.resid not in self:
            raise ValueError(f"mol {molj.resid} not in the cluster")

        self._cluster.remove_edge(moli.resid, molj.resid)

    def get_lifetime(self, time: float) -> float:
        return time - self.initial_time

    def __contains__(self, item) -> bool:
        return item in self._cluster

    def __iter__(self) -> Iterator:
        return iter(self._cluster)

    def __getitem__(self, key: int):
        return self._cluster[key]

    def __eq__(self, other: "Cluster") -> bool:
        if not isinstance(other, Cluster):
            return NotImplemented

        this_molids: List[int] = list(self._cluster)
        other_molids: List[int] = list(other)

        return this_molids.sort() == other_molids.sort()

    def __len__(self) -> int:
        return len(self._cluster)


def parse_pairs(
    pairs: Type[np.ndarray], sel1, sel2, time: float
) -> List[Type[Cluster]]:
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

    clusters: List[Type[Cluster]] = []

    sel1: Type[core.groups.AtomGroup] = universe.select_atoms(
        f"resname {cluster_residues[0]}"
    )
    sel2: Type[core.groups.AtomGroup] = universe.select_atoms(
        f"resname {cluster_residues[1]}"
    )

    for conf in universe.trajectory:
        cm1: Type[np.ndarray] = sel1.center_of_mass(compound="residues")

        pairs: Type[np.ndarray]
        distances: Type[np.ndarray]
        # ? Should I consider only packed clusters or segments are acceptable?
        if cluster_residues[0] != cluster_residues[1]:
            cm2: Type[np.ndarray] = sel2.center_of_mass(compound="residues")
            pairs, distances = mda.lib.distances.capped_distance(
                cm1, cm2, cutoff, box=universe.dimensions
            )
        else:
            pairs, distances = mda.lib.distances.self_capped_distance(
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
