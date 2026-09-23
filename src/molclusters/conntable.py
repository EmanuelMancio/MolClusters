# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Defines the `ConnectionTable` class, which manages molecular connectivity tables for analyzing molecular clusters. It uses NetworkX for graph-based operations and MDAnalysis for molecular dynamics trajectory analysis.

Classes:
--------
    - ConnectionTable: Represents a connectivity table for molecular clusters, allowing operations
      such as retrieving connections, subgraphs, and constructing connectivity graphs.

Dependencies:
-------------
    - MDAnalysis: For molecular dynamics trajectory and structure analysis.
    - NetworkX: For graph-based operations on molecular connectivity.
    - NumPy: For numerical computations.
    - HydrogenBondAnalysis: For analyzing hydrogen bonds between molecules.
"""

import warnings
from typing import Generator, Iterator, overload

import MDAnalysis as mda
import networkx as nx
import numpy as np
from loguru import logger
from MDAnalysis import core
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

from .config import Rule
from .symdict import SymmetricDict

# ConnectionTable drives HydrogenBondAnalysis one frame at a time via these
# undocumented methods (see _check_hb_private_api below), plus a third: it
# also assigns hb._ts directly (read internally by _single_frame() as "the
# current frame") instead of letting run() iterate the trajectory itself,
# which isn't compatible with this per-frame external loop. _ts isn't
# checked below since it doesn't exist as an attribute until that
# assignment happens — a hasattr check on it would always fail regardless
# of whether the underlying mechanism still works.
_HB_PRIVATE_API = ("_prepare", "_single_frame")


def _check_hb_private_api(hb: HydrogenBondAnalysis) -> None:
    """Fail fast if MDAnalysis's private HydrogenBondAnalysis API has changed.

    These methods have no semver guarantee, so an MDAnalysis upgrade could
    silently rename or remove them. Checking eagerly, at construction time,
    turns that into a clear error instead of a cryptic AttributeError raised
    mid-analysis, potentially after a long-running trajectory has already
    been partially processed. This cannot catch every incompatibility (e.g.
    _single_frame changing what it reads instead of self._ts would silently
    produce wrong results, not an error) — only that these hooks still exist.

    Parameters
    ----------
    hb : HydrogenBondAnalysis
        The analysis instance to check.

    Raises
    ------
    RuntimeError
        If any of the private methods this class relies on are missing.
    """
    missing = [name for name in _HB_PRIVATE_API if not hasattr(hb, name)]
    if missing:
        logger.critical(
            f"This MDAnalysis version ({mda.__version__}) is incompatible with "
            "MolClusters' 'hb' rule support and the analysis cannot continue. "
            "This is a MolClusters bug, not something you did — please report "
            "it at https://github.com/EmanuelMancio/MolClusters/issues, "
            "including your MDAnalysis version."
        )
        raise RuntimeError(
            f"HydrogenBondAnalysis no longer exposes {missing} — this "
            "MDAnalysis version is incompatible with ConnectionTable's 'hb' "
            "rule support (see conntable.py)."
        )


class ConnectionTable:
    """Represents a connectivity table for molecular clusters.

    This class manages molecular connectivity data using NetworkX graphs and provides methods
    for retrieving connections, constructing subgraphs, and updating connectivity tables.

    Attributes
    ----------
    uni : MDAnalysis.Universe
        The MDAnalysis Universe object associated with the molecular system.
    clst_args : SymmetricDict[str, Rule]
        The clustering arguments specifying connectivity rules.
    sels : dict[str, core.groups.AtomGroup]
        The atom groups for each residue type.
    conntab : nx.MultiGraph
        The connectivity graph representing molecular connections.
    cms : dict[str, np.ndarray]
        The center of mass for each residue type.
    hbs : SymmetricDict[str, HydrogenBondAnalysis]
        The hydrogen bond analysis objects for residue pairs.
    """

    __slots__ = ["uni", "clst_args", "sels", "conntab", "cms", "hbs"]

    class _SubConnTable:
        """Represents a subgraph of the main connectivity table.

        Attributes
        ----------
        conntab : ConnectionTable
            The parent connectivity table.
        _graph : nx.MultiGraph
            The subgraph representing a subset of the connectivity table.
        _cm : np.ndarray
            The center of mass of the subgraph.
        """

        __slots__ = ["conntab", "_graph", "_cm", "_rg"]

        def __init__(
            self, subgraph: nx.MultiGraph, conntable: "ConnectionTable"
        ) -> None:
            """Initialize a subgraph of the connectivity table.

            Parameters
            ----------
            subgraph : nx.MultiGraph
                The subgraph representing a subset of the connectivity table.
                MultiGraph is used to allow multiple connections between the same nodes
            conntable : ConnectionTable
                The parent connectivity table.
            """
            self.conntab = conntable
            self._graph = subgraph
            self._rg = core.groups.ResidueGroup(
                [m - 1 for m in self.graph], self.conntab.uni
            )
            self._cm = self._rg.center_of_mass()

        def __getitem__(self, key: int) -> float | list[int]:
            """Get the attributes of a connection or molecule in the subgraph.

            Parameters
            ----------
            key : int
                The molecule or connection to retrieve.

            Returns
            -------
            float | list[int]
                The attributes of the connection or molecule.
            """
            return self.conntab[key]

        def __len__(self) -> int:
            """Get the number of molecules in the subgraph.

            Returns
            -------
            int
                The number of molecules in the subgraph.
            """
            return len(self._graph)

        @property
        def cm(self) -> np.ndarray:
            """The center of mass of the subgraph.

            Returns
            -------
            np.ndarray
                The center of mass of the subgraph.
            """
            return self._cm

        @property
        def graph(self) -> nx.Graph:
            """The graph representation of the subgraph.

            Returns
            -------
            nx.Graph
                The graph representation of the subgraph.
            """
            return self._graph

        @property
        def resnames(self) -> list[str]:
            """The residue names in the subgraph.

            Returns
            -------
            list[str]
                The residue names in the subgraph.
            """
            return self._rg.resnames

        def __iter__(self) -> Iterator[int]:
            """Iterate over the molecules in the subgraph.

            Returns
            -------
            Iterator[int]
                An iterator over the molecules in the subgraph.
            """
            return iter(self._graph)

    def __init__(
        self,
        universe: mda.Universe,
        cluster_args: SymmetricDict[str, Rule],
        selections: dict[str, core.groups.AtomGroup],
    ) -> None:
        """Initialize the ConnectionTable.

        Parameters
        ----------
        universe : mda.Universe
            The MDAnalysis Universe object associated with the molecular system.
        cluster_args : dict[str, dict[str, tuple[str, float]]]
            The clustering arguments specifying connectivity rules.
        selections : dict[str, core.groups.AtomGroup]
            The atom groups for each residue type.
        """
        self.uni = universe
        self.clst_args = cluster_args
        self.sels = selections

        self.cms: dict[str, np.ndarray] = {}
        self.hbs: SymmetricDict[str, HydrogenBondAnalysis] = SymmetricDict()
        self.__start_hbonds()
        self.update()

    def __get_mass_centers(self) -> None:
        """Calculate the center of mass for each residue type."""
        for res in self.clst_args.all_keys():
            self.cms[res] = self.sels[res].center_of_mass(compound="residues")

    def __start_hbonds(self) -> None:
        """Initialize hydrogen bond analysis for residue pairs."""
        for resi, resj in self.clst_args:
            if self.clst_args[resi, resj].type == "cm":
                continue

            if (resi, resj) not in self.hbs:
                # TODO: activate supported backend
                hb = HydrogenBondAnalysis(
                    self.uni,
                    between=[f"resname {resi}", f"resname {resj}"],
                    d_a_cutoff=self.clst_args[resi, resj].dist,
                    d_h_a_angle_cutoff=self.clst_args[resi, resj].ang,
                    update_selections=False,
                )
                _check_hb_private_api(hb)

                hb._prepare()

                self.hbs[resi, resj] = hb

    # TODO: break into two methods for cm and hb
    def __get_connections_and_attributes(
        self,
        resi: str,
        resj: str,
        box: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Private method to compute connections and their attributes between two residues.

        This method calculates the connections and associated attributes (e.g., distances, angles)
        between two residues (`resi` and `resj`) based on the provided cutoff criteria or hydrogen
        bonding information. The method supports both center-of-mass (CM) distance calculations
        and hydrogen bond (HB) analysis.

        Parameters
        ----------
            resi (str): The identifier for the first residue.
            resj (str): The identifier for the second residue.
            box (np.ndarray): The simulation box dimensions, used for periodic boundary conditions.

        Returns
        -------
            tuple[np.ndarray, np.ndarray]:
                - A 2D numpy array of connections, where each row represents a pair of residue IDs.
                - A list of dictionaries containing attributes for each connection, such as distance
                  and angle (if applicable).

        Notes
        -----
            - If the connection type is "cm" (center-of-mass), the method computes distances between
              the centers of mass of the residues.
            - If the connection type is "hb" (hydrogen bond), the method computes hydrogen bond
              distances and angles using precomputed hydrogen bond data.
            - The method updates internal state variables (e.g., `self.hbs`) to track progress
              through the hydrogen bond data.
            - Warnings are suppressed when no hydrogen bonds are found during the computation.
        """
        if self.clst_args[resi, resj].type == "cm":
            cm1: np.ndarray = self.cms[resi]
            cutoff = self.clst_args[resi, resj].dist
            if resi != resj:
                cm2: np.ndarray = self.cms[resj]
                connections, distances = mda.lib.distances.capped_distance(
                    cm1, cm2, cutoff, box=box
                )
            else:
                connections, distances = mda.lib.distances.self_capped_distance(
                    cm1, cutoff, box=box
                )

            for k, (moli, molj) in enumerate(connections):
                connections[k, 0] = self.sels[resi].residues[moli].resid  # noqa: B909
                connections[k, 1] = self.sels[resj].residues[molj].resid  # noqa: B909

            attributes = [{"distance": dist} for dist in distances]
        else:
            # TODO: implement own HB analysis as HydrogenBondAnalysis from mda repeats
            # distance and angle calculations. Until then, this depends on the private
            # API guarded by _check_hb_private_api (see its docstring for why).
            hb = self.hbs[resi, resj]
            hb._ts = self.uni.trajectory.ts

            # suppress warnings when there are no HBonds
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                hb._single_frame()

            distances = hb.results.hbonds[-2]
            angles = hb.results.hbonds[-1]

            attributes = [
                {"distance": dist, "angle": ang}
                for dist, ang in zip(distances, angles, strict=True)
            ]

            connections = np.empty((0, 2), int)
            for h_ati, a_ati in zip(
                hb.results.hbonds[2], hb.results.hbonds[3], strict=True
            ):
                h_ati, a_ati = int(h_ati), int(a_ati)
                moli = self.uni.atoms[h_ati].resid
                molj = self.uni.atoms[a_ati].resid
                connections = np.append(connections, [[moli, molj]], axis=0)

            hb._prepare()

        return connections, attributes

    def __construct_table(self) -> None:
        """Construct the connectivity table as a graph."""
        self.conntab = nx.Graph()

        for resi, resj in self.clst_args:
            pairs, attribs = self.__get_connections_and_attributes(
                resi, resj, self.uni.dimensions
            )

            for k, (ri, rj) in enumerate(pairs):
                self.conntab.add_edge(ri, rj, **attribs[k])

    def update(self) -> None:
        """Update the connectivity table."""
        self.__get_mass_centers()
        self.__construct_table()

    @overload
    def __getitem__(self, key: tuple[int, int]) -> int: ...

    @overload
    def __getitem__(self, key: int) -> list[int]: ...

    def __getitem__(self, key: tuple[int, int] | int) -> float | list[int]:
        """Get the attributes of a connection or molecule in the connectivity table.

        Parameters
        ----------
        key : tuple[int, int], int
            The molecule or connection to retrieve.

        Returns
        -------
        float, list[int]
            The attributes of the connection or molecule.

        Raises
        ------
        KeyError
            If the key is invalid or not found in the connectivity table.
        """
        if isinstance(key, tuple):
            if len(key) > 2:
                raise KeyError(
                    "ConnTable accepts only one or two parameters to access data!"
                )

            if key[0] not in self:
                raise KeyError(f"{key[0]} not in ConnTable")

            if key[1] not in self:
                raise KeyError(f"{key[1]} not in ConnTable")

            return self.conntab[key[0]][key[1]]["distance"]

        if key not in self:
            raise KeyError(f"{key} not in ConnTable")

        return self.mols_connected_to(key)

    def __contains__(self, item: int) -> bool:
        """Check if a molecule is in the connectivity table.

        Parameters
        ----------
        item : int
            The molecule to check.

        Returns
        -------
        bool
            True if the molecule is in the connectivity table, False otherwise.
        """
        return item in self.conntab

    def __graph_from_mol(self, mol: int) -> nx.Graph:
        desc = nx.node_connected_component(self.conntab, mol)
        return self.conntab.subgraph(desc)

    def __iter__(self) -> Iterator[int]:
        """Iterate over the molecules in the connectivity table.

        Returns
        -------
        Iterator[int]
            An iterator over the molecules in the connectivity table.
        """
        return iter(self.conntab)

    def connections_from(self, mol: int) -> list[tuple[int, int]]:
        """Get the connections from a molecule.

        Parameters
        ----------
        mol : int
            The molecule to retrieve connections from.

        Returns
        -------
        list[tuple[int, int]]
            A list of connections from the molecule.
        """
        return list(self.conntab.edges(mol))

    def connection_tree_from(self, mol: int) -> list[tuple[int, int]]:
        """Get the connection tree from a molecule.

        Parameters
        ----------
        mol : int
            The molecule to retrieve the connection tree from.

        Returns
        -------
        list[tuple[int, int]]
            A list of connections in the tree.
        """
        g = self.__graph_from_mol(mol)
        return list(g.edges)

    def mols_connected_to(self, mol: int) -> list[int]:
        """Get the molecules connected to a given molecule.

        Parameters
        ----------
        mol : int
            The molecule to retrieve connected molecules for.

        Returns
        -------
        list[int]
            A list of connected molecules.
        """
        return list(self.conntab[mol])

    def mols_connected_tree_to(self, mol: int) -> list[int]:
        """Get the molecules in the connection tree of a given molecule.

        Parameters
        ----------
        mol : int
            The molecule to retrieve the connection tree for.

        Returns
        -------
        list[int]
            A list of molecules in the connection tree.
        """
        g = self.__graph_from_mol(mol)
        return list(g.nodes)

    def __subgraphs(self) -> Generator[nx.Graph, None, None]:
        """Generate subgraphs of the connectivity table.

        Yields
        ------
        nx.Graph
            A subgraph of the connectivity table.
        """
        for c in nx.connected_components(self.conntab):
            yield self.conntab.subgraph(c)

    def subconntables(self) -> Generator["_SubConnTable", None, None]:
        """Generate sub-connectivity tables from the connectivity table.

        Yields
        ------
        _SubConnTable
            A sub-connectivity table.
        """
        # sorts to guarantee that in case of separation the biggest cluster keeps the id
        for s in sorted(self.__subgraphs(), key=lambda x: len(x), reverse=True):
            if len(s) == 1:
                # TODO: find a way to not change how HB calculations are done to avoid this check
                continue
            yield self._SubConnTable(s, self)
