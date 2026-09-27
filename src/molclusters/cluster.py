# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides classes for analyzing molecular clusters using MDAnalysis and NetworkX.

Classes:
--------
    - MDAResidueGroupAnalyzer: A utility class for analyzing MDAnalysis AtomGroups, providing
      properties and methods to calculate various molecular properties such as mass, sphericity,
      dipole moment, radius of gyration, and density.

    - Cluster: A subclass of MDAResidueGroupAnalyzer that represents a molecular cluster. It extends
      the functionality to manage molecular clusters as graphs using NetworkX, allowing for
      operations such as adding/removing molecules and connections, merging clusters, separating
      sub-clusters, and updating clusters from connectivity tables.

Constants:
----------
    - EA2D: Conversion factor for dipole moment from atomic units to Debye.

Dependencies:
-------------
    - MDAnalysis: For molecular dynamics trajectory and structure analysis.
    - NetworkX: For graph-based operations on molecular clusters.
    - NumPy: For numerical computations.
"""

from contextlib import contextmanager
from functools import wraps
from typing import Callable, Iterable, Iterator, Self

import MDAnalysis as mda
import networkx as nx
import numpy as np
from MDAnalysis import core

from .conntable import ConnectionTable

EA2D = 1 / 0.3934303


def _on_whole[T](method: Callable[..., T]) -> Callable[..., T]:
    """Run `method` with the analyzer's atoms temporarily made whole (see `whole`).

    Returns
    -------
    Callable
        The wrapped method.
    """

    @wraps(method)
    def wrapper(self: "MDAResidueGroupAnalyzer", *args: object, **kwargs: object) -> T:
        with self.whole():
            return method(self, *args, **kwargs)

    return wrapper


class MDAResidueGroupAnalyzer:
    """A utility class for analyzing MDAnalysis ResidueGroups.

    This class provides properties and methods to calculate various molecular properties
    such as mass, sphericity, dipole moment, radius of gyration, and density.

    Attributes
    ----------
    uni : MDAnalysis.Universe
        The MDAnalysis Universe object associated with the AtomGroup.
    _rg : MDAnalysis.core.groups.ResidueGroup
        The ResidueGroup being analyzed.
    __whole_key : tuple[int, bytes] | None
        Frame and residue indices that `__whole_positions` was computed for.
    __whole_positions : np.ndarray | None
        Cached whole positions of the group's atoms (see `whole`).

    Methods
    -------
    Various properties and methods to compute molecular properties.
    """

    __slots__ = ["uni", "_rg", "__whole_key", "__whole_positions"]

    def __init__(
        self, universe: mda.Universe, residues: Iterable[int] | core.groups.ResidueGroup
    ) -> None:
        """Initialize the MDAResidueGroupAnalyzer with a universe and residues.

        Parameters
        ----------
        universe : mda.Universe
            The MDAnalysis Universe object associated with the ResidueGroup.
        residues : Union[Iterable[int], core.groups.ResidueGroup]
            The residues to be analyzed, either as an iterable with residue IDs (0-based) or a ResidueGroup.
        """
        self.uni = universe

        if isinstance(residues, core.groups.ResidueGroup):
            self._rg = residues
        else:
            self._rg = core.groups.ResidueGroup(np.array(residues) - 1, self.uni)

        self.__whole_key: tuple[int, bytes] | None = None
        self.__whole_positions: np.ndarray | None = None

    def __compute_whole_positions(self) -> np.ndarray:
        """Compute the group's atom positions with the group made whole.

        The group is translated so that its first residue sits at the box center,
        then each residue is wrapped back into the box by its center of geometry,
        which gathers the group around the center (assuming it spans less than half
        the box). The Universe's positions are left exactly as they were.

        Returns
        -------
        np.ndarray
            Whole positions, one row per atom of ``self._rg.atoms``.
        """
        atoms = self._rg.atoms
        original = atoms.positions  # a copy
        try:
            boxcenter = np.sum(self.uni.trajectory.ts.triclinic_dimensions, axis=0) / 2
            ref_mol_cm = self._rg[:1].center_of_mass(unwrap=True)
            atoms.positions += boxcenter - ref_mol_cm
            atoms.unwrap(compound="residues", reference="cog", inplace=True)
            return atoms.positions
        finally:
            atoms.positions = original

    @contextmanager
    def whole(self) -> Iterator[core.groups.AtomGroup]:
        """Temporarily place the group's atoms in whole, centered positions.

        Inside the ``with`` block the group's atoms are whole across periodic
        boundaries, so geometric properties (and writers) see the real shape. On
        exit, even on error, the atoms go back to their original positions.

        Positions are never changed for good because they belong to the shared
        Universe: other analyzers over overlapping residues (e.g. a nucleus inside a
        cluster) would otherwise move atoms under each other and corrupt each
        other's geometry.

        The whole positions are cached per frame and residue set, so the unwrap
        runs once per frame no matter how many properties are read.

        Yields
        ------
        core.groups.AtomGroup
            The group's atoms, in whole positions.
        """
        key = (self.uni.trajectory.ts.frame, self._rg.ix.tobytes())
        if self.__whole_key != key:
            self.__whole_positions = self.__compute_whole_positions()
            self.__whole_key = key

        atoms = self._rg.atoms
        original = atoms.positions  # a copy
        atoms.positions = self.__whole_positions
        try:
            yield atoms
        finally:
            atoms.positions = original

    def __add__(self, other: core.groups.ResidueGroup | Self) -> Self:
        """Combine this ResidueGroupAnalyzer with another ResidueGroup or ResidueGroupAnalyzer.

        Parameters
        ----------
        other : Union[core.groups.ResidueGroup, Self]
            The other ResidueGroup or ResidueGroupAnalyzer to combine with.

        Returns
        -------
        Self
            A new MDAResidueGroupAnalyzer instance representing the combined ResidueGroups.
        """
        if isinstance(other, core.groups.ResidueGroup):
            return MDAResidueGroupAnalyzer(self.uni, self._rg + other)
        elif isinstance(other, MDAResidueGroupAnalyzer):
            return MDAResidueGroupAnalyzer(self.uni, self._rg + other.ag)

        return NotImplemented

    @property
    def ag(self) -> core.groups.ResidueGroup:
        """The ResidueGroup being analyzed.

        Returns
        -------
        core.groups.ResidueGroup
            The ResidueGroup associated with this analyzer.
        """
        return self._rg

    @property
    def size(self) -> int:
        """The size of the ResidueGroup.

        Returns
        -------
        int
            The number of residues in the ResidueGroup.
        """
        return len(self)

    @property
    def resnames(self) -> list[str]:
        """The residue names of the ResidueGroup.

        Returns
        -------
        list[str]
            A list of residue names in the ResidueGroup.
        """
        return self._rg.resnames

    @property
    def resids(self) -> list[int]:
        """The residue IDs of the ResidueGroup.

        Returns
        -------
        list[int]
            A list of residue IDs in the ResidueGroup.
        """
        return self._rg.resids

    @property
    def mass(self) -> float:
        """Calculate the total mass of the ResidueGroup.

        Returns
        -------
        float
            The total mass of the ResidueGroup.
        """
        return self._rg.total_mass()

    @property
    @_on_whole
    def sphericity(self) -> float:
        """Calculate the sphericity of the ResidueGroup.

        Returns
        -------
        float
            The sphericity of the ResidueGroup, a measure of how spherical the shape is.
            A value near zero represents a spherical shape.
        """
        return 1 - self._rg.asphericity()

    @property
    @_on_whole
    def dipole_moment(self) -> float:
        """Calculate the dipole moment of the ResidueGroup.

        Returns
        -------
        float
            The dipole moment of the ResidueGroup in Debye (D).
        """
        return self._rg.atoms.dipole_moment() * EA2D

    @property
    @_on_whole
    def dipole(self) -> np.ndarray:
        """Calculate the dipole vector of the ResidueGroup.

        Returns
        -------
        np.ndarray
            The dipole vector of the ResidueGroup in Debye (D).
        """
        return self._rg.atoms.dipole_vector() * EA2D

    @property
    @_on_whole
    def shape_parameter(self) -> float:
        """Calculate the shape parameter of the ResidueGroup.

        Returns
        -------
        float
            The shape parameter of the ResidueGroup, a measure of its geometric anisotropy.
        """
        return self._rg.shape_parameter()

    @property
    @_on_whole
    def bsphere(self) -> tuple[float, np.ndarray]:
        """Calculate the bounding sphere of the ResidueGroup.

        Returns
        -------
        tuple[float, np.ndarray,]
            The radius and center of the bounding sphere.
        """
        return self._rg.bsphere()

    @property
    @_on_whole
    def radius_of_gyration(self) -> float:
        """Calculate the radius of gyration of the ResidueGroup.

        Returns
        -------
        float
            The radius of gyration of the ResidueGroup.
        """
        return self._rg.radius_of_gyration()

    @property
    def radius(self) -> float:
        """The radius of the ResidueGroup.

        Returns
        -------
        float
            The radius of the ResidueGroup, equivalent to the radius of gyration.
        """
        return self.radius_of_gyration

    @property
    def diameter(self) -> float:
        """Calculate the diameter of the ResidueGroup.

        Returns
        -------
        float
            The diameter of the ResidueGroup, calculated as twice the radius of gyration.
        """
        return 2 * self.radius_of_gyration

    @property
    def volume(self) -> float:
        """Calculate the ResidueGroup volume assuming a spherical shape and using the radius of gyration.

        Returns
        -------
        float
            Volume in cubic angstroms.
        """
        r = self.radius_of_gyration
        return 4 * np.pi * r**3 / 3  # angstrom^3

    @property
    def density(self) -> float:
        """Calculate the density of the ResidueGroup.

        Returns
        -------
        float
            The density of the ResidueGroup in g/cm^3.
        """
        return (self.mass / self.volume) * 0.602214076  # g/cm^3

    @property
    def charge(self) -> float:
        """Calculate the total charge of the ResidueGroup.

        Returns
        -------
        float
            The total charge of the ResidueGroup.
        """
        return self._rg.total_charge()

    def __len__(self) -> int:
        """Get the number of residues in the ResidueGroup.

        Returns
        -------
        int
            The number of residues in the ResidueGroup.
        """
        return len(self._rg)


class Cluster(MDAResidueGroupAnalyzer):
    """A class representing a molecular cluster.

    This class extends `MDAResidueGroupAnalyzer` to manage molecular clusters as graphs using NetworkX.
    It provides methods for adding/removing molecules and connections, merging clusters, separating
    sub-clusters, and updating clusters from connectivity tables.

    Attributes
    ----------
    initial_time : float
        The time at which the cluster was initialized.
    cluster : nx.Graph
        The graph representation of the molecular cluster.
    _cm : np.ndarray
        The center of mass of the cluster.
    _id : int
        A unique identifier for the cluster.
    """

    __cls_id = 1

    __slots__ = [
        "initial_time",
        "cluster",
        "_cm",
        "_id",
    ]

    def __init__(
        self,
        universe: mda.Universe,
        subconntab: ConnectionTable._SubConnTable | None = None,
        *,
        cluster_id: int | None = None,
    ) -> None:
        """Initialize a Cluster instance.

        Parameters
        ----------
        universe : mda.Universe
            The MDAnalysis Universe object associated with the cluster.
        subconntab : ConnTable._SubConnTable | None
            An optional connectivity table for initializing the cluster graph.
        cluster_id : int | None
            The cluster's id, as handed out by whoever tracks it (e.g. a
            `ClusterTracker`). If None, the next of a counter shared by every
            `Cluster` is used.
        """
        self.initial_time: float = universe.coord.time
        self.cluster: nx.Graph = (
            nx.Graph(subconntab.graph) if subconntab else nx.Graph()
        )
        self._cm: np.ndarray = (
            subconntab._cm if subconntab else np.empty(3)
        )  # TODO: Change to keep the sum of center of masses of molecules

        super().__init__(universe, self.cluster)

        if cluster_id is None:
            cluster_id = Cluster.__cls_id
            Cluster.__cls_id += 1
        self._id = cluster_id

    @classmethod
    def _from_graph(cls, uni: mda.Universe, graph: nx.Graph) -> "Cluster":
        """Create a Cluster instance from a graph.

        Parameters
        ----------
        uni : mda.Universe
            The MDAnalysis Universe object.
        graph : nx.Graph
            The graph representation of the cluster.

        Returns
        -------
        Cluster
            A new Cluster instance.
        """
        tmp_cls = cls(uni, None)
        tmp_cls.cluster = graph
        tmp_cls.__update_rg()
        tmp_cls.__recalculate_cm()
        return tmp_cls

    @property
    def id(self) -> int:
        """The unique identifier of the cluster.

        Returns
        -------
        int
            The unique identifier of the cluster.
        """
        return self._id

    @property
    def cm(self) -> np.ndarray:
        """The center of mass of the cluster.

        Returns
        -------
        np.ndarray
            The center of mass of the cluster.
        """
        return self._cm

    @cm.setter
    def cm(self, value: np.ndarray) -> None:
        """Set the center of mass of the cluster.

        Parameters
        ----------
        value : np.ndarray
            The new center of mass.
        """
        self._cm = value

    def __recalculate_cm(self) -> None:
        """Recalculate the center of mass of the cluster."""
        self._cm = self._rg.center_of_mass()

    def __update_rg(self) -> None:
        """Update the ResidueGroup associated with the cluster."""
        self._rg = core.groups.ResidueGroup(np.array(self.cluster) - 1, self.uni)

    def add_mol(self, ref_mol: int, mol: int, resname: str, dist: float) -> None:
        """Add a molecule to the cluster.

        Parameters
        ----------
        ref_mol : int
            The reference molecule already in the cluster.
        mol : int
            The molecule to add.
        resname : str
            The residue name of the molecule.
        dist : float
            The distance between the reference molecule and the new molecule.

        Raises
        ------
        ValueError
            If ``ref_mol`` not in the cluster, or ``mol`` already in cluster.
        """
        if ref_mol not in self:
            raise ValueError(f"ref_mol {ref_mol} not in the cluster")

        if mol in self:
            raise ValueError(f"mol {mol} already in the cluster, use add_con instead")

        self.cluster.add_node(mol, name=resname)
        self.add_con(ref_mol, mol, dist)
        self._rg += core.groups.ResidueGroup([mol - 1], self.uni)
        self.__recalculate_cm()

    def add_con(self, moli: int, molj: int, dist: float) -> None:
        """Add a connection between two molecules in the cluster.

        Parameters
        ----------
        moli : int
            The first molecule.
        molj : int
            The second molecule.
        dist : float
            The distance between the two molecules.

        Raises
        ------
        ValueError
            If ``moli`` or ``molj`` not in cluster.
        """
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster, use add_mol instead")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster, use add_mol instead")

        # TODO:  use epsilon around 0
        if dist == 0.0:
            raise ValueError("dist is zero. Check your trajectory")

        self.cluster.add_edge(moli, molj, distance=dist, weight=np.exp(1 / dist))

    def remove_mol(self, mol: int) -> None:
        """Remove a molecule from the cluster.

        Parameters
        ----------
        mol : int
            The molecule to remove.

        Raises
        ------
        ValueError
            If ``mol``not in the cluster.
        """
        if mol not in self:
            raise ValueError(f"mol {mol} not in the cluster")

        self.cluster.remove_node(mol)
        self._rg -= core.groups.ResidueGroup([mol - 1], self.uni)
        self.__recalculate_cm()

    def remove_con(self, moli: int, molj: int) -> None:
        """Remove a connection between two molecules in the cluster.

        Parameters
        ----------
        moli : int
            The first molecule.
        molj : int
            The second molecule.

        Raises
        ------
        ValueError
            If ``moli`` or ``molj`` not in cluster.
        """
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster")

        self.cluster.remove_edge(moli, molj)

    def neighbors(self, ref: int, *, level: int | None = None) -> set:
        """Get the neighbors of a molecule in the cluster.

        Parameters
        ----------
        ref : int
            The reference molecule.
        level : int | None
            The level of neighbors to retrieve.

        Returns
        -------
        set
            A set of neighboring molecules.
        """
        if level:
            return set(
                nx.single_target_shortest_path(self._graph, ref, cutoff=level)
            ) - {ref}

        return set(self._graph.neighbors(ref))

    def remove_cons(self, ref_mol: int, cons: Iterable[int]) -> None:
        """Remove multiple connections from a molecule.

        Parameters
        ----------
        ref_mol : int
            The reference molecule.
        cons : Iterable[int]
            The connections to remove.
        """
        for con in cons:
            self.remove_con(ref_mol, con)

    def get_age(self) -> float:
        """Get the age of the cluster.

        Returns
        -------
        float
            The age of the cluster.
        """
        return self.uni.coord.time - self.initial_time

    def merge(self, other: "Cluster") -> None:
        """Merge another cluster into this cluster.

        Parameters
        ----------
        other : Cluster
            The cluster to merge.
        """
        self.cluster = nx.compose(self.cluster, other._graph)
        self.__update_rg()
        self.__recalculate_cm()

    def separate(self) -> list["Cluster"]:
        """Separate the cluster into sub-clusters.

        Returns
        -------
        list[Cluster]
            A list of new Cluster instances representing the sub-clusters.
        """
        sub_clusters = [
            self.cluster.subgraph(c).copy()
            for c in sorted(
                nx.connected_components(self.cluster), key=len, reverse=True
            )
        ]

        self.cluster = sub_clusters[0]
        self.__update_rg()
        self.__recalculate_cm()

        return [Cluster._from_graph(self.uni, nx.Graph(g)) for g in sub_clusters[1:]]

    @property
    def _graph(self) -> nx.Graph:
        """The graph representation of the cluster.

        Returns
        -------
        nx.Graph
            The graph representation of the cluster.
        """
        return self.cluster

    @property
    def size(self) -> int:
        """The size of the cluster.

        Returns
        -------
        int
            The number of molecules in the cluster.
        """
        return len(self)

    def update_from_conntable(self, conn: ConnectionTable._SubConnTable) -> None:
        """Update the cluster from a connectivity table.

        Parameters
        ----------
        conn : ConnTable._SubConnTable
            The connectivity table to update from.
        """
        self.cluster = nx.Graph(conn.graph)
        self._cm = conn.cm
        self.__update_rg()

    def get_dist(self, moli: int, molj: int) -> float:
        """Get the distance between two molecules in the cluster.

        Parameters
        ----------
        moli : int
            The first molecule.
        molj : int
            The second molecule.

        Returns
        -------
        float
            The distance between the two molecules.

        Raises
        ------
        ValueError
            If ``moli`` or ``molj`` not in cluster.
        """
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster")

        return self.cluster[moli][molj]["distance"]

    def set_dist(self, moli: int, molj: int, dist: float) -> None:
        """Set the distance between two molecules in the cluster.

        Parameters
        ----------
        moli : int
            The first molecule.
        molj : int
            The second molecule.
        dist : float
            The new distance between the two molecules.

        Raises
        ------
        ValueError
            If ``moli`` or ``molj`` not in cluster, or if ``dist`` is 0.0
        """
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster")

        # TODO: use eps around 0.0
        if dist == 0.0:
            raise ValueError("dist is zero")

        self.cluster[moli][molj]["distance"] = dist

    def __contains__(self, item: int) -> bool:
        """Check if a molecule is in the cluster.

        Parameters
        ----------
        item : int
            The residue ID of the molecule to check.

        Returns
        -------
        bool
            True if the molecule is in the cluster, False otherwise.
        """
        return item in self.cluster

    def __iter__(self) -> Iterator[int]:
        """Iterate over the residue IDs of the molecules in the cluster.

        Returns
        -------
        Iterator[int]
            An iterator over the residue IDs of the molecules in the cluster.
        """
        return iter(self.cluster.copy())

    def __getitem__(self, key: int) -> dict:
        """Get the attributes of a molecule in the cluster.

        Parameters
        ----------
        key : int
            The molecule to retrieve.

        Returns
        -------
        dict
            The attributes of the molecule.
        """
        return self.cluster[key]

    def __eq__(self, other: object) -> bool:
        """Check if two clusters are equal.

        Parameters
        ----------
        other : object
            The other cluster to compare.

        Returns
        -------
        bool
            True if the clusters are equal, False otherwise.
        """
        if not isinstance(other, Cluster):
            return False

        return nx.utils.graphs_equal(self.cluster, other.cluster)

    # Cluster is mutable (add_mol/remove_mol/merge/... change self.cluster in
    # place), and __eq__ above is value-based. Hashing by value would break if
    # a Cluster's contents changed after being placed in a set/dict; hashing by
    # identity would violate "equal objects must hash equal". So it stays
    # unhashable — made explicit here rather than relying on Python's implicit
    # __hash__ = None whenever __eq__ is defined without __hash__.
    __hash__ = None

    def __len__(self) -> int:
        """Get the number of molecules in the cluster.

        Returns
        -------
        int
            The number of molecules in the cluster.
        """
        return len(self.cluster)

    def __str__(self) -> str:
        """Get a string representation of the cluster.

        Returns
        -------
        str
            A string representation of the cluster.
        """
        return self.cluster.edges.data().__str__()
