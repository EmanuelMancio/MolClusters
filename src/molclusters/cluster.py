# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides the groups of molecules that MolClusters tracks and analyses.

Classes:
--------
    - MolGroup: A group of molecules (residues) of a Universe, with properties such as
      mass, charge, dipole moment, radius of gyration, shape and density, computed on
      the group made whole across periodic boundaries. Nuclei are MolGroups.

    - Cluster: A MolGroup that is a tracked molecular cluster: it has an id, a birth
      time and a read-only graph of the connections between its molecules. Clusters
      are created and updated by a `ClusterTracker`; everyone else only reads them.

Constants:
----------
    - EA2D: Conversion factor for dipole moment from atomic units to Debye.

Dependencies:
-------------
    - MDAnalysis: For molecular dynamics trajectory and structure analysis.
    - NetworkX: For the graph of a cluster's connections.
    - NumPy: For numerical computations.
"""

from collections import Counter
from contextlib import contextmanager
from functools import wraps
from typing import Callable, Iterable, Iterator, Self

import MDAnalysis as mda
import networkx as nx
import numpy as np
from MDAnalysis import core
from MDAnalysis.lib.distances import apply_PBC

from .conntable import ConnectionTable

EA2D = 1 / 0.3934303


def _on_whole[T](method: Callable[..., T]) -> Callable[..., T]:
    """Run `method` with the group's atoms temporarily made whole (see `whole`).

    Returns
    -------
    Callable
        The wrapped method.
    """

    @wraps(method)
    def wrapper(self: "MolGroup", *args: object, **kwargs: object) -> T:
        with self.whole():
            return method(self, *args, **kwargs)

    return wrapper


class MolGroup:
    """A group of molecules (residues) of a Universe, and its properties.

    Geometric properties are computed on the group made whole across periodic
    boundaries (see `whole`), and always for the Universe's current frame.

    Attributes
    ----------
    uni : MDAnalysis.Universe
        The MDAnalysis Universe object associated with the group.
    _rg : MDAnalysis.core.groups.ResidueGroup
        The group's residues.
    __whole_key : tuple[int, bytes] | None
        Frame and residue indices that `__whole_positions` was computed for.
    __whole_positions : np.ndarray | None
        Cached whole positions of the group's atoms (see `whole`).
    __whole_shift : np.ndarray | None
        The translation `__whole_positions` applied to the group, besides whole
        periodic images (see `center_of_mass`).
    """

    __slots__ = ["uni", "_rg", "__whole_key", "__whole_positions", "__whole_shift"]

    def __init__(
        self, universe: mda.Universe, residues: Iterable[int] | core.groups.ResidueGroup
    ) -> None:
        """Initialize the MolGroup with a universe and residues.

        Parameters
        ----------
        universe : mda.Universe
            The MDAnalysis Universe object associated with the group.
        residues : Union[Iterable[int], core.groups.ResidueGroup]
            The residues of the group, either as an iterable with residue IDs (1-based)
            or a ResidueGroup.
        """
        self.uni = universe

        if isinstance(residues, core.groups.ResidueGroup):
            self._rg = residues
        else:
            self._rg = core.groups.ResidueGroup(np.array(residues) - 1, self.uni)

        self.__whole_key: tuple[int, bytes] | None = None
        self.__whole_positions: np.ndarray | None = None
        self.__whole_shift: np.ndarray | None = None

    def __compute_whole_positions(self) -> tuple[np.ndarray, np.ndarray]:
        """Compute the group's atom positions with the group made whole.

        The group is translated so that its first residue sits at the box center,
        then each residue is wrapped back into the box by its center of geometry,
        which gathers the group around the center (assuming it spans less than half
        the box). The Universe's positions are left exactly as they were.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Whole positions, one row per atom of ``self._rg.atoms``, and the
            translation applied to the group (the wrapping only adds whole box
            vectors on top of it).
        """
        atoms = self._rg.atoms
        original = atoms.positions  # a copy
        try:
            boxcenter = np.sum(self.uni.trajectory.ts.triclinic_dimensions, axis=0) / 2
            ref_mol_cm = self._rg[:1].center_of_mass(unwrap=True)
            shift = boxcenter - ref_mol_cm
            atoms.positions += shift
            atoms.unwrap(compound="residues", reference="cog", inplace=True)
            return atoms.positions, shift
        finally:
            atoms.positions = original

    @contextmanager
    def whole(self) -> Iterator[core.groups.AtomGroup]:
        """Temporarily place the group's atoms in whole, centered positions.

        Inside the ``with`` block the group's atoms are whole across periodic
        boundaries, so geometric properties (and writers) see the real shape. On
        exit, even on error, the atoms go back to their original positions.

        Positions are never changed for good because they belong to the shared
        Universe: other groups over overlapping residues (e.g. a nucleus inside a
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
            self.__whole_positions, self.__whole_shift = (
                self.__compute_whole_positions()
            )
            self.__whole_key = key

        atoms = self._rg.atoms
        original = atoms.positions  # a copy
        atoms.positions = self.__whole_positions
        try:
            yield atoms
        finally:
            atoms.positions = original

    def __add__(self, other: core.groups.ResidueGroup | Self) -> "MolGroup":
        """Combine this group with a ResidueGroup or another MolGroup.

        Parameters
        ----------
        other : Union[core.groups.ResidueGroup, MolGroup]
            The ResidueGroup or MolGroup to combine with.

        Returns
        -------
        MolGroup
            A new MolGroup with the residues of both.
        """
        if isinstance(other, core.groups.ResidueGroup):
            return MolGroup(self.uni, self._rg + other)
        elif isinstance(other, MolGroup):
            return MolGroup(self.uni, self._rg + other.residues)

        return NotImplemented

    @property
    def residues(self) -> core.groups.ResidueGroup:
        """The group's residues.

        Returns
        -------
        core.groups.ResidueGroup
            The residues of the group.
        """
        return self._rg

    @property
    def atoms(self) -> core.groups.AtomGroup:
        """The atoms of the group's residues.

        Returns
        -------
        core.groups.AtomGroup
            The atoms of the group, in the Universe's current positions (see
            `whole` for whole positions).
        """
        return self._rg.atoms

    @property
    def size(self) -> int:
        """The number of molecules in the group.

        Returns
        -------
        int
            The number of residues in the group.
        """
        return len(self)

    @property
    def resnames(self) -> list[str]:
        """The residue names of the group.

        Returns
        -------
        list[str]
            A list of residue names in the group.
        """
        return self._rg.resnames

    @property
    def resids(self) -> list[int]:
        """The residue IDs of the group.

        Returns
        -------
        list[int]
            A list of residue IDs in the group.
        """
        return self._rg.resids

    @property
    def composition(self) -> Counter[str]:
        """The number of molecules of each residue name in the group.

        Returns
        -------
        Counter[str]
            Residue name -> number of molecules with that name.
        """
        return Counter(self.resnames)

    @property
    def mass(self) -> float:
        """Calculate the total mass of the group.

        Returns
        -------
        float
            The total mass of the group.
        """
        return self._rg.total_mass()

    @property
    def center_of_mass(self) -> np.ndarray:
        """Calculate the center of mass of the group, made whole.

        Unlike MDAnalysis's own `center_of_mass`, a group split across periodic
        boundaries gets the center of its whole shape, not a point between its
        pieces.

        Returns
        -------
        np.ndarray
            The center of mass, wrapped into the primary unit cell.
        """
        with self.whole() as atoms:
            center = atoms.center_of_mass()
        return apply_PBC(center - self.__whole_shift, self.uni.dimensions)

    @property
    @_on_whole
    def sphericity(self) -> float:
        """Calculate the sphericity of the group.

        Returns
        -------
        float
            The sphericity of the group, a measure of how spherical the shape is.
            A value near zero represents a spherical shape.
        """
        return 1 - self._rg.asphericity()

    @property
    @_on_whole
    def dipole_moment(self) -> float:
        """Calculate the dipole moment of the group.

        Returns
        -------
        float
            The dipole moment of the group in Debye (D).
        """
        return self._rg.atoms.dipole_moment() * EA2D

    @property
    @_on_whole
    def dipole(self) -> np.ndarray:
        """Calculate the dipole vector of the group.

        Returns
        -------
        np.ndarray
            The dipole vector of the group in Debye (D).
        """
        return self._rg.atoms.dipole_vector() * EA2D

    @property
    @_on_whole
    def shape_parameter(self) -> float:
        """Calculate the shape parameter of the group.

        Returns
        -------
        float
            The shape parameter of the group, a measure of its geometric anisotropy.
        """
        return self._rg.shape_parameter()

    @property
    @_on_whole
    def bsphere(self) -> tuple[float, np.ndarray]:
        """Calculate the bounding sphere of the group.

        Returns
        -------
        tuple[float, np.ndarray,]
            The radius and center of the bounding sphere, the center in the
            group's whole positions (see `whole`).
        """
        return self._rg.bsphere()

    @property
    @_on_whole
    def radius_of_gyration(self) -> float:
        """Calculate the radius of gyration of the group.

        Returns
        -------
        float
            The radius of gyration of the group.
        """
        return self._rg.radius_of_gyration()

    @property
    def radius(self) -> float:
        """The radius of the group.

        Returns
        -------
        float
            The radius of the group, equivalent to the radius of gyration.
        """
        return self.radius_of_gyration

    @property
    def diameter(self) -> float:
        """Calculate the diameter of the group.

        Returns
        -------
        float
            The diameter of the group, calculated as twice the radius of gyration.
        """
        return 2 * self.radius_of_gyration

    @property
    def volume(self) -> float:
        """Calculate the group volume assuming a spherical shape and using the radius of gyration.

        Returns
        -------
        float
            Volume in cubic angstroms.
        """
        r = self.radius_of_gyration
        return 4 * np.pi * r**3 / 3  # angstrom^3

    @property
    def density(self) -> float:
        """Calculate the density of the group.

        Returns
        -------
        float
            The density of the group in g/cm^3.
        """
        return (self.mass / self.volume) * 0.602214076  # g/cm^3

    @property
    def charge(self) -> float:
        """Calculate the total charge of the group.

        Returns
        -------
        float
            The total charge of the group.
        """
        return self._rg.total_charge()

    def __len__(self) -> int:
        """Get the number of molecules in the group.

        Returns
        -------
        int
            The number of residues in the group.
        """
        return len(self._rg)

    def __iter__(self) -> Iterator[int]:
        """Iterate over the residue IDs of the group's molecules.

        Returns
        -------
        Iterator[int]
            An iterator over the residue IDs of the group's molecules.
        """
        return iter(self._rg.resids.tolist())

    def __contains__(self, item: int) -> bool:
        """Check if a molecule is in the group.

        Parameters
        ----------
        item : int
            The residue ID of the molecule to check.

        Returns
        -------
        bool
            True if the molecule is in the group, False otherwise.
        """
        return bool(np.isin(item, self._rg.resids))


class Cluster(MolGroup):
    """A tracked molecular cluster.

    A `MolGroup` with an id, a birth time and the graph of the connections between
    its molecules. Clusters are created and updated by a `ClusterTracker` and are
    otherwise read-only: the graph is frozen, and a new one replaces it every frame.

    A Cluster is the same object for as long as it lives, and its properties are
    those of the Universe's current frame, so it describes a given frame only while
    the Universe is at that frame.

    Clusters compare and hash by identity, which is the tracked cluster: the tracker
    keeps one object per cluster for its whole life. To keep data per cluster, key
    it by `id` rather than by the Cluster, so a dead cluster (and its cached whole
    positions) isn't kept in memory.

    Attributes
    ----------
    _graph : nx.Graph
        The frozen graph of the cluster's molecules and their connections.
    _id : int
        The cluster's id.
    _birth_time : float
        The time of the frame the cluster was created at.
    """

    __slots__ = ["_graph", "_id", "_birth_time"]

    def __init__(
        self,
        universe: mda.Universe,
        subconntab: ConnectionTable._SubConnTable,
        *,
        cluster_id: int,
    ) -> None:
        """Initialize a Cluster from a connected group of the current frame.

        Parameters
        ----------
        universe : mda.Universe
            The MDAnalysis Universe object associated with the cluster.
        subconntab : ConnTable._SubConnTable
            The connected group of molecules the cluster is made of.
        cluster_id : int
            The cluster's id, as handed out by whoever tracks it (e.g. a
            `ClusterTracker`).
        """
        self._birth_time: float = universe.coord.time
        self._graph: nx.Graph = nx.freeze(nx.Graph(subconntab.graph))
        self._id = cluster_id

        super().__init__(universe, self._graph)

    def _update(self, subconntab: ConnectionTable._SubConnTable) -> None:
        """Continue the cluster as a connected group of the current frame.

        Only the cluster's tracker calls this.

        Parameters
        ----------
        subconntab : ConnTable._SubConnTable
            The connected group of molecules the cluster is now made of.
        """
        self._graph = nx.freeze(nx.Graph(subconntab.graph))
        self._rg = core.groups.ResidueGroup(np.array(self._graph) - 1, self.uni)

    @property
    def id(self) -> int:
        """The cluster's id, unique among the clusters of its tracker.

        Returns
        -------
        int
            The id of the cluster.
        """
        return self._id

    @property
    def graph(self) -> nx.Graph:
        """The molecules of the cluster and their connections, as a frozen graph.

        Nodes are residue IDs; each edge has a ``distance`` and a ``weight``.

        Returns
        -------
        nx.Graph
            The graph of the cluster, which can't be modified.
        """
        return self._graph

    @property
    def birth_time(self) -> float:
        """The time of the frame the cluster was created at.

        Returns
        -------
        float
            The birth time of the cluster.
        """
        return self._birth_time

    @property
    def age(self) -> float:
        """The time since the cluster was created.

        Returns
        -------
        float
            The age of the cluster at the Universe's current frame.
        """
        return self.uni.coord.time - self._birth_time

    def neighbors(self, ref: int, *, level: int | None = None) -> set[int]:
        """Get the neighbors of a molecule in the cluster.

        Parameters
        ----------
        ref : int
            The reference molecule.
        level : int | None
            How many connections away to look (1 when None).

        Returns
        -------
        set[int]
            The molecules at most `level` connections away from `ref`.
        """
        if level:
            return set(
                nx.single_target_shortest_path(self._graph, ref, cutoff=level)
            ) - {ref}

        return set(self._graph.neighbors(ref))

    def distance(self, moli: int, molj: int) -> float:
        """Get the distance between two connected molecules of the cluster.

        Parameters
        ----------
        moli : int
            The first molecule.
        molj : int
            The second molecule.

        Returns
        -------
        float
            The distance of the connection between the two molecules.

        Raises
        ------
        ValueError
            If ``moli`` or ``molj`` is not in the cluster, or they aren't connected.
        """
        for mol in (moli, molj):
            if mol not in self:
                raise ValueError(f"mol {mol} not in the cluster")
        if not self._graph.has_edge(moli, molj):
            raise ValueError(f"mols {moli} and {molj} are not connected")

        return self._graph[moli][molj]["distance"]

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
        return item in self._graph

    def __iter__(self) -> Iterator[int]:
        """Iterate over the residue IDs of the molecules in the cluster.

        Returns
        -------
        Iterator[int]
            An iterator over the residue IDs of the molecules in the cluster.
        """
        return iter(self._graph)

    def __len__(self) -> int:
        """Get the number of molecules in the cluster.

        Returns
        -------
        int
            The number of molecules in the cluster.
        """
        return len(self._graph)

    def __str__(self) -> str:
        """Get a string representation of the cluster.

        Returns
        -------
        str
            The cluster's connections, with their attributes.
        """
        return self._graph.edges.data().__str__()
