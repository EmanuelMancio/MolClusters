# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

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
    - EA2D: Conversion factor for dipole moments from e·Å (MDAnalysis's unit) to
      Debye.

Dependencies:
-------------
    - MDAnalysis: For molecular dynamics trajectory and structure analysis.
    - NetworkX: For the graph of a cluster's connections.
    - NumPy: For numerical computations.
"""

import weakref
from collections import Counter
from contextlib import contextmanager
from functools import wraps
from types import MappingProxyType
from typing import Callable, Iterable, Iterator, Self

import MDAnalysis as mda
import networkx as nx
import numpy as np
from loguru import logger
from MDAnalysis import core
from MDAnalysis.exceptions import NoDataError
from MDAnalysis.guesser.tables import vdwradii
from MDAnalysis.lib.distances import apply_PBC, distance_array, minimize_vectors

from .conntable import ConnectionTable, _whole_residue_offsets

# 1 D = 0.2081943 e·Å (e·Å, not the atomic unit e·a0: 1 D = 0.3934303 e·a0)
EA2D = 1 / 0.2081943

# Universe -> van der Waals radius of each of its atoms (see `_vdw_radii`)
_VDW_RADII: "weakref.WeakKeyDictionary[mda.Universe, np.ndarray]" = (
    weakref.WeakKeyDictionary()
)

# Universe -> topology attribute -> whether the topology has it (see `_has`)
_TOPOLOGY_HAS: "weakref.WeakKeyDictionary[mda.Universe, dict[str, bool]]" = (
    weakref.WeakKeyDictionary()
)

# what a topology without the attribute gets instead, logged once (see `_has`)
_WITHOUT = {
    # groups are made whole along bonds when there are some (see `MolGroup.whole`)
    "bonds": (
        "INFO",
        "The topology has no bonds: molecules are made whole by minimum image "
        "around their first atom, so each must span less than half the box.",
    ),
    "charges": (
        "WARNING",
        "The topology has no partial charges: charges and dipole moments are NaN "
        "(null in the report).",
    ),
}


def _has(universe: mda.Universe, attr: str) -> bool:
    """Check whether a Universe's topology has an attribute, noting once if not.

    Coordinate-only topologies (``.gro``, ``.pdb``) have no bonds or charges.

    Parameters
    ----------
    universe : mda.Universe
        The Universe to check.
    attr : str
        A topology attribute listed in `_WITHOUT`.

    Returns
    -------
    bool
        Whether the topology has it (bonds possibly none at all, e.g. only ions).
    """
    known = _TOPOLOGY_HAS.setdefault(universe, {})
    if attr not in known:
        # as MDAnalysis' own unwrap checks bonds, on one atom to keep it cheap
        known[attr] = hasattr(universe.atoms[:1], attr)
        if not known[attr]:
            logger.log(*_WITHOUT[attr])
    return known[attr]


def _vdw_radii(universe: mda.Universe) -> np.ndarray:
    """Look up the van der Waals radius of every atom of a Universe, by element.

    Computed once per Universe, since the groups read it every frame.

    Parameters
    ----------
    universe : mda.Universe
        The Universe whose atoms to look up.

    Returns
    -------
    np.ndarray
        One radius per atom, in angstroms, NaN for an element without one.

    Raises
    ------
    ValueError
        If the atoms have no elements.
    """
    if universe not in _VDW_RADII:
        try:
            elements = universe.atoms.elements
        except NoDataError as err:
            raise ValueError(
                "The equivalent sphere's radius needs the atoms' elements, which the "
                "topology doesn't have; guess them with "
                "`universe.guess_TopologyAttrs(to_guess=['elements'])`."
            ) from err

        names, which = np.unique(elements, return_inverse=True)
        radii = np.array([vdwradii.get(name.upper(), np.nan) for name in names])
        _VDW_RADII[universe] = radii[which]

    return _VDW_RADII[universe]


def _readonly[T](value: T) -> T:
    """Make an array, or the arrays of a tuple, read-only, in place.

    A group's arrays describe it, so they can't be changed in place: a caller
    wanting to change one works on a copy (``array.copy()``).

    Parameters
    ----------
    value : T
        An array, a tuple possibly holding arrays, or anything else (left as is).

    Returns
    -------
    T
        `value` itself.
    """
    if isinstance(value, np.ndarray):
        value.flags.writeable = False
    elif isinstance(value, tuple):
        for item in value:
            _readonly(item)
    return value


def _own_residues(ix: np.ndarray, universe: mda.Universe) -> core.groups.ResidueGroup:
    """Build a ResidueGroup over its own copy of the residues' indices.

    An MDAnalysis group's `ix` is its own index array (which can't be made
    read-only: MDAnalysis needs it writable), so a group sharing it with another
    would have its molecules swapped by ``other.ix[0] = ...``.

    Parameters
    ----------
    ix : np.ndarray
        The residues' indices in the topology (copied).
    universe : mda.Universe
        The Universe the residues belong to.

    Returns
    -------
    core.groups.ResidueGroup
        The residues.
    """
    return core.groups.ResidueGroup(np.array(ix, dtype=np.intp), universe)


def _frozen_graph(graph: nx.Graph) -> nx.Graph:
    """Copy a graph into one that can't be modified, attributes included.

    `nx.freeze` only stops nodes and edges from being added or removed, and leaves
    their attribute dicts (and the graph's own) writable, so those are swapped for
    read-only views. Reading, subgraph views and copies (``nx.Graph(graph)``,
    ``graph.copy()``, which are modifiable) work as usual.

    Parameters
    ----------
    graph : nx.Graph
        The graph to copy.

    Returns
    -------
    nx.Graph
        The frozen copy.
    """
    frozen = nx.Graph(graph)
    frozen.graph = MappingProxyType(frozen.graph)
    for node, attrs in list(frozen._node.items()):
        frozen._node[node] = MappingProxyType(attrs)
    for u, v, attrs in list(frozen.edges(data=True)):
        # one dict per edge, shared by both directions
        frozen._adj[u][v] = frozen._adj[v][u] = MappingProxyType(attrs)
    return nx.freeze(frozen)


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


def _per_frame[T](method: Callable[["MolGroup"], T]) -> Callable[["MolGroup"], T]:
    """Compute `method` once per frame and residue set (see `MolGroup._frame_cache`).

    Every caller gets the same value, so arrays (also inside a tuple) are made
    read-only (see `_readonly`): one caller can't change what the next one reads.

    Returns
    -------
    Callable
        The wrapped method.
    """
    name = method.__name__

    @wraps(method)
    def wrapper(self: "MolGroup") -> T:
        cache = self._frame_cache()
        if name not in cache:
            cache[name] = _readonly(method(self))
        return cache[name]

    return wrapper


class MolGroup:
    """A group of molecules (residues) of a Universe, and its properties.

    Geometric properties are computed on the group made whole across periodic
    boundaries (see `whole`), and always for the Universe's current frame.

    The group describes its residues, so it can't be changed through what it hands
    out: its arrays are read-only (to change one, work on a copy:
    ``array.copy()``), and `residues` is a new group every time. Positions and
    topology attributes belong to the shared Universe, not the group, and changing
    them changes every group.

    Attributes
    ----------
    _uni : MDAnalysis.Universe
        The MDAnalysis Universe object associated with the group.
    _rg : MDAnalysis.core.groups.ResidueGroup
        The group's residues, over indices no one else holds.
    __cache_key : tuple[int, bytes] | None
        Frame and residue indices that `__cache` holds values for.
    __cache : dict[str, object]
        Values computed for the current frame and residues (see `_frame_cache`).
    """

    __slots__ = ["_uni", "_rg", "__cache_key", "__cache"]

    def __init__(
        self, universe: mda.Universe, residues: Iterable[int] | core.groups.ResidueGroup
    ) -> None:
        """Initialize the MolGroup with a universe and residues.

        Parameters
        ----------
        universe : mda.Universe
            The MDAnalysis Universe object associated with the group.
        residues : Union[Iterable[int], core.groups.ResidueGroup]
            The residues of the group, either as an iterable with residue IDs or a
            ResidueGroup. Residue IDs are taken as positions in the topology plus
            one, so they must be numbered 1 to N (see `conntable.check_resids`).
        """
        self._uni = universe

        if isinstance(residues, core.groups.ResidueGroup):
            ix = residues.ix
        else:
            ix = np.array(residues) - 1
        self._rg = _own_residues(ix, universe)

        self.__cache_key: tuple[int, bytes] | None = None
        self.__cache: dict[str, object] = {}

    def _frame_cache(self) -> dict[str, object]:
        """Get the values computed for the current frame and residue set.

        Positions are taken to change only from one frame to the next, so values
        computed from them (whole positions, and `_per_frame` properties) are
        kept until the frame, or the group's residues, change. A frame whose
        positions are modified in place keeps its earlier values.

        Returns
        -------
        dict[str, object]
            The cache, emptied if the frame or residues changed since last call.
        """
        key = (self._uni.trajectory.ts.frame, self._rg.ix.tobytes())
        if self.__cache_key != key:
            self.__cache = {}
            self.__cache_key = key
        return self.__cache

    def __whole(self) -> tuple[np.ndarray, np.ndarray]:
        """Get the group's whole positions and shift, computed once per frame.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            See `__compute_whole_positions`.
        """
        cache = self._frame_cache()
        if "whole" not in cache:
            cache["whole"] = self.__compute_whole_positions()
        return cache["whole"]

    def __compute_whole_positions(self) -> tuple[np.ndarray, np.ndarray]:
        """Compute the group's atom positions with the group made whole.

        The group is translated so that its first residue sits at the box center,
        and each residue is made whole and wrapped into the box by its center of
        geometry, which gathers a group spanning less than half the box around the
        center. A larger group would be torn apart that way, so the residues are
        then placed one after another along a spanning tree of the group (see
        `_placement_order`), each at the periodic image nearest its tree neighbour.
        The Universe's positions are left exactly as they were.

        Residues are made whole along their bonds, or, in a topology without
        bonds, by minimum image around their first atom (see `_has`).

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            Whole positions, one row per atom of ``self._rg.atoms``, and the
            translation applied to the group (the wrapping only adds whole box
            vectors on top of it).
        """
        atoms = self._rg.atoms
        boxcenter = np.sum(self._uni.trajectory.ts.triclinic_dimensions, axis=0) / 2
        if _has(self._uni, "bonds"):
            positions, shift = self.__unwrap_along_bonds(atoms, boxcenter)
        else:
            positions, shift = self.__unwrap_by_minimum_image(atoms, boxcenter)

        if len(self._rg) > 1:
            self.__place_residues(atoms, positions)
        return positions, shift

    def __unwrap_along_bonds(
        self, atoms: core.groups.AtomGroup, boxcenter: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Center the group's first residue and make each residue whole by its bonds.

        Parameters
        ----------
        atoms : core.groups.AtomGroup
            The group's atoms.
        boxcenter : np.ndarray
            The center of the box.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            See `__compute_whole_positions`, before the residues are placed.
        """
        original = atoms.positions  # a copy
        try:
            ref_mol_cm = self._rg[:1].center_of_mass(unwrap=True)
            shift = boxcenter - ref_mol_cm
            atoms.positions += shift
            atoms.unwrap(compound="residues", reference="cog", inplace=True)
            return atoms.positions, shift
        finally:
            atoms.positions = original

    def __unwrap_by_minimum_image(
        self, atoms: core.groups.AtomGroup, boxcenter: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Center the group's first residue and make each residue whole, bond-free.

        Each residue is made whole by minimum image around its first atom (see
        `conntable._whole_residue_offsets`), which gives what `__unwrap_along_bonds`
        does, up to rounding, for residues spanning less than half the box.

        Parameters
        ----------
        atoms : core.groups.AtomGroup
            The group's atoms.
        boxcenter : np.ndarray
            The center of the box.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            See `__compute_whole_positions`, before the residues are placed.
        """
        box = self._uni.dimensions
        anchors, offsets, residue = _whole_residue_offsets(atoms, box)
        positions = anchors[residue] + offsets

        first = atoms.resindices == self._rg.ix[0]
        masses = atoms.masses[first]
        ref_mol_cm = masses @ positions[first] / masses.sum()
        shift = boxcenter - ref_mol_cm
        positions += shift.astype(positions.dtype)

        # each residue's center of geometry into the box, as unwrap's reference="cog"
        counts = np.bincount(residue)
        centers = (
            np.stack([np.bincount(residue, positions[:, k]) for k in range(3)], axis=1)
            / counts[:, None]
        )
        wrapped = apply_PBC(centers.astype(np.float32), box)
        positions += (wrapped - centers)[residue].astype(positions.dtype)
        return positions, shift

    def __place_residues(
        self, atoms: core.groups.AtomGroup, positions: np.ndarray
    ) -> None:
        """Move each residue to the periodic image nearest its placed neighbour.

        Changes `positions` in place, by whole box vectors per residue: none at all
        for a group spanning less than half the box.

        Parameters
        ----------
        atoms : core.groups.AtomGroup
            The group's atoms.
        positions : np.ndarray
            Their positions, each residue whole and wrapped into the box.
        """
        box = self._uni.dimensions
        if box is None or not np.any(box[:3]):
            return
        # each atom's residue, as its row in self._rg
        by_ix = np.argsort(self._rg.ix)
        rows = by_ix[np.searchsorted(self._rg.ix[by_ix], atoms.resindices)]
        counts = np.bincount(rows, minlength=len(self._rg))
        centers = (
            np.stack(
                [np.bincount(rows, positions[:, k], len(self._rg)) for k in range(3)],
                axis=1,
            )
            / counts[:, None]
        )

        # in a rectangular box, residues spanning less than half of it along every
        # axis are all each other's nearest images already: nothing to place, and
        # nothing wrapping around the box (a connection torn by the wrapping into
        # the box would stretch the span past half of it)
        # (plain comparisons: np.allclose and np.ptp cost more than the rest here)
        rectangular = (np.abs(box[3:] - 90.0) < 1e-3).all()
        span = centers.max(axis=0) - centers.min(axis=0)
        if rectangular and (span < box[:3] / 2).all():
            return

        order, parent = self._placement_order(centers, box)
        steps = minimize_vectors(
            (centers[order[1:]] - centers[parent[order[1:]]]).astype(np.float32), box
        )
        placed = centers.copy()
        for row, step in zip(order[1:], steps, strict=True):
            placed[row] = placed[parent[row]] + step

        self._check_placement(placed, centers, box)
        positions += (placed - centers)[rows].astype(positions.dtype)

    def _placement_order(
        self, centers: np.ndarray, box: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Choose the order and tree in which `whole` places the residues.

        A minimum spanning tree of the residues' centers under periodic boundaries,
        grown from the first residue (Prim's algorithm), so each residue is placed
        next to its nearest already-placed one.

        Parameters
        ----------
        centers : np.ndarray
            Each residue's center of geometry, one row per residue of the group.
        box : np.ndarray
            The box dimensions.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            The rows in placement order (the first residue first), and each row's
            parent row in the tree (-1 for the first residue).
        """
        n_res = len(centers)
        distances = distance_array(centers, centers, box=box)
        parent = np.full(n_res, -1)
        placed = np.zeros(n_res, dtype=bool)
        nearest = distances[0].copy()  # to the nearest placed residue
        nearest_to = np.zeros(n_res, dtype=int)
        order = [0]
        placed[0] = True
        for _ in range(n_res - 1):
            row = int(np.argmin(np.where(placed, np.inf, nearest)))
            placed[row] = True
            parent[row] = nearest_to[row]
            order.append(row)
            closer = distances[row] < nearest
            nearest[closer] = distances[row][closer]
            nearest_to[closer] = row
        return np.array(order), parent

    def _check_placement(
        self, placed: np.ndarray, centers: np.ndarray, box: np.ndarray
    ) -> None:
        """Check the placed residues, a hook for subclasses (see `Cluster`).

        Parameters
        ----------
        placed : np.ndarray
            Each residue's center where `whole` placed it.
        centers : np.ndarray
            Each residue's center before.
        box : np.ndarray
            The box dimensions.
        """

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

        The whole positions are cached per frame and residue set (see
        `_frame_cache`), so the unwrap runs once per frame no matter how many
        properties are read.

        Yields
        ------
        core.groups.AtomGroup
            The group's atoms, in whole positions.
        """
        positions, _ = self.__whole()
        atoms = self._rg.atoms
        original = atoms.positions  # a copy
        atoms.positions = positions
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
            return MolGroup(self._uni, self._rg + other)
        elif isinstance(other, MolGroup):
            return MolGroup(self._uni, self._rg + other.residues)

        return NotImplemented

    @property
    def universe(self) -> mda.Universe:
        """The Universe the group belongs to.

        Returns
        -------
        mda.Universe
            The group's Universe.
        """
        return self._uni

    @property
    def residues(self) -> core.groups.ResidueGroup:
        """The group's residues.

        Returns
        -------
        core.groups.ResidueGroup
            The residues of the group, as a new group: changing its indices
            (``ix``) in place leaves the group's own residues as they are.
        """
        return _own_residues(self._rg.ix, self._uni)

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
    def resnames(self) -> np.ndarray:
        """The residue names of the group.

        Returns
        -------
        np.ndarray
            The residue names in the group (read-only).
        """
        return _readonly(self._rg.resnames)

    @property
    def resids(self) -> np.ndarray:
        """The residue IDs of the group.

        Returns
        -------
        np.ndarray
            The residue IDs in the group (read-only).
        """
        return _readonly(self._rg.resids)

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
    @_per_frame
    def center_of_mass(self) -> np.ndarray:
        """Calculate the center of mass of the group, made whole.

        Unlike MDAnalysis's own `center_of_mass`, a group split across periodic
        boundaries gets the center of its whole shape, not a point between its
        pieces.

        Returns
        -------
        np.ndarray
            The center of mass, wrapped into the primary unit cell (read-only).
        """
        with self.whole() as atoms:
            center = atoms.center_of_mass()
        _, shift = self.__whole()
        return apply_PBC(center - shift, self._uni.dimensions)

    @property
    @_per_frame
    @_on_whole
    def sphericity(self) -> float:
        """Calculate how spherical the group is, from its gyration tensor.

        This is one minus the asphericity (Dima & Thirumalai), not Wadell's
        surface-area sphericity: 1 when the three principal moments of the gyration
        tensor are equal (a sphere, but also e.g. a cube), 0 for a rod (all mass on
        a line), and 0.75 for a flat disk or ring.

        Returns
        -------
        float
            1 - asphericity, between 0 (rod) and 1 (spherical).
        """
        return 1 - self._rg.asphericity()

    @property
    @_per_frame
    @_on_whole
    def dipole_moment(self) -> float:
        """Calculate the dipole moment of the group, about its center of mass.

        For a neutral group the dipole is the same about any point. For a group
        with a net charge q (ions, charged nuclei or clusters) it isn't: moving the
        reference point by d changes the dipole by -q d, so its value is the dipole
        about the center of mass specifically, and only comparable between groups
        with that convention in mind.

        Returns
        -------
        float
            The dipole moment of the group in Debye (D), NaN if the topology has no
            partial charges.
        """
        if not _has(self._uni, "charges"):
            return np.nan
        return self._rg.atoms.dipole_moment() * EA2D

    @property
    @_per_frame
    @_on_whole
    def dipole(self) -> np.ndarray:
        """Calculate the dipole vector of the group, about its center of mass.

        See `dipole_moment` for what that reference point means for charged groups.

        Returns
        -------
        np.ndarray
            The dipole vector of the group in Debye (D) (read-only), NaN if the
            topology has no partial charges.
        """
        if not _has(self._uni, "charges"):
            return np.full(3, np.nan)
        return self._rg.atoms.dipole_vector() * EA2D

    @property
    @_per_frame
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
    @_per_frame
    @_on_whole
    def bsphere(self) -> tuple[float, np.ndarray]:
        """Calculate the bounding sphere of the group.

        Returns
        -------
        tuple[float, np.ndarray,]
            The radius and center of the bounding sphere, the center (read-only)
            in the group's whole positions (see `whole`).
        """
        return self._rg.bsphere()

    @property
    @_per_frame
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
    def radius_buffer(self) -> float:
        """Half the mean van der Waals radius of the group's atoms, by element.

        The radius of gyration only sees the atoms' centers, so a sphere built from
        it ends at the outer atoms' centers and leaves their size out. Adding half,
        not all, of the atoms' van der Waals radius to it roughly accounts for
        that, since bonded atoms' van der Waals spheres overlap: on methanol-malic
        acid mixtures, the full radius puts a lone malic acid at 0.70 g/cm^3 and
        half at 1.21 (1.61 for the solid), and barely changes large clusters.

        Returns
        -------
        float
            The buffer added to the equivalent sphere's radius, in angstroms.

        Raises
        ------
        ValueError
            If the atoms have no elements, or an element has no van der Waals
            radius in MDAnalysis' table.
        """
        atoms = self._rg.atoms
        radii = _vdw_radii(self._uni)[atoms.ix]
        if np.isnan(radii).any():
            unknown = sorted(set(atoms.elements[np.isnan(radii)]))
            raise ValueError(
                f"No van der Waals radius is known for element(s) {unknown}."
            )
        return 0.5 * radii.mean()

    @property
    @_per_frame
    def radius(self) -> float:
        """The radius of the group's equivalent sphere.

        That is the uniform sphere with the group's radius of gyration, which has a
        radius of sqrt(5/3) times it (a uniform sphere of radius R has a radius of
        gyration of sqrt(3/5) R), grown by `radius_buffer` to roughly account for
        the size of the atoms.

        Returns
        -------
        float
            The group's equivalent sphere radius, in angstroms.
        """
        return np.sqrt(5 / 3) * self.radius_of_gyration + self.radius_buffer

    @property
    def diameter(self) -> float:
        """Calculate the diameter of the group's equivalent sphere (see `radius`).

        Returns
        -------
        float
            Twice the group's radius, in angstroms.
        """
        return 2 * self.radius

    @property
    def volume(self) -> float:
        """Calculate the volume of the group's equivalent sphere (see `radius`).

        Returns
        -------
        float
            Volume in cubic angstroms.
        """
        return 4 * np.pi * self.radius**3 / 3  # angstrom^3

    @property
    def density(self) -> float:
        """Calculate the density of the group.

        Returns
        -------
        float
            The density of the group in g/cm^3.
        """
        # amu/A^3 -> g/cm^3: 1 amu = 1.66053906660e-24 g, 1 A^3 = 1e-24 cm^3
        return (self.mass / self.volume) * 1.66053906660  # g/cm^3

    @property
    def charge(self) -> float:
        """Calculate the total charge of the group.

        Returns
        -------
        float
            The total charge of the group, NaN if the topology has no partial
            charges.
        """
        if not _has(self._uni, "charges"):
            return np.nan
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
    otherwise read-only: the graph is frozen, attributes included, and a new one
    replaces it every frame; its arrays are read-only (see `MolGroup`).

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
    _wrap_reported : bool
        Whether the cluster was already reported wrapping around the box.
    """

    __slots__ = ["_graph", "_id", "_birth_time", "_wrap_reported"]

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
        self._graph: nx.Graph = _frozen_graph(subconntab.graph)
        self._id = cluster_id
        self._wrap_reported = False

        super().__init__(universe, self._graph)

    def _placement_order(
        self, centers: np.ndarray, box: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Place the residues along the cluster's connections (see `MolGroup`).

        A breadth-first walk of the graph from its first molecule: every step is
        a connection, much shorter than half the box.

        Parameters
        ----------
        centers : np.ndarray
            Each residue's center of geometry, in graph order (as ``self._rg``).
        box : np.ndarray
            The box dimensions.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            The rows in placement order, and each row's parent row (-1 for the
            first molecule).
        """
        nodes = list(self._graph)
        row = {node: i for i, node in enumerate(nodes)}
        parent = np.full(len(nodes), -1)
        order = [0]
        for predecessor, node in nx.bfs_edges(self._graph, nodes[0]):
            order.append(row[node])
            parent[row[node]] = row[predecessor]
        return np.array(order), parent

    def _check_placement(
        self, placed: np.ndarray, centers: np.ndarray, box: np.ndarray
    ) -> None:
        """Report a cluster connected to its own periodic image, once.

        Placing the molecules along a spanning tree leaves out some connections;
        each of those must then join its two molecules by the shortest vector
        between them too. If one spans a box vector instead, the cluster wraps
        around the box, an endless structure that can't be made whole.

        Parameters
        ----------
        placed : np.ndarray
            Each residue's center where `whole` placed it, in graph order.
        centers : np.ndarray
            Each residue's center before.
        box : np.ndarray
            The box dimensions.
        """
        if self._wrap_reported or self._graph.number_of_edges() == 0:
            return

        row = {node: i for i, node in enumerate(self._graph)}
        ends = np.array([(row[u], row[v]) for u, v in self._graph.edges])
        shortest = minimize_vectors(
            (centers[ends[:, 1]] - centers[ends[:, 0]]).astype(np.float32), box
        )
        placed_apart = placed[ends[:, 1]] - placed[ends[:, 0]]
        if np.any(np.linalg.norm(placed_apart - shortest, axis=1) > 1e-2):
            self._wrap_reported = True
            logger.warning(
                f"cluster {self._id} wraps around the periodic box: it is connected "
                "to its own periodic image, so it has no whole shape, and its radius, "
                "volume, density, shape and dipole mean little while it does. "
                "(Reported once per cluster.)"
            )

    def _update(self, subconntab: ConnectionTable._SubConnTable) -> None:
        """Continue the cluster as a connected group of the current frame.

        Only the cluster's tracker calls this.

        Parameters
        ----------
        subconntab : ConnTable._SubConnTable
            The connected group of molecules the cluster is now made of.
        """
        self._graph = _frozen_graph(subconntab.graph)
        self._rg = _own_residues(np.array(self._graph) - 1, self._uni)

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

        Nodes are residue IDs, and each edge joins two connected molecules. A
        "cm" rule's edge has the ``distance`` between their centers of mass; an
        "hb" rule's has the donor-acceptor ``distance`` and D-H-A ``angle`` of the
        shortest of their H-bonds, and ``n_hbonds``, how many they share.

        Returns
        -------
        nx.Graph
            The graph of the cluster, which can't be modified, nor can the
            attributes of its edges, nodes or itself (``nx.Graph(graph)`` gives a
            modifiable copy).
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
        return self._uni.coord.time - self._birth_time

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
