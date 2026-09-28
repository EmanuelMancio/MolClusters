# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

from collections import Counter
from contextlib import AbstractContextManager
from typing import Callable

import networkx as nx
import numpy as np
import pytest
from MDAnalysis import Universe
from MDAnalysis.lib.distances import apply_PBC

from molclusters import cluster
from molclusters.cluster import Cluster, MolGroup
from molclusters.config import MolClsConfig
from molclusters.conntable import ConnectionTable

from .conftest import BOND_STEP, CUTOFF, UniverseFactory


def connected_groups(uni: Universe, cutoff: float = CUTOFF) -> list:
    """The connected groups of the current frame under a ``cm cutoff`` MOL rule.

    Returns
    -------
    list
        The frame's connected groups, largest first.
    """
    config = MolClsConfig(rules={"MOL": {"MOL": f"cm {cutoff}"}})
    conntab = ConnectionTable(uni, config._rules, {"MOL": uni.atoms})
    return list(conntab.subconntables())


@pytest.fixture
def populated_cluster() -> Cluster:
    """Build a real two-molecule Cluster via the same path production code uses.

    Returns
    -------
    Cluster
        A cluster containing the two connected MOL residues from the fixture.
    """
    uni = Universe("tests/data/met-mal/met-mal.tpr", "tests/data/met-mal/start.pdb")
    config = MolClsConfig(rules={"MOL": {"MOL": "cm 15.0"}})
    sels = {res: uni.select_atoms(f"resname {res}") for res in config._rules.all_keys()}
    conntab = ConnectionTable(uni, config._rules, sels)

    return Cluster(uni, next(conntab.subconntables()), cluster_id=1)


@pytest.fixture
def chain(make_universe: UniverseFactory) -> Cluster:
    """A Cluster over the chain 1-2-3, with residues 4 and 5 left outside it.

    Returns
    -------
    Cluster
        Cluster built from the synthetic system's only connected component.
    """
    uni = make_universe([[[1, 2, 3]]], 5)
    (group,) = connected_groups(uni)

    return Cluster(uni, group, cluster_id=1)


def members(group: MolGroup) -> list[int]:
    return sorted(int(r) for r in group.resids)


class TestCluster:
    def test_id_is_the_given_one(self, chain: Cluster):
        assert chain.id == 1

    def test_id_is_required(self, make_universe: UniverseFactory):
        uni = make_universe([[[1, 2]]], 2)
        (group,) = connected_groups(uni)

        with pytest.raises(TypeError, match="cluster_id"):
            Cluster(uni, group)  # type: ignore[call-arg]

    def test_add_unsupported_operand(self, chain: Cluster):
        with pytest.raises(TypeError):
            chain + "not a residue group"

    def test_adding_to_a_cluster_gives_a_plain_group(self, chain: Cluster):
        combined = chain + chain.universe.residues[[3]]

        assert type(combined) is MolGroup
        assert members(combined) == [1, 2, 3, 4]


class TestPopulatedCluster:
    def test_distance_between_mols(self, populated_cluster: Cluster):
        moli = next(iter(populated_cluster))
        molj = next(iter(populated_cluster.neighbors(moli)))

        expected = populated_cluster.graph[moli][molj]["distance"]
        assert populated_cluster.distance(moli, molj) == expected

    def test_contains(self, populated_cluster: Cluster):
        mol = next(iter(populated_cluster))

        assert mol in populated_cluster
        assert -1 not in populated_cluster


class TestMolGroup:
    def test_built_from_resids(self, make_universe: UniverseFactory):
        uni = make_universe([[]], 3)

        group = MolGroup(uni, [1, 3])

        assert members(group) == [1, 3]
        assert group.size == len(group) == 2
        assert list(group.resnames) == ["MOL", "MOL"]
        assert list(group) == [1, 3]
        assert 3 in group
        assert 2 not in group

    def test_residues_and_atoms(self, make_universe: UniverseFactory):
        uni = make_universe([[]], 3)

        group = MolGroup(uni, [1, 3])

        assert group.residues == uni.residues[[0, 2]]
        assert group.atoms == uni.residues[[0, 2]].atoms

    def test_composition_counts_resnames(self, make_universe: UniverseFactory):
        uni = make_universe([[]], 3)

        assert MolGroup(uni, [1, 2, 3]).composition == Counter({"MOL": 3})

    def test_add_residue_group_and_group(self, make_universe: UniverseFactory):
        uni = make_universe([[]], 3)
        first = MolGroup(uni, [1])

        with_residues = first + uni.residues[[1]]
        with_group = first + MolGroup(uni, [3])

        assert members(with_residues) == [1, 2]
        assert members(with_group) == [1, 3]

    def test_physical_properties_are_consistent(self, chain: Cluster):
        assert chain.mass == pytest.approx(3 * (12.011 + 15.999))
        assert chain.charge == pytest.approx(0.0)
        assert chain.radius_of_gyration > 0
        assert chain.radius == pytest.approx(
            np.sqrt(5 / 3) * chain.radius_of_gyration + chain.radius_buffer
        )
        assert chain.diameter == pytest.approx(2 * chain.radius)
        assert chain.volume == pytest.approx(4 / 3 * np.pi * chain.radius**3)
        assert chain.density == pytest.approx(chain.mass / chain.volume * 1.66053906660)

    def test_radius_buffer_averages_over_the_atoms(self, chain: Cluster):
        # C-O residues: van der Waals radii of 1.7 and 1.52 A
        assert chain.radius_buffer == pytest.approx((1.7 + 1.52) / 4)

    def test_shape_properties(self, chain: Cluster):
        assert 0.0 <= chain.sphericity <= 1.0
        assert chain.shape_parameter > 0
        radius, center = chain.bsphere
        assert radius > 0
        assert center.shape == (3,)

    @pytest.mark.parametrize(
        ("points", "sphericity"),
        [
            # a rod: all mass on a line
            ([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], 0.0),
            # a flat square ring: two equal moments, the third zero
            ([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]], 0.75),
            # an octahedron: three equal moments, as for a sphere
            (
                [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]],
                1.0,
            ),
        ],
    )
    def test_sphericity_is_one_for_spherical_and_zero_for_rods(
        self, points: list, sphericity: float
    ):
        n_atoms = len(points)
        uni = Universe.empty(
            n_atoms, n_residues=1, atom_resindex=[0] * n_atoms, trajectory=True
        )
        uni.add_TopologyAttr("resids", [1])
        uni.add_TopologyAttr("masses", [1.0] * n_atoms)
        uni.add_TopologyAttr("bonds", [(i, i + 1) for i in range(n_atoms - 1)])
        uni.atoms.positions = np.array(points, dtype=float) + 10.0
        uni.dimensions = [20.0, 20.0, 20.0, 90.0, 90.0, 90.0]

        assert MolGroup(uni, [1]).sphericity == pytest.approx(sphericity, abs=1e-6)

    def test_dipole(self, chain: Cluster):
        # each residue is a +0.3/-0.3 pair 1.2 A apart along z; 1 e·A = 4.8032 D
        expected = 3 * 0.3 * 1.2 * 4.80320
        assert chain.dipole_moment == pytest.approx(expected, rel=1e-4)
        np.testing.assert_allclose(np.abs(chain.dipole), [0, 0, expected], atol=1e-3)

    def test_dipole_of_a_charged_group_is_about_its_center_of_mass(self):
        # +1 e on a light atom at x = 0 and a neutral heavy one at x = 2: the center
        # of mass is at x = 1.5, so the dipole about it is 1 e x -1.5 A
        uni = Universe.empty(2, n_residues=1, atom_resindex=[0, 0], trajectory=True)
        uni.add_TopologyAttr("resids", [1])
        uni.add_TopologyAttr("masses", [1.0, 3.0])
        uni.add_TopologyAttr("charges", [1.0, 0.0])
        uni.add_TopologyAttr("bonds", [(0, 1)])
        uni.atoms.positions = [[10.0, 10.0, 10.0], [12.0, 10.0, 10.0]]
        uni.dimensions = [20.0, 20.0, 20.0, 90.0, 90.0, 90.0]

        ion = MolGroup(uni, [1])

        np.testing.assert_allclose(ion.dipole, [-1.5 * 4.80320, 0, 0], rtol=1e-4)
        assert ion.dipole_moment == pytest.approx(1.5 * 4.80320, rel=1e-4)

    def test_center_of_mass_of_a_group_inside_the_box(self, chain: Cluster):
        np.testing.assert_allclose(
            chain.center_of_mass, chain.atoms.center_of_mass(), rtol=1e-6
        )


def read_geometry(group: MolGroup) -> None:
    """Read every property that needs whole positions."""
    _ = (
        group.sphericity,
        group.dipole_moment,
        group.dipole,
        group.shape_parameter,
        group.bsphere,
        group.radius_of_gyration,
        group.center_of_mass,
    )


def split_across_the_box_edge(
    make_universe: UniverseFactory,
) -> tuple[MolGroup, MolGroup, np.ndarray]:
    """The chain 1-2-3, whole in one Universe and split by the box edge in another.

    Returns
    -------
    tuple[MolGroup, MolGroup, np.ndarray]
        The whole chain, the split chain, and the split Universe's positions.
    """
    reference = MolGroup(make_universe([[[1, 2, 3]]], 3), [1, 2, 3])
    uni = make_universe([[[1, 2, 3]]], 3)
    box = 20.0
    # move the chain so its residues sit at y = 19, 21, 23, then wrap y into
    # the box: the chain is split, one residue at y = 19 and two at y = 1, 3
    shifted = uni.atoms.positions - [95.0, 81.0, 95.0]
    shifted[:, 1] %= box
    uni.dimensions = [box, box, box, 90.0, 90.0, 90.0]
    uni.atoms.positions = shifted
    return reference, MolGroup(uni, [1, 2, 3]), shifted


def chain_across_a_small_box(
    make_universe: UniverseFactory, n_res: int, box: float = 20.0
) -> tuple[MolGroup, Universe, np.ndarray]:
    """A straight chain of `n_res` residues along y, BOND_STEP apart, in a small box.

    The chain starts at y = 13 and is wrapped into a `box`-wide box, so from
    ``n_res = 5`` on it is split by the box edge, and from ``n_res = 7`` it reaches
    more than half a box from its first residue.

    Returns
    -------
    tuple[MolGroup, Universe, np.ndarray]
        The same chain whole in a huge box, the small-box Universe, and the shift
        from the huge box's positions to the small box's before wrapping.
    """
    residues = list(range(1, n_res + 1))
    reference = MolGroup(make_universe([[residues]], n_res), residues)
    uni = make_universe([[residues]], n_res)
    shift = np.array([-90.0, -87.0, -90.0])
    uni.dimensions = [box, box, box, 90.0, 90.0, 90.0]
    uni.atoms.positions = (uni.atoms.positions + shift) % box
    return reference, uni, shift


class TestWholePositions:
    """Geometric properties see whole groups without moving the shared Universe."""

    @pytest.mark.parametrize("kind", [Cluster, MolGroup])
    def test_group_reaching_past_half_the_box_is_made_whole(
        self, make_universe: UniverseFactory, kind: type
    ):
        # 14 A long in a 20 A box: its far end is more than half a box from the
        # first residue, and would be wrapped back onto the wrong side
        reference, uni, shift = chain_across_a_small_box(make_universe, 8)
        if kind is Cluster:
            (connected,) = connected_groups(uni)
            group = Cluster(uni, connected, cluster_id=1)
        else:
            group = MolGroup(uni, range(1, 9))

        with group.whole() as atoms:
            span = np.ptp(atoms.positions[:, 1])

        assert span == pytest.approx(14.0, abs=1e-4)
        assert group.radius_of_gyration == pytest.approx(
            reference.radius_of_gyration, rel=1e-5
        )
        expected = apply_PBC(reference.atoms.center_of_mass() + shift, uni.dimensions)
        np.testing.assert_allclose(group.center_of_mass, expected, atol=1e-3)

    def test_cluster_wrapping_around_the_box_is_reported(
        self, make_universe: UniverseFactory, captured_logs: list[str]
    ):
        # 10 residues 2 A apart in a 20 A box: the last connects back to the first
        # through the box edge, so the cluster is an endless chain with no shape
        _, uni, _ = chain_across_a_small_box(make_universe, 10)
        (connected,) = connected_groups(uni)
        cluster = Cluster(uni, connected, cluster_id=7)

        _ = cluster.radius_of_gyration
        _ = cluster.sphericity

        warnings = [m for m in captured_logs if "wraps around" in m]
        assert len(warnings) == 1
        assert "cluster 7" in warnings[0]

    def test_compact_clusters_are_not_reported(
        self, make_universe: UniverseFactory, captured_logs: list[str]
    ):
        _, uni, _ = chain_across_a_small_box(make_universe, 8)
        (connected,) = connected_groups(uni)

        _ = Cluster(uni, connected, cluster_id=1).radius_of_gyration

        assert not any("wraps around" in m for m in captured_logs)

    def test_group_reaching_past_half_a_triclinic_box_is_made_whole(
        self, make_universe: UniverseFactory, captured_logs: list[str]
    ):
        # 10 A long along y, in a 60-degree box whose height along y is 17.3 A
        residues = [1, 2, 3, 4, 5, 6]
        reference = MolGroup(make_universe([[residues]], 6), residues)
        uni = make_universe([[residues]], 6)
        dimensions = np.array([20.0, 20.0, 20.0, 90.0, 90.0, 60.0])
        uni.dimensions = dimensions
        uni.atoms.positions = apply_PBC(
            uni.atoms.positions - [90.0, 87.0, 90.0], dimensions
        )
        (connected,) = connected_groups(uni)
        cluster = Cluster(uni, connected, cluster_id=1)

        assert cluster.radius_of_gyration == pytest.approx(
            reference.radius_of_gyration, rel=1e-5
        )
        assert not any("wraps around" in m for m in captured_logs)

    def test_positions_getter_returns_a_copy(self, chain: Cluster):
        # whole() saves `atoms.positions` as the originals to restore, which only
        # works if MDAnalysis's getter copies (it indexes the timestep array with
        # an index array). If an MDAnalysis upgrade returns a view, fail here
        # instead of silently skipping the restore.
        ts_positions = chain.universe.trajectory.ts.positions
        positions = chain.atoms.positions
        before = ts_positions.copy()

        positions += 100.0

        assert not np.shares_memory(positions, ts_positions)
        np.testing.assert_array_equal(ts_positions, before)

    def test_properties_leave_the_universe_positions_untouched(self, chain: Cluster):
        before = chain.universe.atoms.positions

        read_geometry(chain)

        np.testing.assert_array_equal(chain.universe.atoms.positions, before)

    def test_overlapping_group_does_not_change_the_cluster_geometry(
        self, make_universe: UniverseFactory
    ):
        uni = make_universe([[[1, 2, 3, 4, 5, 6]]], 6)
        cluster = MolGroup(uni, [1, 2, 3, 4, 5, 6])
        rg = cluster.radius_of_gyration

        # a nucleus-like subset whose first residue isn't the cluster's
        read_geometry(MolGroup(uni, [4, 5, 6]))

        assert cluster.radius_of_gyration == pytest.approx(rg)
        fresh = MolGroup(uni, [1, 2, 3, 4, 5, 6])
        assert fresh.radius_of_gyration == pytest.approx(rg)

    def test_group_split_across_the_box_edge_is_made_whole(
        self, make_universe: UniverseFactory
    ):
        reference, split, shifted = split_across_the_box_edge(make_universe)

        with split.whole() as atoms:
            span = np.ptp(atoms.positions[:, 1])

        assert span < 20.0 / 2
        # float32 positions: wrapping and unwrapping costs a few ulps
        assert split.radius_of_gyration == pytest.approx(
            reference.radius_of_gyration, rel=1e-5
        )
        np.testing.assert_allclose(split.universe.atoms.positions, shifted)

    def test_center_of_mass_of_a_split_group_is_that_of_the_whole_group(
        self, make_universe: UniverseFactory
    ):
        reference, split, _ = split_across_the_box_edge(make_universe)
        # the whole chain's center of mass, moved like the chain and wrapped
        expected = apply_PBC(
            reference.atoms.center_of_mass() - [95.0, 81.0, 95.0],
            split.universe.dimensions,
        )

        np.testing.assert_allclose(split.center_of_mass, expected, atol=1e-4)
        # MDAnalysis's own puts it between the pieces, far from the chain
        assert abs(split.atoms.center_of_mass()[1] - expected[1]) > 5.0

    def test_whole_restores_positions_on_error(self, chain: Cluster):
        before = chain.universe.atoms.positions

        with pytest.raises(RuntimeError), chain.whole():
            raise RuntimeError

        np.testing.assert_array_equal(chain.universe.atoms.positions, before)

    def test_cache_follows_a_change_of_residues(self, make_universe: UniverseFactory):
        uni = make_universe([[[1, 2], [3, 4]]], 4)
        first = Cluster(uni, connected_groups(uni)[0], cluster_id=1)
        rg_before = first.radius_of_gyration

        (everything,) = connected_groups(uni, cutoff=60.0)
        first._update(everything)

        expected = MolGroup(uni, [1, 2, 3, 4]).radius_of_gyration
        assert first.radius_of_gyration == pytest.approx(expected)
        assert first.radius_of_gyration != pytest.approx(rg_before)


class TestFrameCache:
    """Geometric properties are computed once per frame and residue set."""

    @pytest.mark.parametrize(
        "name",
        [
            "radius",
            "radius_of_gyration",
            "sphericity",
            "dipole_moment",
            "shape_parameter",
            "center_of_mass",
            "dipole",
            "bsphere",
        ],
    )
    def test_value_is_computed_once_per_frame(
        self, chain: Cluster, name: str, monkeypatch: pytest.MonkeyPatch
    ):
        entries = []
        whole = MolGroup.whole

        def counted(group: MolGroup) -> AbstractContextManager:
            entries.append(group)
            return whole(group)

        monkeypatch.setattr(MolGroup, "whole", counted)

        first = getattr(chain, name)

        assert getattr(chain, name) is first
        assert len(entries) == 1

    def test_radius_is_computed_once_for_everything_built_on_it(
        self, chain: Cluster, monkeypatch: pytest.MonkeyPatch
    ):
        lookups = []
        vdw_radii = cluster._vdw_radii

        def counted(universe: Universe) -> np.ndarray:
            lookups.append(universe)
            return vdw_radii(universe)

        monkeypatch.setattr(cluster, "_vdw_radii", counted)

        _ = chain.radius, chain.diameter, chain.volume, chain.density, chain.radius

        assert len(lookups) == 1


class TestReadOnly:
    """What a group hands out can't change the group."""

    @pytest.mark.parametrize(
        "read",
        [
            lambda g: g.resids,
            lambda g: g.resnames,
            lambda g: g.center_of_mass,
            lambda g: g.dipole,
            lambda g: g.bsphere[1],
        ],
        ids=["resids", "resnames", "center_of_mass", "dipole", "bsphere"],
    )
    def test_arrays_are_read_only(self, chain: Cluster, read: Callable) -> None:
        array = read(chain)

        with pytest.raises(ValueError, match="read-only"):
            array[0] = array[1]

        assert array.copy().flags.writeable  # a copy is the caller's own

    def test_residue_indices_changed_in_place_leave_the_group_alone(
        self, chain: Cluster
    ):
        residues = chain.residues
        residues.ix[0] = 4

        assert members(chain) == [1, 2, 3]
        assert chain.residues is not residues

    def test_a_residue_group_it_was_built_from_leaves_it_alone(
        self, make_universe: UniverseFactory
    ):
        uni = make_universe([[]], 3)
        residues = uni.residues[[0, 2]]
        group = MolGroup(uni, residues)

        residues.ix[0] = 1

        assert members(group) == [1, 3]

    def test_universe_can_be_read_not_replaced(self, chain: Cluster):
        universe = chain.universe

        with pytest.raises(AttributeError):
            chain.universe = universe  # type: ignore[misc]
        assert chain.residues.universe is universe


class TestClusterGraph:
    def test_built_from_subconntable(self, chain: Cluster):
        assert chain.size == 3
        assert members(chain) == [1, 2, 3]
        assert sorted(chain) == [1, 2, 3]
        assert set(chain.graph[2]) == {1, 3}
        assert chain.distance(1, 2) == pytest.approx(BOND_STEP, abs=1e-3)
        assert "distance" in str(chain)

    def test_graph_is_frozen(self, chain: Cluster):
        with pytest.raises(nx.NetworkXError, match="Frozen"):
            chain.graph.add_edge(1, 5)

        assert members(chain) == [1, 2, 3]

    @pytest.mark.parametrize(
        "attrs",
        [
            lambda g: g[1][2],
            lambda g: g[2][1],
            lambda g: g.edges[1, 2],
            lambda g: g.nodes[1],
            lambda g: g.graph,
        ],
        ids=["edge", "reversed edge", "edges view", "node", "graph"],
    )
    def test_graph_attributes_are_frozen(self, chain: Cluster, attrs: Callable):
        distance = chain.distance(1, 2)

        with pytest.raises(TypeError):
            attrs(chain.graph)["distance"] = 0.0

        assert chain.distance(1, 2) == distance

    def test_graph_copies_can_be_modified(self, chain: Cluster):
        distance = chain.distance(1, 2)

        for copy in (nx.Graph(chain.graph), chain.graph.copy()):
            copy[1][2]["distance"] = 0.0
            copy.nodes[1]["label"] = "first"
            copy.add_edge(1, 5)

        assert chain.distance(1, 2) == distance
        assert sorted(chain) == [1, 2, 3]

    def test_subgraphs_read_the_attributes(self, chain: Cluster):
        sub = nx.induced_subgraph(chain.graph, [1, 2])

        assert sub[1][2]["distance"] == chain.distance(1, 2)

    def test_update_replaces_the_graph_and_residues(self, chain: Cluster):
        old_graph = chain.graph
        (everything,) = connected_groups(chain.universe, cutoff=60.0)

        chain._update(everything)

        assert members(chain) == [1, 2, 3, 4, 5]
        assert sorted(chain.graph) == [1, 2, 3, 4, 5]
        assert nx.is_frozen(chain.graph)
        u, v = next(iter(chain.graph.edges))
        with pytest.raises(TypeError):
            chain.graph[u][v]["distance"] = 0.0
        # a graph read before the update keeps describing the old frame
        assert sorted(old_graph) == [1, 2, 3]

    def test_neighbors_up_to_a_level(self, chain: Cluster):
        assert chain.neighbors(1) == {2}
        assert chain.neighbors(1, level=2) == {2, 3}

    @pytest.mark.parametrize(("moli", "molj"), [(5, 1), (1, 5)])
    def test_distance_rejects_mols_outside(self, chain: Cluster, moli: int, molj: int):
        with pytest.raises(ValueError, match="not in the cluster"):
            chain.distance(moli, molj)

    def test_distance_rejects_unconnected_mols(self, chain: Cluster):
        with pytest.raises(ValueError, match="not connected"):
            chain.distance(1, 3)

    def test_birth_time_and_age(self, make_universe: UniverseFactory):
        uni = make_universe([[[1, 2]], [[1, 2]]], 2)
        cls = Cluster(uni, connected_groups(uni)[0], cluster_id=1)

        uni.trajectory[1]

        assert cls.birth_time == 0.0
        assert cls.age == 1.0

    def test_equality_and_hash_are_by_identity(self, chain: Cluster):
        (group,) = connected_groups(chain.universe)
        twin = Cluster(chain.universe, group, cluster_id=1)
        seen = {chain}

        assert twin != chain
        assert twin not in seen

        # the same object, whatever its later frames hold
        (everything,) = connected_groups(chain.universe, cutoff=60.0)
        chain._update(everything)
        assert chain == chain
        assert chain in seen
