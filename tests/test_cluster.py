# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

from collections import Counter

import networkx as nx
import numpy as np
import pytest
from MDAnalysis import Universe
from MDAnalysis.lib.distances import apply_PBC

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
        combined = chain + chain.uni.residues[[3]]

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


class TestWholePositions:
    """Geometric properties see whole groups without moving the shared Universe."""

    def test_positions_getter_returns_a_copy(self, chain: Cluster):
        # whole() saves `atoms.positions` as the originals to restore, which only
        # works if MDAnalysis's getter copies (it indexes the timestep array with
        # an index array). If an MDAnalysis upgrade returns a view, fail here
        # instead of silently skipping the restore.
        ts_positions = chain.uni.trajectory.ts.positions
        positions = chain.atoms.positions
        before = ts_positions.copy()

        positions += 100.0

        assert not np.shares_memory(positions, ts_positions)
        np.testing.assert_array_equal(ts_positions, before)

    def test_properties_leave_the_universe_positions_untouched(self, chain: Cluster):
        before = chain.uni.atoms.positions

        read_geometry(chain)

        np.testing.assert_array_equal(chain.uni.atoms.positions, before)

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
        np.testing.assert_allclose(split.uni.atoms.positions, shifted)

    def test_center_of_mass_of_a_split_group_is_that_of_the_whole_group(
        self, make_universe: UniverseFactory
    ):
        reference, split, _ = split_across_the_box_edge(make_universe)
        # the whole chain's center of mass, moved like the chain and wrapped
        expected = apply_PBC(
            reference.atoms.center_of_mass() - [95.0, 81.0, 95.0],
            split.uni.dimensions,
        )

        np.testing.assert_allclose(split.center_of_mass, expected, atol=1e-4)
        # MDAnalysis's own puts it between the pieces, far from the chain
        assert abs(split.atoms.center_of_mass()[1] - expected[1]) > 5.0

    def test_whole_restores_positions_on_error(self, chain: Cluster):
        before = chain.uni.atoms.positions

        with pytest.raises(RuntimeError), chain.whole():
            raise RuntimeError

        np.testing.assert_array_equal(chain.uni.atoms.positions, before)

    def test_cache_follows_a_change_of_residues(self, make_universe: UniverseFactory):
        uni = make_universe([[[1, 2], [3, 4]]], 4)
        first = Cluster(uni, connected_groups(uni)[0], cluster_id=1)
        rg_before = first.radius_of_gyration

        (everything,) = connected_groups(uni, cutoff=60.0)
        first._update(everything)

        expected = MolGroup(uni, [1, 2, 3, 4]).radius_of_gyration
        assert first.radius_of_gyration == pytest.approx(expected)
        assert first.radius_of_gyration != pytest.approx(rg_before)


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

    def test_update_replaces_the_graph_and_residues(self, chain: Cluster):
        old_graph = chain.graph
        (everything,) = connected_groups(chain.uni, cutoff=60.0)

        chain._update(everything)

        assert members(chain) == [1, 2, 3, 4, 5]
        assert sorted(chain.graph) == [1, 2, 3, 4, 5]
        assert nx.is_frozen(chain.graph)
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
        (group,) = connected_groups(chain.uni)
        twin = Cluster(chain.uni, group, cluster_id=1)
        seen = {chain}

        assert twin != chain
        assert twin not in seen

        # the same object, whatever its later frames hold
        (everything,) = connected_groups(chain.uni, cutoff=60.0)
        chain._update(everything)
        assert chain == chain
        assert chain in seen
