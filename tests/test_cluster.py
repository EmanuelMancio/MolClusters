# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

import numpy as np
import pytest
from MDAnalysis import Universe

from molclusters.cluster import EA2D, Cluster, MDAResidueGroupAnalyzer
from molclusters.config import MolClsConfig
from molclusters.conntable import ConnectionTable

from .conftest import BOND_STEP, CUTOFF, UniverseFactory


@pytest.fixture
def uni() -> Universe:
    """Fresh Universe for the met-mal fixture data, isolated per test.

    Returns
    -------
    Universe
        A newly loaded Universe from the met-mal topology/structure files.
    """
    return Universe(
        "tests/data/met-mal/met-mal.tpr",
        "tests/data/met-mal/start.pdb",
    )


@pytest.fixture
def empty_cluster(uni: Universe) -> Cluster:
    """A fresh, empty Cluster, independent per test.

    Returns
    -------
    Cluster
        A newly constructed Cluster with no molecules.
    """
    return Cluster(uni)


@pytest.fixture
def populated_cluster(uni: Universe) -> Cluster:
    """Build a real two-molecule Cluster via the same path production code uses.

    Returns
    -------
    Cluster
        A cluster containing the two connected MOL residues from the fixture.
    """
    config = MolClsConfig(rules={"MOL": {"MOL": "cm 15.0"}})
    sels = {res: uni.select_atoms(f"resname {res}") for res in config._rules.all_keys()}
    conntab = ConnectionTable(uni, config._rules, sels)

    return Cluster(uni, next(conntab.subconntables()))


class TestCluster:
    def test_cluster_creation(self, empty_cluster: Cluster):
        assert empty_cluster is not None

    def test_cluster_length(self, empty_cluster: Cluster):
        assert len(empty_cluster) == 0

    def test_id_increments_per_instance(self, uni: Universe):
        first = Cluster(uni)
        second = Cluster(uni)

        assert second.id == first.id + 1

    def test_add_con_exception(self, empty_cluster: Cluster):
        with pytest.raises(ValueError):
            empty_cluster.add_con(0, 2, 2.0)

    def test_set_dist_exception_moli_not_in_cluster(self, empty_cluster: Cluster):
        with pytest.raises(ValueError):
            empty_cluster.set_dist(0, 1, 2.0)

    def test_add_unsupported_operand(self, empty_cluster: Cluster):
        with pytest.raises(TypeError):
            empty_cluster + "not a residue group"

    def test_print_cluster_np(self, empty_cluster: Cluster):
        print(np.array(empty_cluster.cluster))


class TestPopulatedCluster:
    def test_distance_between_mols(self, populated_cluster: Cluster):
        moli = next(iter(populated_cluster))
        molj = next(iter(populated_cluster.neighbors(moli)))

        expected = populated_cluster.cluster[moli][molj]["distance"]
        assert populated_cluster.get_dist(moli, molj) == expected

    def test_contains(self, populated_cluster: Cluster):
        mol = next(iter(populated_cluster))

        assert mol in populated_cluster
        assert -1 not in populated_cluster


@pytest.fixture
def chain(make_universe: UniverseFactory) -> Cluster:
    """A Cluster over the chain 1-2-3, with residues 4 and 5 left outside it.

    Returns
    -------
    Cluster
        Cluster built from the synthetic system's only connected component.
    """
    uni = make_universe([[[1, 2, 3]]], 5)
    config = MolClsConfig(rules={"MOL": {"MOL": f"cm {CUTOFF}"}})
    sels = {"MOL": uni.select_atoms("resname MOL")}
    conntab = ConnectionTable(uni, config._rules, sels)

    return Cluster(uni, next(conntab.subconntables()))


def members(cls: MDAResidueGroupAnalyzer) -> list[int]:
    return sorted(int(r) for r in cls.resids)


class TestResidueGroupAnalyzer:
    def test_built_from_resids(self, make_universe: UniverseFactory):
        uni = make_universe([[]], 3)

        analyzer = MDAResidueGroupAnalyzer(uni, [1, 3])

        assert members(analyzer) == [1, 3]
        assert analyzer.size == len(analyzer) == 2
        assert list(analyzer.resnames) == ["MOL", "MOL"]

    def test_add_residue_group_and_analyzer(self, make_universe: UniverseFactory):
        uni = make_universe([[]], 3)
        first = MDAResidueGroupAnalyzer(uni, [1])

        with_group = first + uni.residues[[1]]
        with_analyzer = first + MDAResidueGroupAnalyzer(uni, [3])

        assert members(with_group) == [1, 2]
        assert members(with_analyzer) == [1, 3]

    def test_physical_properties_are_consistent(self, chain: Cluster):
        assert chain.mass == pytest.approx(3 * (12.011 + 15.999))
        assert chain.charge == pytest.approx(0.0)
        assert chain.radius == chain.radius_of_gyration > 0
        assert chain.diameter == pytest.approx(2 * chain.radius)
        assert chain.volume == pytest.approx(4 / 3 * np.pi * chain.radius**3)
        assert chain.density == pytest.approx(chain.mass / chain.volume * 0.602214076)

    def test_shape_properties(self, chain: Cluster):
        assert 0.0 <= chain.sphericity <= 1.0
        assert chain.shape_parameter > 0
        radius, center = chain.bsphere
        assert radius > 0
        assert center.shape == (3,)

    def test_dipole(self, chain: Cluster):
        # each residue is a +0.3/-0.3 pair 1.2 A apart along z
        expected = 3 * 0.3 * 1.2 * EA2D
        assert chain.dipole_moment == pytest.approx(expected, rel=1e-4)
        np.testing.assert_allclose(np.abs(chain.dipole), [0, 0, expected], atol=1e-3)


class TestClusterGraph:
    def test_built_from_subconntable(self, chain: Cluster):
        assert chain.size == 3
        assert members(chain) == [1, 2, 3]
        assert set(chain[2]) == {1, 3}
        assert chain.get_dist(1, 2) == pytest.approx(BOND_STEP, abs=1e-3)
        assert "distance" in str(chain)

    def test_center_of_mass_setter(self, chain: Cluster):
        chain.cm = np.array([1.0, 2.0, 3.0])

        np.testing.assert_array_equal(chain.cm, [1.0, 2.0, 3.0])

    def test_add_mol(self, chain: Cluster):
        chain.add_mol(3, 4, "MOL", 2.5)

        assert members(chain) == [1, 2, 3, 4]
        assert chain.get_dist(3, 4) == 2.5
        assert chain.cluster[3][4]["weight"] == pytest.approx(np.exp(1 / 2.5))

    @pytest.mark.parametrize(("ref", "mol"), [(5, 4), (3, 2)])
    def test_add_mol_rejects_bad_members(self, chain: Cluster, ref: int, mol: int):
        with pytest.raises(ValueError):
            chain.add_mol(ref, mol, "MOL", 2.0)

    def test_add_con_rejects_zero_distance(self, chain: Cluster):
        with pytest.raises(ValueError, match="zero"):
            chain.add_con(1, 3, 0.0)

    def test_add_con_rejects_second_mol_outside(self, chain: Cluster):
        with pytest.raises(ValueError):
            chain.add_con(1, 5, 2.0)

    def test_remove_mol(self, chain: Cluster):
        chain.remove_mol(3)

        assert members(chain) == [1, 2]
        with pytest.raises(ValueError):
            chain.remove_mol(3)

    def test_remove_con_and_cons(self, chain: Cluster):
        chain.remove_con(1, 2)
        assert chain.neighbors(2) == {3}

        chain.remove_cons(2, [3])
        assert chain.neighbors(2) == set()

    @pytest.mark.parametrize(("moli", "molj"), [(5, 1), (1, 5)])
    def test_remove_con_rejects_mols_outside(
        self, chain: Cluster, moli: int, molj: int
    ):
        with pytest.raises(ValueError):
            chain.remove_con(moli, molj)

    def test_neighbors_up_to_a_level(self, chain: Cluster):
        assert chain.neighbors(1) == {2}
        assert chain.neighbors(1, level=2) == {2, 3}

    @pytest.mark.parametrize(("moli", "molj"), [(5, 1), (1, 5)])
    def test_get_dist_rejects_mols_outside(self, chain: Cluster, moli: int, molj: int):
        with pytest.raises(ValueError):
            chain.get_dist(moli, molj)

    def test_set_dist(self, chain: Cluster):
        chain.set_dist(1, 2, 2.2)

        assert chain.get_dist(1, 2) == 2.2

    @pytest.mark.parametrize(("moli", "molj", "dist"), [(1, 5, 1.0), (1, 2, 0.0)])
    def test_set_dist_rejects_invalid(
        self, chain: Cluster, moli: int, molj: int, dist: float
    ):
        with pytest.raises(ValueError):
            chain.set_dist(moli, molj, dist)

    def test_separate_single_component_is_a_noop(self, chain: Cluster):
        assert chain.separate() == []
        assert members(chain) == [1, 2, 3]

    def test_separate_keeps_largest_component(self, chain: Cluster):
        chain.remove_con(2, 3)

        (split,) = chain.separate()

        assert members(chain) == [1, 2]
        assert members(split) == [3]
        assert split.id > chain.id

    @pytest.fixture
    def pair_of_clusters(self, make_universe: UniverseFactory) -> list[Cluster]:
        uni = make_universe([[[1, 2], [3, 4]]], 4)
        config = MolClsConfig(rules={"MOL": {"MOL": f"cm {CUTOFF}"}})
        conntab = ConnectionTable(uni, config._rules, {"MOL": uni.atoms})

        return [Cluster(uni, sub) for sub in conntab.subconntables()]

    def test_merge_joins_the_graphs(self, pair_of_clusters: list[Cluster]):
        first, second = pair_of_clusters

        first.merge(second)

        assert sorted(first) == [1, 2, 3, 4]
        assert first.size == 4

    def test_merge_updates_residue_group(self, pair_of_clusters: list[Cluster]):
        first, second = pair_of_clusters

        first.merge(second)

        assert members(first) == [1, 2, 3, 4]

    def test_update_from_conntable(self, chain: Cluster):
        uni = chain.uni
        config = MolClsConfig(rules={"MOL": {"MOL": "cm 60.0"}})
        conntab = ConnectionTable(uni, config._rules, {"MOL": uni.atoms})

        chain.update_from_conntable(next(conntab.subconntables()))

        assert members(chain) == [1, 2, 3, 4, 5]

    def test_age(self, make_universe: UniverseFactory):
        uni = make_universe([[], []], 2)
        cls = Cluster(uni)

        uni.trajectory[1]

        assert cls.get_age() == 1.0

    def test_equality_is_by_graph(self, chain: Cluster):
        twin = Cluster(chain.uni)
        twin.cluster = chain.cluster.copy()

        assert twin == chain
        assert chain != Cluster(chain.uni)
        assert chain != "not a cluster"
