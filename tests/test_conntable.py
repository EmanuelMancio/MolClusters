# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

from pathlib import Path

import numpy as np
import pytest
from MDAnalysis import Universe
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

from molclusters.config import MolClsConfig
from molclusters.conntable import ConnectionTable, _check_hb_private_api

from .conftest import BOND_STEP, CUTOFF, UniverseFactory

DATA_DIR = (Path(__file__).parent / "data" / "met-mal").resolve()


def build(uni: Universe, rules: dict) -> ConnectionTable:
    config = MolClsConfig(rules=rules)
    sels = {res: uni.select_atoms(f"resname {res}") for res in config._rules.all_keys()}
    return ConnectionTable(uni, config._rules, sels)


@pytest.fixture
def chain_table(make_universe: UniverseFactory) -> ConnectionTable:
    """Connections for the chain 1-2-3 plus the pair 4-5; residue 6 is isolated.

    Returns
    -------
    ConnectionTable
        Table built with a single MOL-MOL ``cm`` rule.
    """
    uni = make_universe([[[1, 2, 3], [4, 5]]], 6)
    return build(uni, {"MOL": {"MOL": f"cm {CUTOFF}"}})


class TestCenterOfMassRule:
    def test_edges_follow_the_cutoff(self, chain_table: ConnectionTable):
        assert sorted(chain_table) == [1, 2, 3, 4, 5]
        assert 6 not in chain_table
        assert sorted(chain_table.mols_connected_to(2)) == [1, 3]
        assert chain_table[1, 2] == pytest.approx(BOND_STEP, abs=1e-3)

    def test_rule_between_different_resnames(self, make_universe: UniverseFactory):
        resnames = ["MOL", "SOL", "SOL"]
        uni = make_universe([[[1, 2, 3]]], 3, resnames)

        table = build(uni, {"MOL": {"SOL": f"cm {CUTOFF}"}})

        # 2-3 are close too, but no SOL-SOL rule exists
        assert sorted(table.connections_from(1)) == [(1, 2)]
        assert 3 not in table

    def test_update_follows_the_trajectory(self, make_universe: UniverseFactory):
        uni = make_universe([[[1, 2]], [[2, 3]]], 3)
        table = build(uni, {"MOL": {"MOL": f"cm {CUTOFF}"}})

        uni.trajectory[1]
        table.update()

        assert sorted(table) == [2, 3]


class TestLookups:
    def test_item_by_molecule_lists_neighbours(self, chain_table: ConnectionTable):
        assert sorted(chain_table[2]) == [1, 3]

    @pytest.mark.parametrize("key", [6, (6, 1), (1, 6), (1, 2, 3)])
    def test_invalid_keys_raise(
        self, chain_table: ConnectionTable, key: int | tuple[int, ...]
    ):
        with pytest.raises(KeyError):
            chain_table[key]

    def test_connection_trees(self, chain_table: ConnectionTable):
        assert sorted(chain_table.mols_connected_tree_to(1)) == [1, 2, 3]
        assert len(chain_table.connection_tree_from(1)) == 2
        assert chain_table.connections_from(5) == [(5, 4)]

    def test_subconntables_are_sorted_by_size(self, chain_table: ConnectionTable):
        subs = list(chain_table.subconntables())

        assert [sorted(s) for s in subs] == [[1, 2, 3], [4, 5]]
        assert len(subs[0]) == 3
        assert list(subs[0].resnames) == ["MOL"] * 3
        assert sorted(subs[0][2]) == [1, 3]
        np.testing.assert_allclose(
            subs[1].cm, chain_table.uni.residues[[3, 4]].center_of_mass()
        )


class TestHydrogenBondRule:
    # MOL carries no partial charges in this topology, so only MAL can H-bond
    D_A, ANGLE = 3.5, 120.0

    @pytest.fixture
    def table(self) -> ConnectionTable:
        uni = Universe(str(DATA_DIR / "met-mal.tpr"), str(DATA_DIR / "start.pdb"))
        return build(uni, {"MAL": {"MAL": f"hb d {self.D_A} a {self.ANGLE}"}})

    def test_edges_carry_distance_and_angle(self, table: ConnectionTable):
        edges = list(table.conntab.edges(data=True))

        assert edges, "met-mal should have MAL-MAL hydrogen bonds"
        for _, _, data in edges:
            assert 0 < data["distance"] <= self.D_A
            assert self.ANGLE <= data["angle"] <= 180

    def test_matches_mdanalysis_public_run(self, table: ConnectionTable):
        hb = HydrogenBondAnalysis(
            table.uni,
            between=["resname MAL", "resname MAL"],
            d_a_cutoff=self.D_A,
            d_h_a_angle_cutoff=self.ANGLE,
        )
        hb.run()
        atoms = table.uni.atoms
        expected = {
            frozenset((atoms[int(h)].resid, atoms[int(a)].resid))
            for h, a in hb.results.hbonds[:, [2, 3]]
        }

        assert {frozenset(e) for e in table.conntab.edges} == expected

    def test_update_is_repeatable(self, table: ConnectionTable):
        before = sorted(map(sorted, table.conntab.edges))

        table.update()

        assert sorted(map(sorted, table.conntab.edges)) == before


def test_missing_private_hb_api_fails_fast():
    with pytest.raises(RuntimeError, match="_prepare"):
        _check_hb_private_api(object())
