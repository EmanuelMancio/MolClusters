# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

import re
from pathlib import Path

import MDAnalysis as mda
import networkx as nx
import numpy as np
import pytest
from MDAnalysis import Universe
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

from molclusters.config import DistanceBackend, MolClsConfig
from molclusters.conntable import (
    ConnectionTable,
    _check_hb_private_api,
    _warn_if_openmp_unavailable,
    check_resids,
)

from .conftest import BOND_STEP, CUTOFF, UniverseFactory

DATA_DIR = (Path(__file__).parent / "data" / "met-mal").resolve()


def dock_hbonds(uni: Universe, resname: str, n_pairs: int = 3) -> None:
    """Move molecules so that pairs of them H-bond with each other.

    met-mal's start.pdb has no H-bonds between two MAL molecules (only within
    them), so each pair's second molecule is moved to put an acceptor 1.9 A from
    the first molecule's hydrogen, straight along the donor-hydrogen bond.

    Parameters
    ----------
    uni : Universe
        The system, changed in place.
    resname : str
        The residue name of the molecules to pair up, in topology order.
    n_pairs : int
        The number of pairs to dock.
    """
    # in memory, so that the moves survive a trajectory rewind (e.g. run())
    uni.transfer_to_memory()
    guesser = HydrogenBondAnalysis(uni)
    hydrogens = uni.select_atoms(guesser.guess_hydrogens())
    acceptors = uni.select_atoms(guesser.guess_acceptors())
    residues = uni.select_atoms(f"resname {resname}").residues
    pairs = zip(residues[0::2][:n_pairs], residues[1::2][:n_pairs], strict=True)
    for first, second in pairs:
        hydrogen = (hydrogens & first.atoms)[0]
        donor = hydrogen.bonded_atoms[0]
        acceptor = (acceptors & second.atoms)[0]
        bond = hydrogen.position - donor.position
        target = hydrogen.position + 1.9 * bond / np.linalg.norm(bond)
        second.atoms.positions += target - acceptor.position


def build(
    uni: Universe, rules: dict, backend: DistanceBackend = "serial"
) -> ConnectionTable:
    config = MolClsConfig(rules=rules, distance_backend=backend)
    sels = {res: uni.select_atoms(f"resname {res}") for res in config._rules.all_keys()}
    return ConnectionTable(uni, config._rules, sels, backend=config.distance_backend)


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

    def test_molecule_split_across_the_box_still_connects(
        self, make_universe: UniverseFactory
    ):
        uni = make_universe([[]], 2)
        top = uni.dimensions[2]
        # residue 1 sits at the top of the box with its O atom (1.2 A above its C)
        # wrapped to the bottom; residue 2 is whole, just below it
        uni.atoms.positions = [
            [100.0, 100.0, top - 0.5],
            [100.0, 100.0, 0.7],
            [100.0, 102.0, top - 2.0],
            [100.0, 102.0, top - 0.8],
        ]

        table = build(uni, {"MOL": {"MOL": f"cm {CUTOFF}"}})

        # the whole residues' centers of mass are 2 A apart in y and 1.5 A in z
        assert sorted(table) == [1, 2]
        assert table[1, 2] == pytest.approx(2.5, abs=1e-3)

    def test_no_molecules_within_cutoff(self, make_universe: UniverseFactory):
        uni = make_universe([[]], 3)

        table = build(uni, {"MOL": {"MOL": f"cm {CUTOFF}"}})

        assert len(list(table)) == 0
        assert table.rule_connections["MOL", "MOL"] == 0


class TestDistanceBackend:
    def test_defaults_to_serial(self, chain_table: ConnectionTable):
        assert chain_table.backend == "serial"

    def test_backend_is_forwarded_to_self_capped_distance(
        self, make_universe: UniverseFactory, monkeypatch: pytest.MonkeyPatch
    ):
        uni = make_universe([[[1, 2, 3]]], 3)
        seen = {}
        original = mda.lib.distances.self_capped_distance

        def spy(*args: object, **kwargs: object) -> object:
            seen["backend"] = kwargs.get("backend")
            return original(*args, **kwargs)

        monkeypatch.setattr(mda.lib.distances, "self_capped_distance", spy)

        table = build(uni, {"MOL": {"MOL": f"cm {CUTOFF}"}}, backend="OpenMP")

        assert seen["backend"] == "OpenMP"
        assert sorted(table) == [1, 2, 3]

    def test_backend_is_forwarded_to_capped_distance(
        self, make_universe: UniverseFactory, monkeypatch: pytest.MonkeyPatch
    ):
        uni = make_universe([[[1, 2]]], 2, ["MOL", "SOL"])
        seen = {}
        original = mda.lib.distances.capped_distance

        def spy(*args: object, **kwargs: object) -> object:
            seen["backend"] = kwargs.get("backend")
            return original(*args, **kwargs)

        monkeypatch.setattr(mda.lib.distances, "capped_distance", spy)

        table = build(uni, {"MOL": {"SOL": f"cm {CUTOFF}"}}, backend="OpenMP")

        assert seen["backend"] == "OpenMP"
        assert sorted(table) == [1, 2]


class TestOpenmpAvailabilityWarning:
    def test_warns_when_openmp_was_compiled_without_it(
        self, monkeypatch: pytest.MonkeyPatch, captured_logs: list[str]
    ):
        monkeypatch.setattr(mda.lib.distances, "USED_OPENMP", False)

        _warn_if_openmp_unavailable("OpenMP")

        assert any("OpenMP" in m and "serial" in m for m in captured_logs)

    def test_no_warning_when_openmp_is_actually_available(
        self, monkeypatch: pytest.MonkeyPatch, captured_logs: list[str]
    ):
        monkeypatch.setattr(mda.lib.distances, "USED_OPENMP", True)

        _warn_if_openmp_unavailable("OpenMP")

        assert captured_logs == []

    def test_no_warning_for_serial_backend(
        self, monkeypatch: pytest.MonkeyPatch, captured_logs: list[str]
    ):
        monkeypatch.setattr(mda.lib.distances, "USED_OPENMP", False)

        _warn_if_openmp_unavailable("serial")

        assert captured_logs == []

    def test_connection_table_construction_triggers_the_warning(
        self,
        make_universe: UniverseFactory,
        monkeypatch: pytest.MonkeyPatch,
        captured_logs: list[str],
    ):
        monkeypatch.setattr(mda.lib.distances, "USED_OPENMP", False)
        uni = make_universe([[[1, 2]]], 2)

        build(uni, {"MOL": {"MOL": f"cm {CUTOFF}"}}, backend="OpenMP")

        assert any("OpenMP" in m for m in captured_logs)


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
        dock_hbonds(uni, "MAL")
        return build(uni, {"MAL": {"MAL": f"hb d {self.D_A} a {self.ANGLE}"}})

    def test_edges_carry_distance_and_angle(self, table: ConnectionTable):
        edges = list(table.conntab.edges(data=True))

        assert edges, "met-mal should have MAL-MAL hydrogen bonds"
        for _, _, data in edges:
            assert 0 < data["distance"] <= self.D_A
            assert self.ANGLE <= data["angle"] <= 180

    def test_hbonds_within_a_molecule_are_not_connections(self, table: ConnectionTable):
        hb = HydrogenBondAnalysis(
            table.uni,
            between=["resname MAL", "resname MAL"],
            d_a_cutoff=self.D_A,
            d_h_a_angle_cutoff=self.ANGLE,
        )
        hb.run()
        resids = table.uni.atoms.resids
        hbonds = hb.results.hbonds.astype(np.intp)
        within = resids[hbonds[:, 2]] == resids[hbonds[:, 3]]
        assert within.any(), "MAL has an H-bond within the molecule"

        assert nx.number_of_selfloops(table.conntab) == 0
        assert table.rule_connections["MAL", "MAL"] == np.count_nonzero(~within)

    @pytest.fixture
    def mixed_uni(self) -> Universe:
        """met-mal with every other MAL renamed MAX, and some MAL pairs H-bonding.

        Returns
        -------
        Universe
            A system with H-bonding residues outside a MAL-MAL rule.
        """
        uni = Universe(str(DATA_DIR / "met-mal.tpr"), str(DATA_DIR / "start.pdb"))
        resnames = uni.residues.resnames.copy()
        resnames[np.flatnonzero(resnames == "MAL")[::2]] = "MAX"
        uni.residues.resnames = resnames
        dock_hbonds(uni, "MAL")
        return uni

    def public_run_edges(self, uni: Universe) -> set[frozenset[int]]:
        hb = HydrogenBondAnalysis(
            uni,
            between=["resname MAL", "resname MAL"],
            d_a_cutoff=self.D_A,
            d_h_a_angle_cutoff=self.ANGLE,
        )
        hb.run()
        atoms = uni.atoms
        pairs = {
            frozenset((atoms[int(h)].resid, atoms[int(a)].resid))
            for h, a in hb.results.hbonds[:, [2, 3]]
        }
        # H-bonds within one molecule connect it to nothing
        return {pair for pair in pairs if len(pair) == 2}

    def test_matches_mdanalysis_public_run(self, table: ConnectionTable):
        expected = self.public_run_edges(table.uni)

        assert {frozenset(e) for e in table.conntab.edges} == expected

    def test_mixed_system_matches_mdanalysis_public_run(self, mixed_uni: Universe):
        table = build(mixed_uni, {"MAL": {"MAL": f"hb d {self.D_A} a {self.ANGLE}"}})

        expected = self.public_run_edges(mixed_uni)

        assert expected, "the remaining MAL should still H-bond"
        assert {frozenset(e) for e in table.conntab.edges} == expected

    def test_only_searches_the_rules_residues(self, mixed_uni: Universe):
        table = build(mixed_uni, {"MAL": {"MAL": f"hb d {self.D_A} a {self.ANGLE}"}})
        hb = table.hbs["MAL", "MAL"]

        # what MDAnalysis guesses over the whole system, minus the MAX residues
        mal = mixed_uni.select_atoms("resname MAL")
        hydrogens = mixed_uni.select_atoms(hb.guess_hydrogens()) & mal
        acceptors = mixed_uni.select_atoms(hb.guess_acceptors()) & mal

        assert hydrogens and acceptors
        np.testing.assert_array_equal(hb._hydrogens.indices, hydrogens.indices)
        np.testing.assert_array_equal(hb._acceptors.indices, acceptors.indices)
        assert set(hb._donors.resnames) == {"MAL"}

    def test_rule_sets_the_hydrogen_and_acceptor_criteria(self, mixed_uni: Universe):
        spec = f"hb d {self.D_A} a {self.ANGLE} hmin 1.0 hmax 1.2 hq 0.46 aq -0.4"
        hb = build(mixed_uni, {"MAL": {"MAL": spec}}).hbs["MAL", "MAL"]

        mal = mixed_uni.select_atoms("resname MAL")
        hydrogens = (
            mixed_uni.select_atoms(
                hb.guess_hydrogens(min_mass=1.0, max_mass=1.2, min_charge=0.46)
            )
            & mal
        )
        acceptors = mixed_uni.select_atoms(hb.guess_acceptors(max_charge=-0.4)) & mal

        # these thresholds drop HMA and add OMC/OMD (the defaults' OM2/OM9 stay)
        assert len(hydrogens) == 2 * mal.n_residues
        assert len(acceptors) == 4 * mal.n_residues
        np.testing.assert_array_equal(hb._hydrogens.indices, hydrogens.indices)
        np.testing.assert_array_equal(hb._acceptors.indices, acceptors.indices)

    def test_rule_without_hbonding_atoms_finds_nothing(self):
        uni = Universe(str(DATA_DIR / "met-mal.tpr"), str(DATA_DIR / "start.pdb"))

        table = build(uni, {"MOL": {"MOL": f"hb d {self.D_A} a {self.ANGLE}"}})

        assert table.conntab.number_of_edges() == 0
        assert table.rule_connections["MOL", "MOL"] == 0

    def test_update_is_repeatable(self, table: ConnectionTable):
        before = sorted(map(sorted, table.conntab.edges))

        table.update()

        assert sorted(map(sorted, table.conntab.edges)) == before


def test_missing_private_hb_api_fails_fast():
    with pytest.raises(RuntimeError, match="_prepare"):
        _check_hb_private_api(object())


class TestCheckResids:
    def test_residues_numbered_one_to_n_pass(self, make_universe: UniverseFactory):
        check_resids(make_universe([[]], 3))

    @pytest.mark.parametrize(
        ("resids", "message"),
        [
            ([0, 1, 2], "residue 1 has resid 0 (3 residue(s)"),  # an offset
            ([1, 2, 5], "residue 3 has resid 5 (1 residue(s)"),  # a gap
            ([1, 1, 3], "residue 2 has resid 1 (1 residue(s)"),  # a repeat
        ],
    )
    def test_other_numberings_are_refused(
        self, make_universe: UniverseFactory, resids: list[int], message: str
    ):
        uni = make_universe([[]], 3)
        uni.residues.resids = resids

        with pytest.raises(ValueError, match=re.escape(message)):
            check_resids(uni)

    def test_the_table_checks_them(self, make_universe: UniverseFactory):
        uni = make_universe([[[1, 2]]], 2)
        uni.residues.resids = [5, 6]

        with pytest.raises(ValueError, match="numbered 1 to 2"):
            build(uni, {"MOL": {"MOL": f"cm {CUTOFF}"}})
