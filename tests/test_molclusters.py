# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

import json
from collections import Counter
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from MDAnalysis import Universe

from molclusters.cluster import MDAResidueGroupAnalyzer
from molclusters.config import MolClsConfig
from molclusters.molclusters import MolClusters

from .conftest import CUTOFF, Groups, UniverseFactory

DATA_DIR = (Path(__file__).parent / "data" / "met-mal").resolve()

MOL_RULES = {"MOL": {"MOL": f"cm {CUTOFF}"}}

type Analyze = Callable[..., MolClusters]


def step(molcls: MolClusters, frame: int) -> None:
    """Advance to `frame` and reconcile clusters, as one iteration of `run()` does."""
    molcls.uni.trajectory[frame]
    molcls._MolClusters__update_clusters()


def snapshot(molcls: MolClusters) -> dict[int, set[int]]:
    """Cluster id -> set of resids, as plain ints for readable assertions.

    Returns
    -------
    dict[int, set[int]]
        The current clusters.
    """
    return {cid: {int(m) for m in cls} for cid, cls in molcls.clusters.items()}


def id_of(molcls: MolClusters, members: set[int]) -> int:
    """Return the id of the cluster made of exactly `members`.

    Returns
    -------
    int
        The matching cluster id.
    """
    matches = [cid for cid, mols in snapshot(molcls).items() if mols == members]
    assert len(matches) == 1, f"no single cluster {members} in {snapshot(molcls)}"
    return matches[0]


@pytest.fixture
def analyze(make_universe: UniverseFactory) -> Analyze:
    """Build a MolClusters over synthetic frames (frame 0 is loaded on creation).

    Returns
    -------
    Analyze
        ``analyze(frames, n_res, resnames=None, **config)``.
    """

    def factory(
        frames: Sequence[Groups],
        n_res: int,
        resnames: Sequence[str] | None = None,
        **config_kwargs: Any,  # noqa: ANN401
    ) -> MolClusters:
        uni = make_universe(frames, n_res, resnames)
        config_kwargs.setdefault("rules", MOL_RULES)
        return MolClusters(uni, MolClsConfig(**config_kwargs))

    return factory


class TestClusterIdentity:
    """The dominance algorithm: which cluster keeps its id from one frame to the next."""

    def test_initial_clusters_are_the_connected_components(self, analyze: Analyze):
        molcls = analyze([[[1, 2, 3], [4, 5]]], 6)

        assert sorted(snapshot(molcls).values(), key=len) == [{4, 5}, {1, 2, 3}]
        assert molcls.find(1) == id_of(molcls, {1, 2, 3})
        assert molcls.find(6) is False

    def test_unchanged_cluster_keeps_its_id(self, analyze: Analyze):
        molcls = analyze([[[1, 2, 3]], [[1, 2, 3]]], 4)
        cid = id_of(molcls, {1, 2, 3})

        step(molcls, 1)

        assert snapshot(molcls) == {cid: {1, 2, 3}}

    def test_growing_cluster_keeps_its_id(self, analyze: Analyze):
        molcls = analyze([[[1, 2, 3]], [[1, 2, 3, 4]]], 5)
        cid = id_of(molcls, {1, 2, 3})

        step(molcls, 1)

        assert snapshot(molcls) == {cid: {1, 2, 3, 4}}
        assert molcls.find(4) == cid

    def test_shrinking_cluster_keeps_its_id(self, analyze: Analyze):
        molcls = analyze([[[1, 2, 3, 4]], [[1, 2, 3]]], 4)
        cid = id_of(molcls, {1, 2, 3, 4})

        step(molcls, 1)

        assert snapshot(molcls) == {cid: {1, 2, 3}}
        assert molcls.find(4) is False

    def test_formation_creates_a_new_cluster(self, analyze: Analyze):
        molcls = analyze([[], [[1, 2]]], 3)
        assert molcls.clusters == {}

        step(molcls, 1)

        assert list(snapshot(molcls).values()) == [{1, 2}]

    def test_dissolution_removes_the_cluster(self, analyze: Analyze):
        molcls = analyze([[[1, 2]], []], 3)

        step(molcls, 1)

        assert molcls.clusters == {}
        assert molcls.mol_clt == {}

    def test_split_larger_fragment_keeps_the_id(self, analyze: Analyze):
        molcls = analyze([[[1, 2, 3, 4, 5]], [[1, 2, 3], [4, 5]]], 5)
        cid = id_of(molcls, {1, 2, 3, 4, 5})

        step(molcls, 1)

        assert id_of(molcls, {1, 2, 3}) == cid
        assert id_of(molcls, {4, 5}) > cid

    def test_merge_larger_cluster_keeps_its_id_even_if_younger(self, analyze: Analyze):
        molcls = analyze([[[1, 2]], [[1, 2], [3, 4, 5]], [[1, 2, 3, 4, 5]]], 5)
        step(molcls, 1)
        small, big = id_of(molcls, {1, 2}), id_of(molcls, {3, 4, 5})
        assert small < big

        step(molcls, 2)

        assert snapshot(molcls) == {big: {1, 2, 3, 4, 5}}

    def test_merge_of_equal_sizes_keeps_the_older_id(self, analyze: Analyze):
        molcls = analyze([[[4, 5, 6]], [[4, 5, 6], [1, 2, 3]], [[1, 2, 3, 4, 5, 6]]], 6)
        step(molcls, 1)
        older, younger = id_of(molcls, {4, 5, 6}), id_of(molcls, {1, 2, 3})
        assert older < younger

        step(molcls, 2)

        assert snapshot(molcls) == {older: {1, 2, 3, 4, 5, 6}}

    def test_three_way_merge_of_equal_sizes_keeps_the_oldest_id(self, analyze: Analyze):
        # the oldest cluster holds the highest resids, so it is not the first
        # candidate seen when the merged component is scanned
        molcls = analyze(
            [
                [[7, 8]],
                [[7, 8], [1, 2]],
                [[7, 8], [1, 2], [4, 5]],
                [[1, 2, 4, 5, 7, 8]],
            ],
            8,
        )
        step(molcls, 1)
        step(molcls, 2)
        oldest = id_of(molcls, {7, 8})
        middle, youngest = id_of(molcls, {1, 2}), id_of(molcls, {4, 5})
        assert oldest < middle < youngest

        step(molcls, 3)

        assert snapshot(molcls) == {oldest: {1, 2, 4, 5, 7, 8}}

    @pytest.mark.parametrize("order", [[6, 7, 5], [7, 5, 6], [5, 6, 7]])
    def test_older_cluster_is_the_lowest_id_among_all_top_ties(
        self, analyze: Analyze, order: list[int]
    ):
        molcls = analyze([[[1, 2]]], 2)
        origin = Counter({cid: 3 for cid in order})

        older = molcls._MolClusters__get_older_cluster({True: order, False: []}, origin)

        assert older == 5

    def test_older_cluster_ignores_lower_ids_with_smaller_counts(
        self, analyze: Analyze
    ):
        molcls = analyze([[[1, 2]]], 2)
        origin = Counter({6: 4, 7: 4, 5: 2})

        older = molcls._MolClusters__get_older_cluster(
            {True: [6, 7, 5], False: []}, origin
        )

        assert older == 6

    def test_dimer_bridging_two_clusters_is_always_new(self, analyze: Analyze):
        molcls = analyze([[[1, 2], [3, 4]], [[2, 3]]], 4)
        old_ids = set(molcls.clusters)

        step(molcls, 1)

        (new_id,) = molcls.clusters
        assert new_id not in old_ids
        assert snapshot(molcls)[new_id] == {2, 3}

    def test_one_molecule_from_each_cluster_is_new(self, analyze: Analyze):
        molcls = analyze([[[1, 2], [3, 4], [5, 6]], [[1, 3, 5]]], 6)
        old_ids = set(molcls.clusters)

        step(molcls, 1)

        (new_id,) = molcls.clusters
        assert new_id not in old_ids

    def test_fragment_loses_to_a_later_fragment_holding_more_of_the_cluster(
        self, analyze: Analyze
    ):
        # {1,2,10..14} is processed first (larger), but {3,4,5,6} carries more of
        # the original cluster, so the id must go to the latter.
        molcls = analyze(
            [[[1, 2, 3, 4, 5, 6]], [[1, 2, 10, 11, 12, 13, 14], [3, 4, 5, 6]]], 14
        )
        cid = id_of(molcls, {1, 2, 3, 4, 5, 6})

        step(molcls, 1)

        assert id_of(molcls, {3, 4, 5, 6}) == cid
        assert id_of(molcls, {1, 2, 10, 11, 12, 13, 14}) > cid

    def test_fragment_loses_to_a_later_dimer_fragment(self, analyze: Analyze):
        molcls = analyze([[[1, 2, 3, 4, 5]], [[1, 2, 10, 11, 12, 13], [3, 4]]], 13)
        cid = id_of(molcls, {1, 2, 3, 4, 5})

        step(molcls, 1)

        assert id_of(molcls, {3, 4}) == cid
        assert id_of(molcls, {1, 2, 10, 11, 12, 13}) > cid

    def test_unrelated_dimer_does_not_take_the_id(self, analyze: Analyze):
        # {20, 21} pair up elsewhere in the box, sharing nothing with the cluster
        molcls = analyze([[[1, 2, 3, 4, 5, 6]], [[1, 2, 3, 4, 5, 10], [20, 21]]], 21)
        cid = id_of(molcls, {1, 2, 3, 4, 5, 6})

        step(molcls, 1)

        assert id_of(molcls, {1, 2, 3, 4, 5, 10}) == cid
        assert id_of(molcls, {20, 21}) > cid

    def test_larger_fragment_keeps_the_id_over_a_dimer_fragment(self, analyze: Analyze):
        # the cluster splits into seven members (plus a newcomer) and a dimer; the
        # dimer holds less of the cluster, so it must not inherit the id
        molcls = analyze(
            [[[1, 2, 3, 4, 5, 6, 7, 8, 9]], [[1, 2, 3, 4, 5, 6, 7, 20], [8, 9]]], 20
        )
        cid = id_of(molcls, set(range(1, 10)))

        step(molcls, 1)

        assert id_of(molcls, {1, 2, 3, 4, 5, 6, 7, 20}) == cid
        assert id_of(molcls, {8, 9}) > cid

    def test_surviving_fragment_with_new_members_keeps_the_id(self, analyze: Analyze):
        # half the cluster dissolves while the other half picks up a newcomer
        molcls = analyze([[[1, 2, 3, 4, 5, 6]], [[1, 2, 3, 10]]], 10)
        cid = id_of(molcls, {1, 2, 3, 4, 5, 6})

        step(molcls, 1)

        assert snapshot(molcls) == {cid: {1, 2, 3, 10}}

    def test_cluster_already_continued_cannot_claim_a_second_fragment(
        self, analyze: Analyze
    ):
        molcls = analyze(
            [[[1, 2, 3, 4, 5, 6], [20, 21]], [[1, 2, 3, 4], [5, 6, 20]]], 21
        )
        big = id_of(molcls, {1, 2, 3, 4, 5, 6})
        small = id_of(molcls, {20, 21})

        step(molcls, 1)

        assert id_of(molcls, {1, 2, 3, 4}) == big
        assert id_of(molcls, {5, 6, 20}) not in {big, small}
        assert small not in molcls.clusters

    def test_ignored_composition_never_becomes_a_cluster(self, analyze: Analyze):
        resnames = ["MOL"] * 3 + ["SOL"] * 2
        molcls = analyze(
            [[[1, 2, 3], [4, 5]], [[1, 2, 3], [4, 5]]],
            5,
            resnames,
            rules={"MOL": {"MOL": f"cm {CUTOFF}"}, "SOL": {"SOL": f"cm {CUTOFF}"}},
            ignore_composition=[["SOL"]],
        )
        assert list(snapshot(molcls).values()) == [{1, 2, 3}]

        step(molcls, 1)

        assert list(snapshot(molcls).values()) == [{1, 2, 3}]
        assert molcls.find(4) is False

    @pytest.mark.xfail(
        strict=True,
        reason="dominance only scans later (smaller) groups and a merge loser is not "
        "marked as used, so it keeps its id in a remnant only if the remnant's group "
        "is processed after the merge group",
    )
    def test_merge_loser_fate_does_not_depend_on_group_order(self, analyze: Analyze):
        # A loses a merge to B in both cases, leaving a remnant elsewhere; only the
        # remnant's size relative to the merge group differs (free newcomers)
        a, b = list(range(1, 11)), list(range(11, 21))
        remnant_smaller = analyze(
            [[a, b], [list(range(1, 9)) + list(range(11, 20)), [9, 10, 30]]], 30
        )
        a_smaller = id_of(remnant_smaller, set(a))

        a, b = list(range(1, 11)), list(range(11, 18))
        remnant_larger = analyze(
            [[a, b], [[1, 2, 3, 4, *range(30, 42)], list(range(5, 11)) + b]], 41
        )
        a_larger = id_of(remnant_larger, set(a))

        step(remnant_smaller, 1)
        step(remnant_larger, 1)

        survives = {
            "remnant smaller": a_smaller in remnant_smaller.clusters,
            "remnant larger": a_larger in remnant_larger.clusters,
        }
        assert len(set(survives.values())) == 1, f"A survives: {survives}"

    @pytest.mark.parametrize(
        ("frames", "n_res", "pure"),
        [
            pytest.param(
                [[[1, 2, 3, 4]], [[1, 2, 10], [3, 4]]], 10, {3, 4}, id="dimer"
            ),
            pytest.param(
                [[[1, 2, 3, 4, 5, 6]], [[1, 2, 3, 10, 11], [4, 5, 6]]],
                11,
                {4, 5, 6},
                id="trimer",
                marks=pytest.mark.xfail(
                    strict=True,
                    reason="the pure-fragment tie-break in __check_dominance only "
                    "applies to dimers (len(subconn) == 2)",
                ),
            ),
        ],
    )
    def test_pure_fragment_wins_a_tie_whatever_its_size(
        self, analyze: Analyze, frames: list[Groups], n_res: int, pure: set[int]
    ):
        # the cluster splits into two halves; one half picked up newcomers, the
        # other is made only of the cluster's own molecules and must keep the id
        molcls = analyze(frames, n_res)
        cid = id_of(molcls, set(frames[0][0]))

        step(molcls, 1)

        assert id_of(molcls, pure) == cid

    @pytest.mark.parametrize(
        "halves",
        [
            [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]],
            [[6, 7, 8, 9, 10], [1, 2, 3, 4, 5]],
            [[2, 4, 6, 8, 10], [1, 3, 5, 7, 9]],
        ],
        ids=["low-first", "high-first", "interleaved"],
    )
    def test_even_split_keeps_the_id_in_the_half_with_the_lowest_resid(
        self, analyze: Analyze, halves: list[list[int]]
    ):
        # pins the current tie-break: equal-size groups keep networkx's component
        # order (resid order), so the half holding resid 1 wins. Arbitrary rather
        # than physical; update if a spatial tie-break is adopted.
        molcls = analyze([[list(range(1, 11))], halves], 10)
        cid = id_of(molcls, set(range(1, 11)))

        step(molcls, 1)

        (low,) = [set(h) for h in halves if 1 in h]
        (high,) = [set(h) for h in halves if 1 not in h]
        assert id_of(molcls, low) == cid
        assert id_of(molcls, high) > cid

    def test_one_frame_break_mints_a_transient_id(self, analyze: Analyze):
        # pins the current behaviour: no persistence window, so a contact broken for
        # a single frame creates an id that dies when the halves rejoin
        molcls = analyze(
            [[[1, 2, 3, 4, 5, 6]], [[1, 2, 3], [4, 5, 6]], [[1, 2, 3, 4, 5, 6]]], 6
        )
        cid = id_of(molcls, {1, 2, 3, 4, 5, 6})

        step(molcls, 1)
        transient = id_of(molcls, {4, 5, 6})
        assert transient > cid

        step(molcls, 2)

        assert snapshot(molcls) == {cid: {1, 2, 3, 4, 5, 6}}

    def test_ignored_composition_group_does_not_claim_the_id(self, analyze: Analyze):
        # the cluster's solvent shell {3, 4, 5} drifts off as a pure-solvent group,
        # which is ignored; it holds more of the cluster than {1, 2, 6}, but ignored
        # groups are invisible to the dominance check, so the id stays with the solute
        resnames = ["MOL"] + ["SOL"] * 4 + ["MOL"]
        molcls = analyze(
            [[[1, 2, 3, 4, 5]], [[1, 2, 6], [3, 4, 5]]],
            6,
            resnames,
            rules=ALL_PAIRS_RULES,
            ignore_composition=[["SOL"]],
        )
        cid = id_of(molcls, {1, 2, 3, 4, 5})

        step(molcls, 1)

        assert snapshot(molcls) == {cid: {1, 2, 6}}


# resids 1-3 and 7-8 are solute (MOL), 4-6 and 9-10 solvent (SOL); {9, 10} is a
# pure-solvent cluster throughout, which every solute analysis must ignore
RUN_RESNAMES = ["MOL"] * 3 + ["SOL"] * 3 + ["MOL"] * 2 + ["SOL"] * 2
RUN_FRAMES = [
    [[1, 4, 5], [7, 8], [9, 10]],  # one solvated solute, one bare two-solute cluster
    [[1, 4, 5, 6], [7, 8], [9, 10]],
    [[7, 8], [9, 10]],  # no solute-solvent cluster at all
]
ALL_PAIRS_RULES = {
    "MOL": {"MOL": f"cm {CUTOFF}", "SOL": f"cm {CUTOFF}"},
    "SOL": {"SOL": f"cm {CUTOFF}"},
}


class TestRun:
    def test_minimal_run_writes_evolution_and_json(
        self, analyze: Analyze, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        molcls = analyze([[[1, 2, 3]], []], 3)
        monkeypatch.chdir(tmp_path)

        molcls.run()

        evo = np.loadtxt(tmp_path / "evo.txt")
        np.testing.assert_allclose(evo, [[0, 1, 3, 3, 3], [1, 0, 0, 0, 0]])

        data = json.loads((tmp_path / "molclusters.json").read_text())
        assert data["Config"]["rules"] == MOL_RULES
        assert [f["NClusters"] for f in data["MolClusters"]] == [1, 0]
        (cls,) = data["MolClusters"][0]["Clusters"]
        assert cls["Size"] == 3
        assert cls["ResIDs"] == [1, 2, 3]
        assert cls["Composition"] == [{"resname": "MOL", "n": 3, "resids": [1, 2, 3]}]
        assert {(a, b) for a, b, _ in cls["Connections"]} == {(1, 2), (2, 3)}

        assert not (tmp_path / "solute_solvent.csv").exists()
        assert not (tmp_path / "nucleus_data.csv").exists()

    @pytest.fixture
    def full_run(
        self, analyze: Analyze, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> Path:
        molcls = analyze(
            RUN_FRAMES,
            10,
            RUN_RESNAMES,
            rules=ALL_PAIRS_RULES,
            solute=["MOL"],
            solvent=["SOL"],
            nucleus=["MOL"],
            follow=["solute"],
        )
        monkeypatch.chdir(tmp_path)
        molcls.run()
        return tmp_path

    def test_solute_solvent_table(self, full_run: Path):
        df = pd.read_csv(full_run / "solute_solvent.csv")

        assert list(df.columns) == [
            "Time", "NCls", "NSolt", "NSolv", "Radius",
            "Density", "Charge", "Dipole", "Spher", "Shape",
        ]  # fmt: skip
        assert df["Time"].tolist() == [0, 1, 2]
        assert df["NCls"].tolist() == [1, 1, 0]
        assert df["NSolt"].tolist() == [1, 1, 0]
        assert df["NSolv"].tolist() == [2, 3, 0]
        assert df.loc[:1, "Radius"].gt(0).all()
        assert np.isnan(df.loc[2, "Radius"])

    def test_nucleus_table(self, full_run: Path):
        df = pd.read_csv(full_run / "nucleus_data.csv")

        # frame 0: nuclei {1} (inside {1,4,5}) and {7,8}; frame 2: only {7,8}
        assert df["NNuc"].tolist() == [1, 1, 1]
        assert df["Size"].tolist() == [1.5, 1.5, 2]

    def test_json_records_nuclei(self, full_run: Path):
        data = json.loads((full_run / "molclusters.json").read_text())

        frame0 = {tuple(c["ResIDs"]): c for c in data["MolClusters"][0]["Clusters"]}
        assert [n["ResIDs"] for n in frame0[(1, 4, 5)]["Nucleus"]] == [[1]]
        assert [n["ResIDs"] for n in frame0[(7, 8)]["Nucleus"]] == [[7, 8]]
        assert "NucleiDipole" in frame0[(7, 8)]

    def test_coordinates_are_written_per_size_id_and_followed_solute(
        self, full_run: Path
    ):
        # coordinates are only written for frames after the first
        assert {p.name for p in full_run.glob("cls-n*.gro")} == {
            "cls-n2.gro",
            "cls-n4.gro",
        }
        assert len(list(full_run.glob("cls-id*.gro"))) == 2

        # {1,4,5,6} has exactly one solute to follow; {7,8} has two, so it's skipped
        assert [p.name for p in full_run.glob("solute-*.gro")] == ["solute-1.gro"]
        frames = (full_run / "solute-1.gro").read_text().count("Cluster-")
        assert frames == 1

    def test_nucleus_analysis_does_not_distort_the_cluster_geometry(
        self,
        analyze: Analyze,
        make_universe: UniverseFactory,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ):
        # the nucleus (MOL 3, 4) is not the cluster's first residue, so centering
        # the nucleus in place used to shift it away from the rest of the cluster
        frames = [[[1, 2, 3, 4]], [[1, 2, 3, 4]]]
        resnames = ["SOL", "SOL", "MOL", "MOL"]
        molcls = analyze(frames, 4, resnames, rules=ALL_PAIRS_RULES, nucleus=["MOL"])
        pristine = MDAResidueGroupAnalyzer(
            make_universe(frames, 4, resnames), [1, 2, 3, 4]
        )
        monkeypatch.chdir(tmp_path)

        molcls.run()

        data = json.loads((tmp_path / "molclusters.json").read_text())
        for frame in data["MolClusters"]:
            (cluster,) = frame["Clusters"]
            assert cluster["Radius"] == pytest.approx(pristine.radius_of_gyration)
            assert cluster["Shape"] == pytest.approx(pristine.shape_parameter)


class TestWriteCoordinates:
    def test_writes_both_size_and_id_grouped_files(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        uni = Universe(str(DATA_DIR / "met-mal.tpr"), str(DATA_DIR / "start.pdb"))
        config = MolClsConfig(rules={"MOL": {"MOL": "cm 15.0"}}, solute=["MOL"])
        molcls = MolClusters(uni, config)

        solute_clusters = [
            cls
            for cls in molcls.clusters.values()
            if set(config.solute).intersection(cls.resnames)
        ]
        assert solute_clusters, "fixture/rule setup should yield a solute cluster"

        monkeypatch.chdir(tmp_path)
        molcls._MolClusters__write_coordinates()

        size_files = list(tmp_path.glob("cls-n*.gro"))
        id_files = list(tmp_path.glob("cls-id*.gro"))

        assert size_files, "expected a size-grouped .gro output"
        assert len(id_files) == len(solute_clusters), (
            "expected one id-grouped .gro output per solute cluster"
        )
