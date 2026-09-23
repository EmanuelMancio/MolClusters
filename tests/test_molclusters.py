# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

import json
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from MDAnalysis import Universe

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
