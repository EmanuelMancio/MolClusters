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

from molclusters.cluster import MDAResidueGroupAnalyzer
from molclusters.config import MolClsConfig
from molclusters.molclusters import MolClusters
from molclusters.output import RunOutput

from .conftest import ALL_PAIRS_RULES, CUTOFF, MOL_RULES, Groups, UniverseFactory

DATA_DIR = (Path(__file__).parent / "data" / "met-mal").resolve()

type Analyze = Callable[..., MolClusters]


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


class TestConfigResolution:
    def test_solvent_defaults_to_every_non_solute_resname(
        self, analyze: Analyze, captured_logs: list[str]
    ):
        molcls = analyze(
            [[]], 3, ["MOL", "SOL", "ION"], rules=ALL_PAIRS_RULES, solute=["MOL"]
        )

        assert molcls.config.solvent == ["ION", "SOL"]
        assert any("No 'solvent' configured" in m for m in captured_logs)

    def test_effective_configuration_is_logged(
        self, analyze: Analyze, captured_logs: list[str]
    ):
        molcls = analyze([[]], 2, ["MOL", "SOL"], rules=ALL_PAIRS_RULES, solute=["MOL"])

        (logged,) = [m for m in captured_logs if m.startswith("Effective config")]
        assert molcls.config.describe() in logged
        assert "solvent: SOL" in logged

    def test_resnames_missing_from_the_topology_warn(
        self, analyze: Analyze, captured_logs: list[str]
    ):
        analyze([[]], 2, ["MOL", "SOL"], solute=["MOL"], solvent=["WAT"])

        (warning,) = [m for m in captured_logs if "not in the topology" in m]
        assert warning.startswith("'solvent' names residue(s) ['WAT']")


# resids 1-3 and 7-8 are solute (MOL), 4-6 and 9-10 solvent (SOL); {9, 10} is a
# pure-solvent cluster throughout, which every solute analysis must ignore
RUN_RESNAMES = ["MOL"] * 3 + ["SOL"] * 3 + ["MOL"] * 2 + ["SOL"] * 2
RUN_FRAMES = [
    [[1, 4, 5], [7, 8], [9, 10]],  # one solvated solute, one bare two-solute cluster
    [[1, 4, 5, 6], [7, 8], [9, 10]],
    [[7, 8], [9, 10]],  # no solute-solvent cluster at all
]


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

    def test_skipped_solute_following_is_summarized_once(
        self,
        analyze: Analyze,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        captured_logs: list[str],
    ):
        # {7,8} holds two solutes in frames 1 and 2, so it can't be followed
        molcls = analyze(
            RUN_FRAMES,
            10,
            RUN_RESNAMES,
            rules=ALL_PAIRS_RULES,
            solute=["MOL"],
            follow=["solute"],
        )
        monkeypatch.chdir(tmp_path)

        molcls.run()

        (warning,) = [m for m in captured_logs if m.startswith("Solutes were not")]
        assert "in 2 frame(s) of 1 cluster(s)" in warning
        assert any("Results written to" in m for m in captured_logs)

    def test_progress_and_duration_are_logged(
        self,
        analyze: Analyze,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        captured_logs: list[str],
    ):
        molcls = analyze([[[1, 2]]] * 20, 2)
        monkeypatch.chdir(tmp_path)

        molcls.run()

        progress = [m for m in captured_logs if m.startswith("Frame ")]
        # every 10% of the 20 frames, except the last one
        assert [m.split(" (")[0] for m in progress] == [
            f"Frame {n}/20" for n in range(2, 20, 2)
        ]
        assert any(m.startswith("Tracked 20 frame(s) in ") for m in captured_logs)

    def test_connections_are_summarized_by_rule(
        self,
        analyze: Analyze,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        captured_logs: list[str],
    ):
        # MOL-MOL connects 1-2 in both frames; nothing ever connects to SOL
        molcls = analyze(
            [[[1, 2]], [[1, 2]]], 3, ["MOL", "MOL", "SOL"], rules=ALL_PAIRS_RULES
        )
        monkeypatch.chdir(tmp_path)

        molcls.run()

        (summary,) = [m for m in captured_logs if m.startswith("Connections per")]
        assert f"MOL - MOL (cm {CUTOFF}) 1.0" in summary
        (warning,) = [m for m in captured_logs if m.startswith("Rule(s)")]
        assert f"MOL - SOL (cm {CUTOFF}), SOL - SOL (cm {CUTOFF}) never" in warning

    def test_outputs_of_an_earlier_run_are_reported(
        self,
        analyze: Analyze,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        captured_logs: list[str],
    ):
        monkeypatch.chdir(tmp_path)
        for name in ["evo.txt", "cls-n2.gro", "solute-1.gro", "cls-n2-mine.gro"]:
            (tmp_path / name).write_text("")
        molcls = analyze([[[1, 2]], [[1, 2]]], 2, solute=["MOL"])

        molcls.run()

        assert "Overwriting results of an earlier run: ['evo.txt']\n" in captured_logs
        # solute-1.gro is only written when following solutes
        (warning,) = [m for m in captured_logs if ".gro file(s) from an" in m]
        assert warning.startswith("1 .gro file(s)")
        assert "(cls-n2.gro)" in warning

    def test_a_fresh_directory_reports_no_earlier_outputs(
        self,
        analyze: Analyze,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        captured_logs: list[str],
    ):
        monkeypatch.chdir(tmp_path)
        molcls = analyze([[[1, 2]], [[1, 2]]], 2, solute=["MOL"])

        molcls.run()

        assert not any("earlier run" in m for m in captured_logs)

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
    def test_writes_both_size_and_id_grouped_files(self, tmp_path: Path):
        uni = Universe(str(DATA_DIR / "met-mal.tpr"), str(DATA_DIR / "start.pdb"))
        config = MolClsConfig(rules={"MOL": {"MOL": "cm 15.0"}}, solute=["MOL"])
        molcls = MolClusters(uni, config)

        solute_clusters = [
            cls
            for cls in molcls.clusters.values()
            if set(config.solute).intersection(cls.resnames)
        ]
        assert solute_clusters, "fixture/rule setup should yield a solute cluster"

        molcls.output = RunOutput(tmp_path)  # as run() does
        molcls._MolClusters__start_solute_solvent()
        molcls._MolClusters__write_coordinates()
        molcls.output.flush()  # as run() does when it's done

        size_files = list(tmp_path.glob("cls-n*.gro"))
        id_files = list(tmp_path.glob("cls-id*.gro"))

        assert size_files, "expected a size-grouped .gro output"
        assert len(id_files) == len(solute_clusters), (
            "expected one id-grouped .gro output per solute cluster"
        )
