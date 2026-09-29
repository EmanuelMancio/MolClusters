# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

import re
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from molclusters.analysis import (
    ClusterCoordinates,
    Frame,
    FrameAnalysis,
    JsonReport,
    Lineage,
    Nucleus,
    Run,
    SizeEvolution,
    SoluteSolvent,
)
from molclusters.cluster import MolGroup
from molclusters.config import MolClsConfig
from molclusters.molclusters import MolClusters
from molclusters.output import OutputFile
from molclusters.report import read_report

from .conftest import ALL_PAIRS_RULES, CUTOFF, MOL_RULES, Groups, UniverseFactory

DATA_DIR = (Path(__file__).parent / "data" / "met-mal").resolve()

type Analyze = Callable[..., MolClusters]


@pytest.fixture
def analyze(make_universe: UniverseFactory) -> Analyze:
    """Build a MolClusters over synthetic frames.

    Returns
    -------
    Analyze
        ``analyze(frames, n_res, resnames=None, analyses=(), **config)``.
    """

    def factory(
        frames: Sequence[Groups],
        n_res: int,
        resnames: Sequence[str] | None = None,
        analyses: Sequence[FrameAnalysis] = (),
        **config_kwargs: Any,  # noqa: ANN401
    ) -> MolClusters:
        uni = make_universe(frames, n_res, resnames)
        config_kwargs.setdefault("rules", MOL_RULES)
        return MolClusters(uni, MolClsConfig(**config_kwargs), analyses)

    return factory


class TestConfigResolution:
    def test_effective_configuration_is_logged(
        self, analyze: Analyze, captured_logs: list[str]
    ):
        molcls = analyze([[]], 2, ["MOL", "SOL"], rules=ALL_PAIRS_RULES, solute=["MOL"])

        (logged,) = [m for m in captured_logs if m.startswith("Effective config")]
        assert molcls.config.describe() in logged
        assert "solvent: SOL" in logged

    def test_residues_not_numbered_one_to_n_are_refused_upfront(
        self, make_universe: UniverseFactory
    ):
        uni = make_universe([[[1, 2]]], 2)
        uni.residues.resids = [0, 1]

        with pytest.raises(ValueError, match="residue 1 has resid 0"):
            MolClusters(uni, MolClsConfig(rules=MOL_RULES))

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
    def test_minimal_run_writes_evolution_and_report(
        self, analyze: Analyze, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        molcls = analyze([[[1, 2, 3]], []], 3)
        monkeypatch.chdir(tmp_path)

        molcls.run()

        evo = np.loadtxt(tmp_path / "evo.txt")
        np.testing.assert_allclose(evo, [[0, 1, 3, 3, 3], [1, 0, 0, 0, 0]])

        report = read_report(tmp_path / "molclusters.jsonl.zst")
        assert report.header["config"]["rules"] == MOL_RULES
        frames = list(report.frames())
        assert [len(f["clusters"]) for f in frames] == [1, 0]
        (cls,) = frames[0]["clusters"]
        assert cls["size"] == 3
        assert cls["resids"] == [1, 2, 3]
        assert cls["composition"] == {"MOL": 3}
        assert {(i, j) for i, j, *_ in cls["connections"]} == {(1, 2), (2, 3)}

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

    # MDAnalysis' GRO writer writes a zero box then, as GROMACS does
    @pytest.mark.filterwarnings("ignore:missing dimension")
    def test_a_trajectory_without_a_box_runs_like_one_in_a_big_box(
        self,
        make_universe: UniverseFactory,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ):
        # the synthetic clusters sit far from the big box's faces, so the box
        # changes nothing but whether there is one, and where the groups are
        # made whole: at its center, ~1000 A, float32 positions keep ~1e-4 A
        config = MolClsConfig(
            rules=ALL_PAIRS_RULES,
            solute=["MOL"],
            solvent=["SOL"],
            nucleus=["MOL"],
            follow=["solute"],
        )
        boxed, vacuum = tmp_path / "boxed", tmp_path / "vacuum"
        for folder, box in ((boxed, True), (vacuum, False)):
            folder.mkdir()
            monkeypatch.chdir(folder)
            uni = make_universe(RUN_FRAMES, 10, RUN_RESNAMES, box=box)
            MolClusters(uni, config).run()

        tables = [
            read_report(folder / "molclusters.jsonl.zst").clusters()
            for folder in (vacuum, boxed)
        ]
        # records, not numbers: compared by the molecules they hold
        nuclei = [
            [[n["resids"] for n in row] for row in table.pop("Nucleus.nuclei")]
            for table in tables
        ]
        assert nuclei[0] == nuclei[1]
        pd.testing.assert_frame_equal(*tables, rtol=1e-4)
        for name in ("solute_solvent.csv", "nucleus_data.csv"):
            pd.testing.assert_frame_equal(
                pd.read_csv(vacuum / name), pd.read_csv(boxed / name), rtol=1e-4
            )
        assert (vacuum / "evo.txt").read_text() == (boxed / "evo.txt").read_text()

    def test_report_records_nuclei(self, full_run: Path):
        first = next(read_report(full_run / "molclusters.jsonl.zst").frames())

        frame0 = {tuple(c["resids"]): c["Nucleus"] for c in first["clusters"]}
        assert [n["resids"] for n in frame0[(1, 4, 5)]["nuclei"]] == [[1]]
        assert [n["resids"] for n in frame0[(7, 8)]["nuclei"]] == [[7, 8]]
        assert frame0[(7, 8)]["combined_dipole_moment"] is not None

    def test_coordinates_are_written_per_size_id_and_followed_solute(
        self, full_run: Path
    ):
        assert {p.name for p in full_run.glob("coordinates/cls-n*.gro")} == {
            "cls-n2.gro",
            "cls-n3.gro",
            "cls-n4.gro",
        }
        assert len(list(full_run.glob("coordinates/cls-id*.gro"))) == 2

        # {1,4,5(,6)} has exactly one solute to follow; {7,8} has two, so it's skipped
        assert [p.name for p in full_run.glob("coordinates/solute-*.gro")] == [
            "solute-1.gro"
        ]
        frames = (
            (full_run / "coordinates" / "solute-1.gro").read_text().count("Cluster-")
        )
        assert frames == 2

    def test_skipped_solute_following_is_summarized_once(
        self,
        analyze: Analyze,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        captured_logs: list[str],
    ):
        # {7,8} holds two solutes in every frame, so it can't be followed
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
        assert "in 3 frame(s) of 1 cluster(s)" in warning
        assert any("Results written to" in m for m in captured_logs)

    def test_results_go_to_the_output_dir(
        self,
        analyze: Analyze,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        captured_logs: list[str],
    ):
        molcls = analyze([[[1, 2]], [[1, 2]]], 2)
        monkeypatch.chdir(tmp_path)
        out = tmp_path / "results" / "run1"

        molcls.run(output_dir="results/run1")
        molcls.run(output_dir=out)

        assert sorted(p.name for p in out.iterdir()) == [
            "cluster_events.csv",
            "cluster_lifetimes.csv",
            "evo.txt",
            "molclusters.jsonl.zst",
        ]
        assert not (tmp_path / "evo.txt").exists()
        # the earlier-run check looks there too
        overwriting = [m for m in captured_logs if m.startswith("Overwriting")]
        assert len(overwriting) == 1
        (summary, _) = [m for m in captured_logs if m.startswith("Results written")]
        assert summary.startswith(f"Results written to {out}: ")

    def test_the_tracker_is_created_by_run(
        self, analyze: Analyze, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        molcls = analyze([[[1, 2]], [[1, 2, 3]]], 3)
        monkeypatch.chdir(tmp_path)
        assert molcls.tracker is None

        molcls.run()

        assert [sorted(c) for c in molcls.tracker.clusters.values()] == [[1, 2, 3]]

    def test_running_again_gives_the_same_results(
        self,
        analyze: Analyze,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        captured_logs: list[str],
    ):
        molcls = analyze(
            RUN_FRAMES, 10, RUN_RESNAMES, rules=ALL_PAIRS_RULES, nucleus=["MOL"]
        )
        monkeypatch.chdir(tmp_path)
        files = [
            "evo.txt",
            "cluster_events.csv",
            "cluster_lifetimes.csv",
            "molclusters.jsonl.zst",
            "nucleus_data.csv",
        ]

        molcls.run()
        first = {name: (tmp_path / name).read_bytes() for name in files}
        molcls.uni.trajectory[1]  # wherever the trajectory was left
        molcls.run()

        assert {name: (tmp_path / name).read_bytes() for name in files} == first
        connections = [m for m in captured_logs if m.startswith("Connections per")]
        assert len(connections) == 2
        assert connections[0] == connections[1]

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
        (durations,) = [m for m in captured_logs if m.startswith("Time spent")]
        # then the built-ins the config enables, in order
        assert re.fullmatch(
            r"Time spent tracking the clusters: [\d.]+s; in each analysis: "
            r"SizeEvolution [\d.]+s, Lineage [\d.]+s, JsonReport [\d.]+s\s*",
            durations,
        )

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
        (tmp_path / "coordinates").mkdir()
        for name in [
            "evo.txt",
            "coordinates/cls-n2.gro",
            "coordinates/cls-n5.gro",
            "coordinates/solute-1.gro",
            "coordinates/cls-n2-mine.gro",
            "cls-n3.gro",  # outside the folder
        ]:
            (tmp_path / name).write_text("earlier run\n")
        molcls = analyze([[[1, 2]], [[1, 2]]], 2, solute=["MOL"])

        molcls.run()

        # solute-1.gro is only written when following solutes
        assert (
            "Overwriting results of an earlier run: "
            "['evo.txt', 'coordinates/cls-n<size>.gro (2 file(s))']\n"
        ) in captured_logs
        # started over, not appended to
        gro = (tmp_path / "coordinates/cls-n2.gro").read_text()
        assert gro.startswith("Cluster-1 - Time = 0")
        # this run had no cluster of 5, so the earlier run's file is still there
        (warning,) = [m for m in captured_logs if "file(s) of an earlier run" in m]
        assert warning.startswith("1 file(s)")
        assert "(coordinates/cls-n5.gro)" in warning

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
        pristine = MolGroup(make_universe(frames, 4, resnames), [1, 2, 3, 4])
        monkeypatch.chdir(tmp_path)

        molcls.run()

        for frame in read_report(tmp_path / "molclusters.jsonl.zst").frames():
            (cluster,) = frame["clusters"]
            assert cluster["radius"] == pytest.approx(pristine.radius)
            assert cluster["shape_parameter"] == pytest.approx(pristine.shape_parameter)


class LargestCluster(FrameAnalysis):
    """A user analysis: the largest cluster of every frame, next to the built-ins."""

    outputs = (OutputFile("largest.txt"), OutputFile("largest-<n>.log"))

    def prepare(self, run: Run) -> None:
        self.size = run.analysis(SizeEvolution)
        self.largest: list[int] = []

    def analyse(self, frame: Frame) -> None:
        self.largest.append(max((c.size for c in frame.clusters.values()), default=0))

    def finish(self, run: Run) -> None:
        run.output.path("largest.txt").write_text(f"{self.largest}\n")
        run.output.append("largest-1.log", "done\n")


class TestUserAnalyses:
    def test_run_after_the_builtins_but_the_report_and_see_every_frame(
        self, analyze: Analyze, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        largest = LargestCluster()
        molcls = analyze([[[1, 2, 3]], [], [[1, 2], [3, 4]]], 4, analyses=[largest])
        monkeypatch.chdir(tmp_path)

        molcls.run()

        # the report comes last, for every analysis to add to it
        assert molcls.analyses[-2:] == [largest, molcls.analysis(JsonReport)]
        assert largest.largest == [3, 0, 2]
        # a built-in's results are found, and complete for the frames seen
        assert largest.size is molcls.analysis(SizeEvolution)
        assert (tmp_path / "largest.txt").read_text() == "[3, 0, 2]\n"
        assert (tmp_path / "largest-1.log").read_text() == "done\n"

    def test_their_outputs_are_checked_and_reported(
        self,
        analyze: Analyze,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        captured_logs: list[str],
    ):
        monkeypatch.chdir(tmp_path)
        for name in ["largest.txt", "largest-1.log", "largest-7.log"]:
            (tmp_path / name).write_text("earlier run\n")
        molcls = analyze([[[1, 2]], [[1, 2]]], 2, analyses=[LargestCluster()])

        molcls.run()

        overwriting = (
            "Overwriting results of an earlier run: "
            "['largest.txt', 'largest-<n>.log (2 file(s))']\n"
        )
        assert overwriting in captured_logs
        # started over, not appended to
        assert (tmp_path / "largest-1.log").read_text() == "done\n"
        (warning,) = [m for m in captured_logs if "file(s) of an earlier run" in m]
        assert warning.startswith("1 file(s)")
        assert "(largest-7.log)" in warning
        (summary,) = [m for m in captured_logs if m.startswith("Results written")]
        assert summary.endswith(
            ": evo.txt, cluster_events.csv, cluster_lifetimes.csv, largest.txt, "
            "largest-<n>.log, molclusters.jsonl.zst\n"
        )

    def test_none_are_added_by_default(self, analyze: Analyze):
        molcls = analyze([[[1, 2]]], 2)

        assert [type(a) for a in molcls.analyses] == [
            SizeEvolution,
            Lineage,
            JsonReport,
        ]


class TestBuiltins:
    def test_the_config_enables_them_in_order(self, analyze: Analyze):
        molcls = analyze(
            [[]],
            2,
            ["MOL", "SOL"],
            rules=ALL_PAIRS_RULES,
            solute=["MOL"],
            nucleus=["MOL"],
            follow=["solute"],
        )

        assert [type(a) for a in molcls.analyses] == [
            SizeEvolution,
            Lineage,
            SoluteSolvent,
            ClusterCoordinates,
            Nucleus,
            JsonReport,
        ]
        solute, coordinates = molcls.analyses[2:4]
        assert solute.solutes == {"MOL"}
        assert solute.solvents == {"SOL"}
        assert coordinates.follow is True
        assert molcls.analysis(Nucleus) is molcls.analyses[4]

    def test_the_config_turns_them_off(
        self,
        make_universe: UniverseFactory,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ):
        config = MolClsConfig(
            rules=ALL_PAIRS_RULES,
            solute=["MOL"],
            nucleus=["MOL"],
            analyses={"ClusterCoordinates": False, "Nucleus": False},
        )
        molcls = MolClusters(make_universe(RUN_FRAMES, 10, RUN_RESNAMES), config)
        monkeypatch.chdir(tmp_path)

        molcls.run()

        assert [type(a) for a in molcls.analyses] == [
            SizeEvolution,
            Lineage,
            SoluteSolvent,
            JsonReport,
        ]
        assert not (tmp_path / "coordinates").exists()
        assert not (tmp_path / "nucleus_data.csv").exists()
        # the report leaves the nuclei out, as when no nucleus is configured
        frames = read_report(tmp_path / "molclusters.jsonl.zst").frames()
        assert all(
            "Nucleus" not in cls for frame in frames for cls in frame["clusters"]
        )

    def test_turning_one_on_keeps_the_rest(self, make_universe: UniverseFactory):
        config = MolClsConfig(rules=MOL_RULES, analyses={"SizeEvolution": True})
        molcls = MolClusters(make_universe([[]], 2), config)

        assert [type(a) for a in molcls.analyses] == [
            SizeEvolution,
            Lineage,
            JsonReport,
        ]

    def test_every_one_turned_off_warns(
        self,
        make_universe: UniverseFactory,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        captured_logs: list[str],
    ):
        config = MolClsConfig(
            rules=MOL_RULES,
            analyses={"SizeEvolution": False, "Lineage": False, "JsonReport": False},
        )
        molcls = MolClusters(make_universe([[[1, 2]]], 2), config)
        monkeypatch.chdir(tmp_path)

        molcls.run()

        assert molcls.analyses == []
        (warning,) = [m for m in captured_logs if m.startswith("Every analysis")]
        assert "no results written" in warning
        assert list(tmp_path.iterdir()) == []

    def test_analysis_finds_one_by_type(self, analyze: Analyze):
        largest = LargestCluster()
        molcls = analyze([[]], 2, analyses=[largest])

        assert molcls.analysis(LargestCluster) is largest
        assert molcls.analysis(FrameAnalysis) is molcls.analyses[0]
        # not enabled by the config
        assert molcls.analysis(Nucleus) is None

    def test_a_class_instead_of_an_instance_is_refused(self, analyze: Analyze):
        with pytest.raises(TypeError, match="pass an instance, not the class"):
            analyze([[[1, 2]]], 2, analyses=[LargestCluster])

    def test_other_objects_are_refused(self, analyze: Analyze):
        with pytest.raises(TypeError, match="FrameAnalysis instances, got 3"):
            analyze([[[1, 2]]], 2, analyses=[3])

    def test_an_error_names_the_analysis_and_frame(
        self, analyze: Analyze, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        class Broken(FrameAnalysis):
            def analyse(self, frame: Frame) -> None:
                if frame.index == 1:
                    raise ValueError("boom")

        molcls = analyze([[[1, 2]], [[1, 2]]], 2, analyses=[Broken()])
        monkeypatch.chdir(tmp_path)

        with pytest.raises(ValueError, match="boom") as err_info:
            molcls.run()

        assert err_info.value.__notes__ == ["Raised by Broken.analyse() on frame 1"]
