# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

import json
import time
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from loguru import logger

from molclusters.analysis import (
    ClusterCoordinates,
    Frame,
    FrameAnalysis,
    JsonReport,
    Nucleus,
    Run,
    SizeEvolution,
    SoluteSolvent,
)
from molclusters.config import MolClsConfig
from molclusters.output import OutputFile, RunOutput
from molclusters.tracker import ClusterTracker

from .conftest import ALL_PAIRS_RULES, MOL_RULES, Groups, UniverseFactory


def run_analyses(
    analyses: Sequence[FrameAnalysis],
    make_universe: UniverseFactory,
    frames: Sequence[Groups],
    n_res: int,
    directory: Path,
    resnames: Sequence[str] | None = None,
    rules: dict = MOL_RULES,
) -> Run:
    """Drive `analyses` over every frame and finish them, as `MolClusters.run()` does.

    Returns
    -------
    Run
        The run's context.
    """
    uni = make_universe(frames, n_res, resnames)
    config = MolClsConfig(rules=rules)
    tracker = ClusterTracker(uni, config)
    run = Run(uni, config, len(frames), RunOutput(directory), analyses)

    with run.output:
        run.prepare_analyses()
        for i, _ in enumerate(uni.trajectory):
            if i:
                tracker.update()
            run.analyse_frame(Frame(i, tracker, run.output))
        run.finish_analyses()
    return run


class Recorder(FrameAnalysis):
    """Records what it sees of every frame, and which analyses it can find."""

    def __init__(self) -> None:
        self.seen: list[tuple[int, float, dict[int, set[int]]]] = []
        self.found: dict[str, FrameAnalysis | None] = {}

    def prepare(self, run: Run) -> None:
        self.found = {
            "size": run.analysis(SizeEvolution),
            "recorder": run.analysis(Recorder),
        }

    def analyse(self, frame: Frame) -> None:
        clusters = {cid: {int(m) for m in cls} for cid, cls in frame.clusters.items()}
        self.seen.append((frame.index, frame.time, clusters))


class TestFrameAnalysis:
    def test_only_analyse_must_be_implemented(self):
        class Incomplete(FrameAnalysis):
            pass

        with pytest.raises(TypeError, match="analyse"):
            Incomplete()

    def test_prepare_finish_and_outputs_default_to_nothing(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        class Minimal(FrameAnalysis):
            def analyse(self, frame: Frame) -> None:
                pass

        run_analyses([Minimal()], make_universe, [[[1, 2]], [[1, 2]]], 2, tmp_path)

        assert Minimal.outputs == ()
        assert list(tmp_path.iterdir()) == []


class TestRunAndFrame:
    def test_every_frame_is_seen_with_its_clusters(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        recorder = Recorder()
        frames = [[[1, 2, 3]], [], [[1, 2], [4, 5]]]

        run_analyses([recorder], make_universe, frames, 5, tmp_path)

        assert [(i, t) for i, t, _ in recorder.seen] == [(0, 0), (1, 1), (2, 2)]
        clusters = [c for _, _, c in recorder.seen]
        assert [sorted(map(sorted, c.values())) for c in clusters] == [
            [[1, 2, 3]],
            [],
            [[1, 2], [4, 5]],
        ]

    def test_frame_finds_the_cluster_of_a_molecule(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        uni = make_universe([[[1, 2, 3]]], 4)
        tracker = ClusterTracker(uni, MolClsConfig(rules=MOL_RULES))
        frame = Frame(0, tracker, RunOutput(tmp_path))

        (cid,) = frame.clusters
        assert frame.find(2) == cid
        assert frame.find(4) is None

    def test_frame_clusters_are_read_only(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        uni = make_universe([[[1, 2, 3]]], 4)
        tracker = ClusterTracker(uni, MolClsConfig(rules=MOL_RULES))

        with pytest.raises(TypeError):
            Frame(0, tracker, RunOutput(tmp_path)).clusters[99] = None

    def test_frames_append_through_the_run_output(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        class Appender(FrameAnalysis):
            outputs = (OutputFile("frames.txt"),)

            def analyse(self, frame: Frame) -> None:
                frame.output.append("frames.txt", f"{frame.index}\n")

        run = run_analyses([Appender()], make_universe, [[], [], []], 2, tmp_path)

        assert (tmp_path / "frames.txt").read_text() == "0\n1\n2\n"
        assert run.output.directory == tmp_path

    def test_preparing_finds_only_the_analyses_before(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        before, after = Recorder(), Recorder()
        size = SizeEvolution()

        run = run_analyses(
            [before, size, after], make_universe, [[[1, 2]]], 2, tmp_path
        )

        assert before.found == {"size": None, "recorder": None}
        assert after.found == {"size": size, "recorder": before}
        # once prepared, every analysis can be found
        assert run.analysis(SizeEvolution) is size


class TestSizeEvolution:
    def test_records_time_count_and_sizes_per_frame(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        size = SizeEvolution()

        run_analyses(
            [size], make_universe, [[[1, 2, 3], [4, 5]], [], [[1, 2]]], 5, tmp_path
        )

        np.testing.assert_allclose(
            size.data,
            [[0, 2, 2, 2.5, 3], [1, 0, 0, 0, 0], [2, 1, 2, 2, 2]],
        )

    def test_finish_writes_evo_and_logs_a_summary(
        self,
        make_universe: UniverseFactory,
        tmp_path: Path,
        captured_logs: list[str],
    ):
        size = SizeEvolution()

        run_analyses(
            [size], make_universe, [[[1, 2, 3], [4, 5]], [[1, 2, 3, 4]]], 5, tmp_path
        )

        text = (tmp_path / "evo.txt").read_text()
        assert text.startswith("# Time NClusters MinSize AvgSize MaxSize\n")
        np.testing.assert_allclose(np.loadtxt(tmp_path / "evo.txt"), size.data)
        (summary,) = [m for m in captured_logs if m.startswith("Found ")]
        assert summary.startswith(
            "Found 1.5 cluster(s) per frame on average (1-2); the largest held 4 "
        )

    def test_declares_evo(self):
        assert [out.name for out in SizeEvolution.outputs] == ["evo.txt"]


class TestSoluteSolvent:
    # frame 0: {1, 2, 3} holds a MOL and two SOL, {4, 5} only MOL; frame 1: none
    FRAMES = [[[1, 2, 3], [4, 5]], []]
    RESNAMES = ["MOL", "SOL", "SOL", "MOL", "MOL"]

    def run(self, make_universe: UniverseFactory, directory: Path) -> SoluteSolvent:
        analysis = SoluteSolvent(["MOL"], ["SOL"])
        run_analyses(
            [analysis],
            make_universe,
            self.FRAMES,
            5,
            directory,
            self.RESNAMES,
            ALL_PAIRS_RULES,
        )
        return analysis

    def test_records_only_clusters_with_solutes_and_solvents(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        analysis = self.run(make_universe, tmp_path)

        first, second = analysis.data
        np.testing.assert_allclose(first[:4], [0, 1, 1, 2])
        assert np.isfinite(first[4:]).all()
        np.testing.assert_allclose(second[:4], [1, 0, 0, 0])
        assert np.isnan(second[4:]).all()

    def test_finish_writes_the_table(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        analysis = self.run(make_universe, tmp_path)

        table = pd.read_csv(tmp_path / "solute_solvent.csv")
        assert list(table.columns) == [
            "Time",
            "NCls",
            "NSolt",
            "NSolv",
            "Radius",
            "Density",
            "Charge",
            "Dipole",
            "Spher",
            "Shape",
        ]
        np.testing.assert_allclose(table.to_numpy(), analysis.data)

    def test_declares_the_table(self):
        assert [out.name for out in SoluteSolvent.outputs] == ["solute_solvent.csv"]


class TestClusterCoordinates:
    # {1, 2} holds one solute (1); {3, 4} holds two (3, 4), and gains a SOL in
    # frame 2
    FRAMES = [[[1, 2], [3, 4]], [[1, 2], [3, 4]], [[1, 2], [3, 4, 5]]]
    RESNAMES = ["MOL", "SOL", "MOL", "MOL", "SOL"]

    def run(
        self, make_universe: UniverseFactory, directory: Path, *, follow: bool
    ) -> tuple[ClusterCoordinates, int, int]:
        """Run the analysis.

        Returns
        -------
        tuple[ClusterCoordinates, int, int]
            The analysis, and the ids of the clusters of resids 1 and 3.
        """
        analysis, recorder = ClusterCoordinates(["MOL"], follow=follow), Recorder()
        run_analyses(
            [recorder, analysis],
            make_universe,
            self.FRAMES,
            5,
            directory,
            self.RESNAMES,
            ALL_PAIRS_RULES,
        )
        (_, _, clusters) = recorder.seen[-1]
        (one,) = [cid for cid, mols in clusters.items() if 1 in mols]
        (three,) = [cid for cid, mols in clusters.items() if 3 in mols]
        return analysis, one, three

    @staticmethod
    def frames_in(file: Path) -> list[str]:
        return [
            line
            for line in file.read_text().splitlines()
            if line.startswith("Cluster-")
        ]

    def test_writes_each_frame_by_size_and_by_id(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        _, one, three = self.run(make_universe, tmp_path, follow=False)
        tmp_path /= "coordinates"

        assert sorted(p.name for p in tmp_path.iterdir()) == [
            "cls-id1.gro",
            "cls-id2.gro",
            "cls-n2.gro",
            "cls-n3.gro",
        ]
        assert self.frames_in(tmp_path / f"cls-id{three}.gro") == [
            f"Cluster-{three} - Time = 0.0",
            f"Cluster-{three} - Time = 1.0",
            f"Cluster-{three} - Time = 2.0",
        ]
        assert sorted(self.frames_in(tmp_path / "cls-n2.gro")) == sorted(
            [
                f"Cluster-{one} - Time = 0.0",
                f"Cluster-{three} - Time = 0.0",
                f"Cluster-{one} - Time = 1.0",
                f"Cluster-{three} - Time = 1.0",
                f"Cluster-{one} - Time = 2.0",
            ]
        )
        assert self.frames_in(tmp_path / "cls-n3.gro") == [
            f"Cluster-{three} - Time = 2.0"
        ]

    def test_follows_clusters_with_one_solute_and_warns_about_the_others(
        self,
        make_universe: UniverseFactory,
        tmp_path: Path,
        captured_logs: list[str],
    ):
        analysis, one, three = self.run(make_universe, tmp_path, follow=True)
        tmp_path /= "coordinates"

        assert self.frames_in(tmp_path / "solute-1.gro") == [
            f"Cluster-{one} - Time = 0.0",
            f"Cluster-{one} - Time = 1.0",
            f"Cluster-{one} - Time = 2.0",
        ]
        assert not list(tmp_path.glob("solute-[34].gro"))
        assert analysis.follow_skipped == {three: 3}
        (warning,) = [m for m in captured_logs if "were not followed" in m]
        assert warning.startswith("Solutes were not followed in 3 frame(s) of 1 ")

    @pytest.mark.parametrize(
        ("follow", "names"),
        [
            (False, ["coordinates/cls-n<size>.gro", "coordinates/cls-id<id>.gro"]),
            (
                True,
                [
                    "coordinates/cls-n<size>.gro",
                    "coordinates/cls-id<id>.gro",
                    "coordinates/solute-<resid>.gro",
                ],
            ),
        ],
    )
    def test_declares_its_files(self, follow: bool, names: list[str]):
        outputs = ClusterCoordinates(["MOL"], follow=follow).outputs

        assert [out.name for out in outputs] == names
        assert all(out.is_pattern for out in outputs)

    @pytest.mark.parametrize(("folder", "where"), [("", "."), ("gro/all", "gro/all")])
    def test_writes_in_the_folder_given(
        self, make_universe: UniverseFactory, tmp_path: Path, folder: str, where: str
    ):
        analysis = ClusterCoordinates(["MOL"], folder=folder)
        run_analyses(
            [analysis],
            make_universe,
            self.FRAMES,
            5,
            tmp_path,
            self.RESNAMES,
            ALL_PAIRS_RULES,
        )

        assert (tmp_path / where / "cls-n3.gro").exists()
        assert analysis.outputs[0].name == str(
            Path(where, "cls-n<size>.gro").as_posix()
        )


class TestNucleus:
    # the chain 1-2-3-4 has a SOL at 3, so its MOL nuclei are {1, 2} and {4}
    FRAMES = [[[1, 2, 3, 4]], []]
    RESNAMES = ["MOL", "MOL", "SOL", "MOL"]

    def test_finds_the_nuclei_of_each_cluster(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        seen = []

        class Reader(FrameAnalysis):
            def prepare(self, run: Run) -> None:
                self.nucleus = run.analysis(Nucleus)

            def analyse(self, frame: Frame) -> None:
                seen.append(
                    {
                        cid: sorted(sorted(int(r) for r in n.resids) for n in nuclei)
                        for cid, nuclei in self.nucleus.nuclei.items()
                    }
                )

        nucleus = Nucleus(["MOL"])
        run_analyses(
            [nucleus, Reader()],
            make_universe,
            self.FRAMES,
            4,
            tmp_path,
            self.RESNAMES,
            ALL_PAIRS_RULES,
        )

        ((cid, nuclei),) = seen[0].items()
        assert nuclei == [[1, 2], [4]]
        # a later frame's nuclei replace the earlier ones
        assert seen[1] == {}
        assert nucleus.nuclei == {}

    def test_records_the_averages_and_writes_the_table(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        nucleus = Nucleus(["MOL"])

        run_analyses(
            [nucleus],
            make_universe,
            self.FRAMES,
            4,
            tmp_path,
            self.RESNAMES,
            ALL_PAIRS_RULES,
        )

        first, second = nucleus.data
        np.testing.assert_allclose(first[:3], [0, 2, 1.5])
        assert np.isfinite(first[3:]).all()
        np.testing.assert_allclose(second[:2], [1, 0])
        assert np.isnan(second[2:]).all()
        table = pd.read_csv(tmp_path / "nucleus_data.csv")
        assert list(table.columns) == [
            "Time",
            "NNuc",
            "Size",
            "Radius",
            "Density",
            "Charge",
            "Dipole",
            "Spher",
            "Shape",
        ]
        np.testing.assert_allclose(table.to_numpy(), nucleus.data)

    def test_declares_the_table(self):
        assert [out.name for out in Nucleus.outputs] == ["nucleus_data.csv"]


class TestJsonReport:
    FRAMES = [[[1, 2, 3, 4]], []]
    RESNAMES = ["MOL", "MOL", "SOL", "MOL"]

    def run(
        self,
        analyses: list[FrameAnalysis],
        make_universe: UniverseFactory,
        directory: Path,
    ) -> dict:
        run_analyses(
            analyses,
            make_universe,
            self.FRAMES,
            4,
            directory,
            self.RESNAMES,
            ALL_PAIRS_RULES,
        )
        return json.loads((directory / "molclusters.json").read_text())

    def test_writes_every_frame_and_cluster(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        report = self.run([JsonReport()], make_universe, tmp_path)

        assert report["Software"].startswith("MolClusters ")
        assert report["Config"]["rules"]
        first, second = report["MolClusters"]
        assert (first["Time"], first["Frame"], first["NClusters"]) == (0, 0, 1)
        (cluster,) = first["Clusters"]
        assert cluster["ResIDs"] == [1, 2, 3, 4]
        assert "Nucleus" not in cluster
        assert (second["Time"], second["NClusters"], second["Clusters"]) == (1, 0, [])

    def test_includes_the_nuclei_of_a_nucleus_analysis_before_it(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        report = self.run([Nucleus(["MOL"]), JsonReport()], make_universe, tmp_path)

        (cluster,) = report["MolClusters"][0]["Clusters"]
        assert sorted(n["ResIDs"] for n in cluster["Nucleus"]) == [[1, 2], [4]]
        assert "NucleiDipole" in cluster

    def test_declares_the_report(self):
        assert [out.name for out in JsonReport.outputs] == ["molclusters.json"]


class TestErrors:
    @pytest.mark.parametrize(
        ("hook", "note"),
        [
            ("prepare", "Raised by Failing.prepare()"),
            ("analyse", "Raised by Failing.analyse() on frame 0"),
            ("finish", "Raised by Failing.finish()"),
        ],
    )
    def test_an_error_notes_the_analysis_and_hook(
        self, make_universe: UniverseFactory, tmp_path: Path, hook: str, note: str
    ):
        class Failing(FrameAnalysis):
            def prepare(self, run: Run) -> None:
                if hook == "prepare":
                    raise RuntimeError(hook)

            def analyse(self, frame: Frame) -> None:
                if hook == "analyse":
                    raise RuntimeError(hook)

            def finish(self, run: Run) -> None:
                if hook == "finish":
                    raise RuntimeError(hook)

        with pytest.raises(RuntimeError, match=hook) as err_info:
            run_analyses([Failing()], make_universe, [[[1, 2]]], 2, tmp_path)

        assert err_info.value.__notes__ == [note]


class Chatty(FrameAnalysis):
    """Logs from every hook, and takes a while to finish."""

    def prepare(self, run: Run) -> None:
        logger.info("prepare")

    def analyse(self, frame: Frame) -> None:
        logger.info("analyse {}", frame.index)

    def finish(self, run: Run) -> None:
        time.sleep(0.02)
        logger.info("finish")


class TestLogging:
    def test_what_an_analysis_logs_is_tagged_with_its_name(
        self,
        make_universe: UniverseFactory,
        tmp_path: Path,
        captured_logs: list[str],
    ):
        tags: dict[str, str | None] = {}
        handler_id = logger.add(
            lambda msg: tags.setdefault(
                msg.record["message"], msg.record["extra"].get("analysis")
            )
        )
        try:
            run_analyses(
                [Chatty(), SizeEvolution()], make_universe, [[[1, 2]]] * 2, 2, tmp_path
            )
            logger.info("after the run")
        finally:
            logger.remove(handler_id)

        assert tags["prepare"] == tags["analyse 1"] == tags["finish"] == "Chatty"
        (summary,) = [m for m in tags if m.startswith("Found ")]
        assert tags[summary] == "SizeEvolution"
        assert tags["after the run"] is None

    def test_the_time_of_each_analysis_is_added_up(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        run = run_analyses(
            [SizeEvolution(), Chatty()], make_universe, [[[1, 2]]] * 3, 2, tmp_path
        )

        assert len(run.durations) == 2
        assert run.durations[0] > 0
        assert run.durations[1] >= 0.02  # its finish sleeps that long
