# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pytest

from molclusters.analysis import Frame, FrameAnalysis, Run, SizeEvolution
from molclusters.config import MolClsConfig
from molclusters.output import RunOutput
from molclusters.tracker import ClusterTracker

from .conftest import MOL_RULES, Groups, UniverseFactory


def run_analyses(
    analyses: Sequence[FrameAnalysis],
    make_universe: UniverseFactory,
    frames: Sequence[Groups],
    n_res: int,
    directory: Path,
) -> Run:
    """Drive `analyses` over every frame and finish them, as `MolClusters.run()` does.

    Returns
    -------
    Run
        The run's context.
    """
    uni = make_universe(frames, n_res)
    config = MolClsConfig(rules=MOL_RULES)
    tracker = ClusterTracker(uni, config)
    run = Run(uni, config, len(frames), RunOutput(directory), analyses)

    with run.output:
        run.prepare_analyses()
        for i, _ in enumerate(uni.trajectory):
            if i:
                tracker.update()
            run.analyse_frame(Frame(i, tracker))
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
        self, make_universe: UniverseFactory
    ):
        uni = make_universe([[[1, 2, 3]]], 4)
        tracker = ClusterTracker(uni, MolClsConfig(rules=MOL_RULES))
        frame = Frame(0, tracker)

        (cid,) = frame.clusters
        assert frame.find(2) == cid
        assert frame.find(4) is False

    def test_frame_clusters_are_read_only(self, make_universe: UniverseFactory):
        uni = make_universe([[[1, 2, 3]]], 4)
        tracker = ClusterTracker(uni, MolClsConfig(rules=MOL_RULES))

        with pytest.raises(TypeError):
            Frame(0, tracker).clusters[99] = None

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
