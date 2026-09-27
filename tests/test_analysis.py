# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pytest

from molclusters.analysis import FrameAnalysis, SizeEvolution
from molclusters.config import MolClsConfig
from molclusters.tracker import ClusterTracker

from .conftest import MOL_RULES, Groups, UniverseFactory


def run_analysis(
    analysis: FrameAnalysis,
    make_universe: UniverseFactory,
    frames: Sequence[Groups],
    n_res: int,
) -> None:
    """Drive `analysis` over every frame, as `MolClusters.run()` does (no finish)."""
    uni = make_universe(frames, n_res)
    tracker = ClusterTracker(uni, MolClsConfig(rules=MOL_RULES))

    analysis.prepare(tracker, len(frames))
    analysis.analyse(tracker, 0)
    for i, _ in enumerate(uni.trajectory[1:], start=1):
        tracker.update()
        analysis.analyse(tracker, i)


class TestSizeEvolution:
    def test_records_time_count_and_sizes_per_frame(
        self, make_universe: UniverseFactory
    ):
        size = SizeEvolution()

        run_analysis(size, make_universe, [[[1, 2, 3], [4, 5]], [], [[1, 2]]], 5)

        np.testing.assert_allclose(
            size.data,
            [[0, 2, 2, 2.5, 3], [1, 0, 0, 0, 0], [2, 1, 2, 2, 2]],
        )

    def test_finish_writes_evo_and_logs_a_summary(
        self,
        make_universe: UniverseFactory,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        captured_logs: list[str],
    ):
        size = SizeEvolution()
        run_analysis(size, make_universe, [[[1, 2, 3], [4, 5]], [[1, 2, 3, 4]]], 5)
        monkeypatch.chdir(tmp_path)

        written = size.finish()

        assert written == ["evo.txt"]
        text = (tmp_path / "evo.txt").read_text()
        assert text.startswith("# Time NClusters MinSize AvgSize MaxSize\n")
        np.testing.assert_allclose(np.loadtxt(tmp_path / "evo.txt"), size.data)
        (summary,) = [m for m in captured_logs if m.startswith("Found ")]
        assert summary.startswith(
            "Found 1.5 cluster(s) per frame on average (1-2); the largest held 4 "
        )
