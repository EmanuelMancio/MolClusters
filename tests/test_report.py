# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

import gzip
import json
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
import pytest

from molclusters.analysis import Frame, FrameAnalysis, JsonReport, Nucleus
from molclusters.cluster import MolGroup
from molclusters.config import MolClsConfig
from molclusters.molclusters import MolClusters
from molclusters.output import RunOutput
from molclusters.report import (
    DECIMALS,
    FORMAT,
    VERSION,
    connection_rows,
    dumps,
    group_record,
    read_report,
    report_name,
)
from molclusters.tracker import ClusterTracker

from .conftest import ALL_PAIRS_RULES, MOL_RULES, UniverseFactory
from .test_analysis import run_analyses

FRAMES = [[[1, 2, 3]], [[1, 2], [3, 4]], []]


class Interrupt(FrameAnalysis):
    """Interrupts the run on a frame, as Ctrl+C would."""

    def __init__(self, frame: int) -> None:
        self.frame = frame

    def analyse(self, frame: Frame) -> None:
        if frame.index == self.frame:
            raise KeyboardInterrupt


class TestEncoding:
    @pytest.mark.parametrize(
        ("compression", "name"),
        [
            ("zstd", "molclusters.jsonl.zst"),
            ("gzip", "molclusters.jsonl.gz"),
            ("none", "molclusters.jsonl"),
        ],
    )
    def test_the_name_says_the_compression(self, compression: str, name: str):
        assert report_name(compression) == name

    def test_a_record_is_one_line(self):
        record = {"a": np.float64(1.5), "b": np.arange(3), "c": [np.int64(2)]}

        assert dumps(record) == b'{"a":1.5,"b":[0,1,2],"c":[2]}\n'

    def test_nan_is_written_as_null(self):
        # a bare NaN isn't JSON: strict parsers reject the whole file
        line = dumps({"x": float("nan"), "y": np.float64("nan")})

        assert json.loads(line) == {"x": None, "y": None}

    def test_floats_are_rounded(self, make_universe: UniverseFactory):
        uni = make_universe([[[1, 2]]], 2)

        record = group_record(MolGroup(uni, [1, 2]))

        floats = [v for v in record.values() if isinstance(v, float)]
        assert len(floats) == 9
        assert all(v == round(v, DECIMALS) for v in floats)
        # as the property is, to the decimals kept
        assert record["radius"] == pytest.approx(MolGroup(uni, [1, 2]).radius, abs=1e-6)

    def test_a_group_counts_its_residue_names(self, make_universe: UniverseFactory):
        uni = make_universe([[[1, 2, 3]]], 3, ["SOL", "MOL", "SOL"])

        record = group_record(MolGroup(uni, [3, 1, 2]))

        assert record["resids"] == [1, 2, 3]
        assert record["composition"] == {"SOL": 2, "MOL": 1}

    def test_cm_connections_have_no_angle_or_count(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        run_analyses([JsonReport()], make_universe, [[[1, 2]]], 2, tmp_path)

        (frame,) = read_report(tmp_path / "molclusters.jsonl.zst").frames()
        ((i, j, distance, angle, n_hbonds),) = frame["clusters"][0]["connections"]
        assert (i, j, angle, n_hbonds) == (1, 2, None, None)
        assert isinstance(distance, float)

    def test_connection_rows_follow_the_graph(self, make_universe: UniverseFactory):
        uni = make_universe([[[1, 2, 3]]], 3)
        (cluster,) = ClusterTracker(
            uni, MolClsConfig(rules=MOL_RULES)
        ).clusters.values()

        rows = connection_rows(cluster)

        assert [row[:2] for row in rows] == [list(e) for e in cluster.graph.edges]


class TestHeader:
    def test_describes_the_run(self, make_universe: UniverseFactory, tmp_path: Path):
        run_analyses([JsonReport()], make_universe, FRAMES, 4, tmp_path)

        header = read_report(tmp_path / "molclusters.jsonl.zst").header
        assert (header["format"], header["version"]) == (FORMAT, VERSION)
        assert header["n_frames"] == 3
        assert header["decimals"] == DECIMALS
        assert header["trajectory"] == str(Path("synthetic.traj").absolute())
        assert header["topology"] == str(Path("synthetic.top").absolute())
        assert header["connection_columns"] == [
            "i", "j", "distance", "angle", "n_hbonds"
        ]  # fmt: skip
        assert header["config"]["report_compression"] == "zstd"

    def test_an_in_memory_universe_has_no_paths(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        uni = make_universe([[[1, 2]]], 2)
        merged = mda.Merge(uni.atoms)
        merged.dimensions = uni.dimensions  # Merge leaves the box out
        assert merged.filename is None

        MolClusters(merged, MolClsConfig(rules=MOL_RULES)).run(tmp_path)

        header = read_report(tmp_path / "molclusters.jsonl.zst").header
        assert (header["trajectory"], header["topology"]) == (None, None)


class TestReader:
    def test_reads_every_compression(self, tmp_path: Path):
        lines = [
            dumps({"format": FORMAT, "version": VERSION, "n_frames": 1}),
            dumps({"frame": 0, "time": 0.0, "clusters": []}),
        ]
        for compression in ["zstd", "gzip", "none"]:
            with RunOutput(tmp_path) as output:
                for line in lines:
                    output.append(report_name(compression), line)

            report = read_report(tmp_path / report_name(compression))
            assert list(report) == [{"frame": 0, "time": 0.0, "clusters": []}]

    def test_the_compression_is_told_from_the_content(self, tmp_path: Path):
        header = dumps({"format": FORMAT, "version": VERSION})
        misnamed = tmp_path / "report.jsonl"
        misnamed.write_bytes(gzip.compress(header))

        assert read_report(misnamed).header["format"] == FORMAT

    def test_the_old_json_is_not_a_report(self, tmp_path: Path):
        old = tmp_path / "molclusters.json"
        old.write_text(json.dumps({"Software": "MolClusters 0.6.0"}, indent=2))

        with pytest.raises(ValueError, match="single JSON document"):
            read_report(old)

    def test_a_newer_version_is_refused(self, tmp_path: Path):
        newer = tmp_path / "molclusters.jsonl"
        newer.write_bytes(dumps({"format": FORMAT, "version": VERSION + 1}))

        with pytest.raises(ValueError, match="update MolClusters"):
            read_report(newer)

    @pytest.mark.parametrize("compression", ["zstd", "gzip", "none"])
    def test_an_interrupted_run_leaves_the_frames_it_analysed(
        self,
        make_universe: UniverseFactory,
        tmp_path: Path,
        captured_logs: list[str],
        compression: str,
    ):
        analyses = [JsonReport(compression), Interrupt(frame=1)]
        with pytest.raises(KeyboardInterrupt):
            run_analyses(analyses, make_universe, FRAMES, 4, tmp_path)

        frames = list(read_report(tmp_path / report_name(compression)).frames())

        assert [f["frame"] for f in frames] == [0, 1]
        (warning,) = [m for m in captured_logs if "ends early" in m]
        assert "after 2 frame(s) of the 3 of the run" in warning

    @pytest.mark.parametrize("compression", ["zstd", "gzip", "none"])
    def test_a_killed_run_leaves_what_was_flushed(
        self, tmp_path: Path, captured_logs: list[str], compression: str
    ):
        name = report_name(compression)
        output = RunOutput(tmp_path)
        output.append(name, dumps({"format": FORMAT, "version": VERSION}))
        output.append(name, dumps({"frame": 0, "time": 0.0, "clusters": []}))
        # killed while writing the next frame: never closed, the line cut off
        output.append(name, dumps({"frame": 1, "time": 1.0, "clusters": []})[:-9])
        output.flush()

        frames = list(read_report(tmp_path / name).frames())

        assert [f["frame"] for f in frames] == [0]
        assert any("ends early, after 1 frame(s):" in m for m in captured_logs)

    def test_a_complete_report_reads_without_warnings(
        self, make_universe: UniverseFactory, tmp_path: Path, captured_logs: list[str]
    ):
        run_analyses([JsonReport()], make_universe, FRAMES, 4, tmp_path)

        assert len(list(read_report(tmp_path / "molclusters.jsonl.zst"))) == 3
        assert not [m for m in captured_logs if "ends early" in m]


class TestTables:
    @pytest.fixture
    def report_path(self, make_universe: UniverseFactory, tmp_path: Path) -> Path:
        run_analyses(
            [Nucleus(["MOL"]), JsonReport()],
            make_universe,
            FRAMES,
            4,
            tmp_path,
            ["MOL", "MOL", "SOL", "SOL"],
            ALL_PAIRS_RULES,
        )
        return tmp_path / "molclusters.jsonl.zst"

    def test_clusters_have_a_row_per_frame(self, report_path: Path):
        table = read_report(report_path).clusters()

        assert table[["frame", "id", "size"]].values.tolist() == [
            [0, 1, 3],
            [1, 1, 2],
            [1, 2, 2],
        ]
        # a residue name missing from a cluster is NaN
        assert table["composition.MOL"].fillna(0).tolist() == [2, 2, 0]
        assert table["composition.SOL"].fillna(0).tolist() == [1, 0, 2]
        assert table["resids"].tolist() == [[1, 2, 3], [1, 2], [3, 4]]
        assert "connections" not in table
        assert table["Nucleus.combined_dipole_moment"].notna().tolist() == [
            True,
            True,
            False,
        ]

    def test_connections_have_a_row_per_frame(self, report_path: Path):
        table = read_report(report_path).connections()

        assert list(table.columns) == [
            "frame", "time", "cluster", "i", "j", "distance", "angle", "n_hbonds"
        ]  # fmt: skip
        assert table[["frame", "cluster", "i", "j"]].values.tolist() == [
            [0, 1, 1, 2],
            [0, 1, 2, 3],
            [1, 1, 1, 2],
            [1, 2, 3, 4],
        ]
        assert table["n_hbonds"].dtype == pd.Int64Dtype()
        assert table["n_hbonds"].isna().all()
