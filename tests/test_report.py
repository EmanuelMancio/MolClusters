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
from molclusters.cluster import Cluster, MolGroup
from molclusters.config import MolClsConfig
from molclusters.molclusters import MolClusters
from molclusters.output import RunOutput
from molclusters.report import (
    DECIMALS,
    FORMAT,
    VERSION,
    Report,
    cluster_record,
    connection_rows,
    dumps,
    group_record,
    read_report,
    report_name,
    rounded,
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
        assert len(floats) == 10
        assert all(v == round(v, DECIMALS) for v in floats)
        # as the property is, to the decimals kept
        assert record["radius"] == pytest.approx(MolGroup(uni, [1, 2]).radius, abs=1e-6)

    def test_a_group_counts_its_residue_names(self, make_universe: UniverseFactory):
        uni = make_universe([[[1, 2, 3]]], 3, ["SOL", "MOL", "SOL"])

        record = group_record(MolGroup(uni, [3, 1, 2]))

        assert record["resids"] == [1, 2, 3]
        assert record["composition"] == {"SOL": 2, "MOL": 1}

    def test_a_cluster_has_its_birth_time(self, make_universe: UniverseFactory):
        uni = make_universe([[], [[1, 2]]], 2)
        tracker = ClusterTracker(uni, MolClsConfig(rules=MOL_RULES))
        uni.trajectory[1]
        tracker.update()
        (cluster,) = tracker.clusters.values()

        record = cluster_record(cluster)

        assert list(record)[:3] == ["id", "birth_time", "size"]
        assert record["birth_time"] == 1.0

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
        assert header["analyses"] == ["JsonReport"]
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


class Sizes(FrameAnalysis):
    """A user analysis adding to the report: frame counts and cluster centers."""

    def analyse(self, frame: Frame) -> None:
        pass

    def report_frame(self, frame: Frame) -> dict[str, object]:
        sizes = [c.size for c in frame.clusters.values()]
        return {"largest": max(sizes, default=0), "n": np.int64(len(sizes))}

    def report_cluster(self, frame: Frame, cluster: Cluster) -> dict | None:
        if cluster.size < 3:
            return None  # adds nothing
        return {"center": cluster.center_of_mass, "thirds": (1 / 3, [2 / 3])}


class TestContributions:
    def run(
        self,
        make_universe: UniverseFactory,
        directory: Path,
        analyses: list[FrameAnalysis],
        **config: object,
    ) -> Report:
        uni = make_universe(FRAMES, 4)
        config = MolClsConfig(rules=MOL_RULES, **config)
        MolClusters(uni, config, analyses).run(directory)
        return read_report(directory / "molclusters.jsonl.zst")

    def test_go_under_the_analysis_name(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        report = self.run(make_universe, tmp_path, [Sizes()])

        assert report.header["contributors"] == ["Lineage", "Sizes"]
        first, second, third = report.frames()
        assert first["Sizes"] == {"largest": 3, "n": 1}
        assert third["Sizes"] == {"largest": 0, "n": 0}
        (cluster,) = first["clusters"]
        center, thirds = cluster["Sizes"]["center"], cluster["Sizes"]["thirds"]
        assert len(center) == 3
        assert all(value == round(value, DECIMALS) for value in center)
        assert thirds == [0.333333, [0.666667]]
        # None adds nothing
        assert all("Sizes" not in cluster for cluster in second["clusters"])

    def test_come_after_the_builtins(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        report = self.run(make_universe, tmp_path, [Sizes()], nucleus=["MOL"])

        assert report.header["analyses"] == [
            "SizeEvolution",
            "Lineage",
            "Nucleus",
            "Sizes",
            "JsonReport",
        ]
        assert report.header["contributors"] == ["Lineage", "Nucleus", "Sizes"]
        (cluster,) = next(report.frames())["clusters"]
        assert list(cluster)[-2:] == ["Nucleus", "Sizes"]

    def test_an_analysis_not_adding_is_no_contributor(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        report = self.run(make_universe, tmp_path, [Interrupt(frame=-1)])

        # it ran all the same
        assert report.header["analyses"] == [
            "SizeEvolution",
            "Lineage",
            "Interrupt",
            "JsonReport",
        ]
        assert report.header["contributors"] == ["Lineage"]

    def test_what_cant_be_written_is_blamed_on_its_analysis(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        class Unwritable(FrameAnalysis):
            def analyse(self, frame: Frame) -> None:
                pass

            def report_cluster(self, frame: Frame, cluster: Cluster) -> dict:
                return {"members": set(cluster)}

        with pytest.raises(TypeError, match="A set can't be written") as err:
            self.run(make_universe, tmp_path, [Unwritable()])

        assert err.value.__notes__ == [
            "In Unwritable.report_cluster(), for the report",
            "Raised by JsonReport.analyse() on frame 0",
        ]

    def test_two_of_one_name_are_refused(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        with pytest.raises(ValueError, match="More than one Sizes analysis"):
            self.run(make_universe, tmp_path, [Sizes(), Sizes()])


class TestRounded:
    def test_rounds_every_float_however_nested(self):
        value = {"a": [1 / 3, (2 / 3, {"b": np.float32(0.1234567)})], "c": "x"}

        expected = {"a": [0.333333, [0.666667, {"b": 0.123457}]], "c": "x"}
        assert rounded(value) == expected

    def test_keeps_whole_numbers_and_arrays(self):
        ints = np.arange(3)[::2]  # not contiguous

        assert rounded(np.int64(2)) == 2
        assert rounded(None) is None
        assert rounded(True) is True
        np.testing.assert_array_equal(rounded(ints), [0, 2])
        np.testing.assert_array_equal(rounded(np.array([1 / 3])), [0.333333])
        assert dumps({"a": rounded(ints)}) == b'{"a":[0,2]}\n'

    def test_integer_keys_become_strings(self):
        assert rounded({1: "a", np.int64(2): "b"}) == {"1": "a", "2": "b"}

    def test_refuses_what_json_cant_hold(self):
        with pytest.raises(TypeError, match="A set can't be written"):
            rounded({"a": {1, 2}})
        with pytest.raises(TypeError, match="A tuple key can't be written"):
            rounded({(1, 2): 0.5})
        with pytest.raises(TypeError, match="A bool key can't be written"):
            rounded({True: 0.5})
