# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

import io
import time
import warnings
from collections.abc import Sequence
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
import pytest
from loguru import logger
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.lib.util import NamedStream

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
from molclusters.analysis.coordinates import gro_frame
from molclusters.config import MolClsConfig, ReportCompression
from molclusters.output import OutputFile, RunOutput
from molclusters.report import Report, read_report
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

    def test_frame_gives_the_transition_from_the_previous_frame(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        uni = make_universe([[[1, 2, 3]], [[1, 2], [3, 4]]], 4)
        tracker = ClusterTracker(uni, MolClsConfig(rules=MOL_RULES))
        uni.trajectory[1]
        tracker.update()

        frame = Frame(1, tracker, RunOutput(tmp_path))

        assert frame.transition is tracker.transition
        (new,) = frame.transition.born
        assert frame.transition.sources(new) == {1: 1, 0: 1}

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


class TestLineage:
    # ids: {1, 2, 3} is 1 and {4, 5} is 2; 2 merges into 1 as {6, 7} forms (3),
    # then 1 splits {4, 5} off again (4), and both 3 and 4 dissolve
    FRAMES = [
        [[1, 2, 3], [4, 5]],
        [[1, 2, 3, 4, 5], [6, 7]],
        [[1, 2, 3], [4, 5], [6, 7]],
        [[1, 2, 3]],
    ]

    @staticmethod
    def events(directory: Path) -> list[str]:
        return (directory / "cluster_events.csv").read_text().splitlines()

    def test_writes_every_event_of_every_frame_but_the_first(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        run_analyses([Lineage()], make_universe, self.FRAMES, 7, tmp_path)

        assert self.events(tmp_path) == [
            "Frame,Time,Event,Cluster,Other,NMols",
            "1,1.0,formation,3,,2",
            "1,1.0,merge,1,2,2",
            "2,2.0,split,4,1,2",
            "3,3.0,dissolution,3,,2",
            "3,3.0,dissolution,4,,2",
        ]

    def test_records_each_cluster_birth_end_and_lifetime(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        lineage = Lineage()

        run_analyses([lineage], make_universe, self.FRAMES, 7, tmp_path)

        table = lineage.lifetimes()
        assert table.drop(columns=["DeathTime", "MergedInto"]).to_dict("list") == {
            "Id": [1, 2, 3, 4],
            "BirthTime": [0.0, 0.0, 1.0, 2.0],
            # the one alive at the end, up to the last frame
            "Lifetime": [3.0, 1.0, 2.0, 1.0],
            "NFrames": [4, 1, 2, 1],
            "BornAtStart": [True, True, False, False],
            "AliveAtEnd": [True, False, False, False],
            "Origin": ["initial", "initial", "formation", "split"],
            "Parents": ["", "", "", "1"],
            "Fate": ["alive", "merge", "dissolution", "dissolution"],
            "BirthSize": [3, 2, 2, 2],
            "MaxSize": [5, 2, 2, 2],
            "LastSize": [3, 2, 2, 2],
        }
        assert table["DeathTime"].tolist()[1:] == [1.0, 3.0, 3.0]
        assert np.isnan(table["DeathTime"][0])
        assert table["MergedInto"].tolist() == [pd.NA, 1, pd.NA, pd.NA]
        written = pd.read_csv(tmp_path / "cluster_lifetimes.csv")
        assert list(written.columns) == list(table.columns)
        assert written["Lifetime"].tolist() == table["Lifetime"].tolist()

    def test_merges_and_splits_of_several_clusters_give_a_row_each(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        # 2 and 3 merge into 1; 1 splits into three pieces (4 and 5 new); then one
        # molecule of 4 and one of 5 pair up (6), and the rest of both dissolves
        frames = [
            [[1, 2], [3, 4, 5], [6, 7, 8, 9]],
            [list(range(1, 10))],
            [[1, 2, 3, 4], [5, 6, 7], [8, 9]],
            [[1, 2, 3, 4], [7, 9]],
        ]
        lineage = Lineage()

        run_analyses([lineage], make_universe, frames, 9, tmp_path)

        assert self.events(tmp_path)[1:] == [
            "1,1.0,merge,1,2,3",
            "1,1.0,merge,1,3,2",
            "2,2.0,split,4,1,3",
            "2,2.0,split,5,1,2",
            "3,3.0,split,6,4,1",
            "3,3.0,split,6,5,1",
            "3,3.0,dissolution,4,,3",
            "3,3.0,dissolution,5,,2",
        ]
        table = lineage.lifetimes().set_index("Id")
        assert table.loc[6, "Parents"] == "4;5"
        assert table["MergedInto"].tolist()[:3] == [pd.NA, 1, 1]

    def test_finish_logs_a_summary(
        self,
        make_universe: UniverseFactory,
        tmp_path: Path,
        captured_logs: list[str],
    ):
        run_analyses([Lineage()], make_universe, self.FRAMES, 7, tmp_path)

        (summary,) = [m for m in captured_logs if m.startswith("Cluster lineage")]
        assert summary == (
            "Cluster lineage: 1 cluster(s) formed of free molecules and 1 split off "
            "others; 1 merged into another and 2 dissolved. The 2 born and ended "
            "within the run lived 1.5 ps (median); 1 of them for a single frame.\n"
        )

    def test_starts_over_on_every_run(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        lineage = Lineage()

        run_analyses([lineage], make_universe, self.FRAMES, 7, tmp_path)
        first = self.events(tmp_path)
        run_analyses([lineage], make_universe, self.FRAMES, 7, tmp_path)

        assert self.events(tmp_path) == first
        assert len(lineage.lifetimes()) == 4

    def test_adds_each_frame_transition_to_the_report(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        run_analyses([Lineage(), JsonReport()], make_universe, self.FRAMES, 7, tmp_path)

        report = read_report(tmp_path / "molclusters.jsonl.zst")
        assert report.header["contributors"] == ["Lineage"]
        first, second, *_ = report.frames()
        assert first["Lineage"] == {
            "flows": [[0, 1, 3], [0, 2, 2]],
            "born": [1, 2],
            "merged": [],
            "dissolved": [],
        }
        assert second["Lineage"] == {
            "flows": [[0, 3, 2], [1, 1, 3], [2, 1, 2]],
            "born": [3],
            "merged": [[2, 1]],
            "dissolved": [],
        }
        assert [(c["id"], c["birth_time"]) for c in second["clusters"]] == [
            (1, 0.0),
            (3, 1.0),
        ]

    def test_declares_its_files(self):
        assert [out.name for out in Lineage.outputs] == [
            "cluster_events.csv",
            "cluster_lifetimes.csv",
        ]


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


def written_by_mdanalysis(atoms: AtomGroup, title: str) -> str:
    """Write atoms with MDAnalysis' own GRO writer, with `title` for its title line.

    Returns
    -------
    str
        The file's text.
    """
    buf = io.StringIO()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # e.g. missing names, or box
        with mda.Writer(NamedStream(buf, "frame.gro"), multiframe=False) as writer:
            writer.write(atoms)
    _, body = buf.getvalue().split("\n", 1)
    return f"{title}\n{body}"


class TestGroFrame:
    """`gro_frame` renders what MDAnalysis' GRO writer writes, byte for byte."""

    @staticmethod
    def chain(make_universe: UniverseFactory, **kwargs: object) -> mda.Universe:
        uni = make_universe([[[1, 2, 3]]], 3, **kwargs)
        # negative coordinates, and digits past the 3 decimals written
        uni.atoms.positions = uni.atoms.positions - 104.123456
        return uni

    @pytest.mark.parametrize(
        "box",
        [
            [30.0, 40.0, 50.0, 90.0, 90.0, 90.0],
            [30.0, 40.0, 50.0, 80.0, 70.0, 60.0],  # triclinic
            None,
        ],
    )
    def test_matches_the_writer_for_any_box(
        self, make_universe: UniverseFactory, box: list[float] | None
    ):
        uni = self.chain(make_universe, box=box is not None)
        if box is not None:
            uni.dimensions = box

        assert gro_frame(uni.atoms, "t") == written_by_mdanalysis(uni.atoms, "t")

    def test_matches_the_writer_with_velocities(self, make_universe: UniverseFactory):
        uni = self.chain(make_universe)
        uni.trajectory.ts.has_velocities = True
        uni.atoms.velocities = np.linspace(-9.87654, 12.3456, 18).reshape(6, 3)

        text = gro_frame(uni.atoms, "t")

        assert text == written_by_mdanalysis(uni.atoms, "t")
        assert len(text.splitlines()[2]) == 44 + 24

    def test_numbers_the_atoms_in_the_order_given_and_truncates_as_the_writer(
        self, make_universe: UniverseFactory
    ):
        uni = self.chain(make_universe)
        uni.residues.resids = [123_456, 99_999, 100_000]
        uni.residues.resnames = ["LONGNAME", "A", "SOL"]
        uni.atoms.names = ["C1", "OXYGENS"] * 3
        atoms = uni.atoms[[5, 0, 3, 2]]

        assert gro_frame(atoms, "t") == written_by_mdanalysis(atoms, "t")

    def test_atoms_without_names_are_written_as_the_writer_does(self):
        uni = mda.Universe.empty(2, n_residues=1, atom_resindex=[0, 0], trajectory=True)
        uni.add_TopologyAttr("resnames", ["MOL"])
        uni.add_TopologyAttr("resids", [1])
        uni.atoms.positions = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
        uni.dimensions = [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]

        text = gro_frame(uni.atoms, "t")

        assert text == written_by_mdanalysis(uni.atoms, "t")
        assert "    X" in text

    def test_coordinates_out_of_the_format_range_are_rejected(
        self, make_universe: UniverseFactory
    ):
        uni = self.chain(make_universe)
        uni.atoms.positions = uni.atoms.positions + 200_000.0  # 20,000 nm or so

        with pytest.raises(ValueError, match="'Cluster-1' can't be written"):
            gro_frame(uni.atoms, "Cluster-1")


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
    ) -> Report:
        run_analyses(
            analyses,
            make_universe,
            self.FRAMES,
            4,
            directory,
            self.RESNAMES,
            ALL_PAIRS_RULES,
        )
        return read_report(directory / "molclusters.jsonl.zst")

    def test_writes_every_frame_and_cluster(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        report = self.run([JsonReport()], make_universe, tmp_path)

        assert report.header["software"].startswith("MolClusters ")
        assert report.header["config"]["rules"]
        assert report.header["n_frames"] == 2
        first, second = report.frames()
        assert (first["time"], first["frame"]) == (0, 0)
        (cluster,) = first["clusters"]
        assert cluster["resids"] == [1, 2, 3, 4]
        assert cluster["composition"] == {"MOL": 3, "SOL": 1}
        assert "Nucleus" not in cluster
        assert (second["time"], second["clusters"]) == (1, [])

    def test_includes_the_nuclei_of_a_nucleus_analysis_before_it(
        self, make_universe: UniverseFactory, tmp_path: Path
    ):
        report = self.run([Nucleus(["MOL"]), JsonReport()], make_universe, tmp_path)

        (cluster,) = next(report.frames())["clusters"]
        nuclei = cluster["Nucleus"]
        assert sorted(n["resids"] for n in nuclei["nuclei"]) == [[1, 2], [4]]
        assert nuclei["combined_dipole_moment"] is not None

    @pytest.mark.parametrize(
        ("compression", "name"),
        [
            ("zstd", "molclusters.jsonl.zst"),
            ("gzip", "molclusters.jsonl.gz"),
            ("none", "molclusters.jsonl"),
        ],
    )
    def test_declares_the_report_it_compresses(
        self,
        make_universe: UniverseFactory,
        tmp_path: Path,
        compression: ReportCompression,
        name: str,
    ):
        report = JsonReport(compression)
        run_analyses([report], make_universe, self.FRAMES, 4, tmp_path, self.RESNAMES)

        assert [out.name for out in report.outputs] == [name]
        assert [p.name for p in tmp_path.iterdir()] == [name]
        assert len(list(read_report(tmp_path / name).frames())) == 2


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
