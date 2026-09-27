# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

import runpy
import shlex
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from MDAnalysis import Universe
from pydantic import ValidationError

import molclusters.main as cli
from molclusters.config import MolClsConfig
from molclusters.main import (
    _apply_lammps_dump_elements,
    _apply_lammps_resnames,
    _assign_radii,
    _describe_error,
    _lammps_dump_timestep,
    _log_system,
    _traj_format,
)
from molclusters.version import __version__

from .conftest import CUTOFF, UniverseFactory

DATA_DIR = (Path(__file__).parent / "data").resolve()

type RunCli = Callable[..., None]
type FakeUniverse = dict[str, Any]


class TestApplyLammpsResnames:
    def test_noop_when_mapping_not_configured(self):
        uni = Universe(str(DATA_DIR / "lammps_mini.data"))
        config = MolClsConfig(rules={"SOL": {"SOL": "cm 5.0"}})

        _apply_lammps_resnames(uni, config)

        assert not hasattr(uni.residues, "resnames")

    def test_assigns_resnames_from_molecule_id_mapping(self):
        uni = Universe(str(DATA_DIR / "lammps_mini.data"))
        config = MolClsConfig(
            rules={"SOL": {"SOL": "cm 5.0"}},
            lammps_resnames={"SOL": "1-2"},
        )

        _apply_lammps_resnames(uni, config)

        assert list(uni.residues.resnames) == ["SOL", "SOL"]

    def test_overwrites_existing_resnames(self, make_universe: UniverseFactory):
        uni = make_universe([[]], 3)
        config = MolClsConfig(
            rules={"SOL": {"SOL": "cm 5.0"}},
            lammps_resnames={"SOL": "1-2", "NA": 3},
        )

        _apply_lammps_resnames(uni, config)

        assert list(uni.residues.resnames) == ["SOL", "SOL", "NA"]

    def test_raises_when_a_molecule_id_is_not_covered(self):
        uni = Universe(str(DATA_DIR / "lammps_mini.data"))
        config = MolClsConfig(
            rules={"SOL": {"SOL": "cm 5.0"}},
            lammps_resnames={"SOL": 1},
        )

        with pytest.raises(ValueError, match="does not cover"):
            _apply_lammps_resnames(uni, config)


LAMMPS_DATA = str(DATA_DIR / "lammps_mini.data")
LAMMPS_DUMP = str(DATA_DIR / "lammps_mini.lammpsdump")


def write_dump(tmp_path: Path, columns: str, rows: list[str]) -> str:
    path = tmp_path / "traj.lammpsdump"
    path.write_text(
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n"
        f"{len(rows)}\n"
        "ITEM: BOX BOUNDS pp pp pp\n0.0 20.0\n0.0 20.0\n0.0 20.0\n"
        f"ITEM: ATOMS {columns}\n" + "\n".join(rows) + "\n"
    )
    return str(path)


@pytest.mark.parametrize(
    "traj",
    ["a.lammpsdump", "a.lammpstrj", "a.dump", "a.lammpstrj.gz", "run.v2/a.DUMP"],
)
def test_traj_format_reads_lammps_dump_extensions(traj: str):
    assert _traj_format(traj) == "LAMMPSDUMP"


@pytest.mark.parametrize("traj", ["a.xtc", "a.trr", "a.data", "traj"])
def test_traj_format_leaves_other_files_to_mdanalysis(traj: str):
    assert _traj_format(traj) is None


class TestApplyLammpsDumpElements:
    def test_assigns_elements_by_atom_id(self):
        # the dump lists the atoms out of id order
        uni = Universe(LAMMPS_DATA, LAMMPS_DUMP)

        _apply_lammps_dump_elements(uni, LAMMPS_DUMP)

        assert list(uni.atoms.elements) == ["O", "H", "N", "C"]

    def test_noop_when_trajectory_is_not_a_lammps_dump(self):
        uni = Universe(LAMMPS_DATA)

        _apply_lammps_dump_elements(uni, "traj.xtc")

        assert not hasattr(uni.atoms, "elements")

    def test_noop_when_dump_has_no_element_column(self, tmp_path: Path):
        dump = write_dump(
            tmp_path, "id mol type x y z",
            ["1 1 1 0 0 0", "2 1 1 1 0 0", "3 2 1 10 0 0", "4 2 1 11 0 0"],
        )  # fmt: skip
        uni = Universe(LAMMPS_DATA, dump)

        _apply_lammps_dump_elements(uni, dump)

        assert not hasattr(uni.atoms, "elements")

    def test_keeps_existing_elements(self, make_universe: UniverseFactory):
        uni = make_universe([[]], 2)
        elements = list(uni.atoms.elements)

        _apply_lammps_dump_elements(uni, LAMMPS_DUMP)

        assert list(uni.atoms.elements) == elements

    def test_raises_when_atom_ids_do_not_match(self, tmp_path: Path):
        dump = write_dump(
            tmp_path, "id mol type element x y z",
            ["1 1 1 O 0 0 0", "2 1 1 H 1 0 0", "3 2 1 N 10 0 0", "5 2 1 C 11 0 0"],
        )  # fmt: skip
        uni = Universe(LAMMPS_DATA)

        with pytest.raises(ValueError, match="don't match"):
            _apply_lammps_dump_elements(uni, dump)


class TestLammpsDumpTimestep:
    CONFIG = {"rules": {"A": {"A": "cm 5.0"}}}

    def test_is_passed_for_a_dump(self, captured_logs: list[str]):
        config = MolClsConfig(**self.CONFIG, lammps_timestep="2 fs")

        assert _lammps_dump_timestep(config, "traj.lammpsdump") == {
            "dt": pytest.approx(0.002)
        }
        assert not any("lammps_timestep" in m for m in captured_logs)

    def test_without_it_dump_times_are_step_numbers(self, captured_logs: list[str]):
        config = MolClsConfig(**self.CONFIG)

        assert _lammps_dump_timestep(config, "traj.lammpstrj") == {}
        assert "will be a step number, not ps" in captured_logs[-1]

    def test_has_no_effect_on_other_trajectories(self, captured_logs: list[str]):
        config = MolClsConfig(**self.CONFIG, lammps_timestep="2 fs")

        assert _lammps_dump_timestep(config, "traj.xtc") == {}
        assert "'lammps_timestep' has no effect" in captured_logs[-1]

    def test_other_trajectories_need_nothing(self, captured_logs: list[str]):
        config = MolClsConfig(**self.CONFIG)

        assert _lammps_dump_timestep(config, "traj.xtc") == {}
        assert not any("lammps_timestep" in m for m in captured_logs)


def write_config(tmp_path: Path, extra: str = "") -> Path:
    path = tmp_path / "input.yaml"
    path.write_text(
        "rules:\n"
        "  MOL:\n"
        f"    MOL: cm {CUTOFF}\n"
        f"    SOL: cm {CUTOFF}\n"
        "solute: [MOL]\n" + extra
    )
    return path


@pytest.fixture
def config_file(tmp_path: Path) -> Path:
    return write_config(tmp_path, "solvent: [SOL]\n")


@pytest.fixture
def cli_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> RunCli:
    """Run the CLI from `tmp_path`, without installing global log sinks.

    Returns
    -------
    Callable
        ``run(*argv)`` sets ``sys.argv`` and calls `main()`.
    """
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(cli, "start_logging", lambda **_: None)

    def run(*argv: str) -> None:
        monkeypatch.setattr(sys, "argv", ["molclusters", *argv])
        cli.main()

    return run


@pytest.fixture
def fake_universe(
    monkeypatch: pytest.MonkeyPatch, make_universe: UniverseFactory
) -> FakeUniverse:
    """Make `main()` load a synthetic system instead of reading files.

    Returns
    -------
    dict
        Filled with the Universe handed to `main()` and the arguments it was built
        with, once `main()` runs.
    """
    # residue 1 is solute (MOL) bonded to solvent residue 2 (SOL); 3 sits next to 2,
    # but there is no SOL-SOL rule. Atom 0's element is blank, so main() has to
    # derive it from the atom name
    uni = make_universe(
        [[[1, 2, 3]], [[1, 2, 3]]], 4, ["MOL", "SOL", "SOL", "SOL"], blank_elements=[0]
    )
    seen = {"uni": uni}

    def factory(*args: str, **kwargs: object) -> Universe:
        seen.update(args=args, kwargs=kwargs)
        return uni

    monkeypatch.setattr(cli.mda, "Universe", factory)
    return seen


class TestMain:
    def test_runs_the_analysis(
        self,
        cli_env: RunCli,
        fake_universe: FakeUniverse,
        config_file: Path,
        tmp_path: Path,
    ):
        cli_env("traj.xtc", "top.tpr", str(config_file))

        assert fake_universe["args"] == ("top.tpr", "traj.xtc")
        assert fake_universe["kwargs"] == {"format": None}

        uni = fake_universe["uni"]
        assert uni.atoms[0].element == "C"
        assert uni.atoms.radii[:2].tolist() == [1.7, 1.52]

        df = pd.read_csv(tmp_path / "solute_solvent.csv")
        assert df["NSolv"].tolist() == [1, 1]
        assert (tmp_path / "molclusters.json").exists()

    def test_output_dir_holds_the_results_and_the_log(
        self,
        cli_env: RunCli,
        fake_universe: FakeUniverse,
        config_file: Path,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ):
        log_files = []
        monkeypatch.setattr(
            cli, "start_logging", lambda filename, **_: log_files.append(filename)
        )

        cli_env("traj.xtc", "top.tpr", str(config_file), "--output-dir", "out/run1")

        out = tmp_path / "out" / "run1"
        assert (out / "solute_solvent.csv").exists()
        assert (out / "molclusters.json").exists()
        assert not (tmp_path / "molclusters.json").exists()
        (log_file,) = log_files
        assert log_file.parent == Path("out/run1")

    def test_solvent_defaults_to_every_non_solute_resname(
        self, cli_env: RunCli, fake_universe: FakeUniverse, tmp_path: Path
    ):
        cli_env("traj.xtc", "top.tpr", str(write_config(tmp_path)))

        df = pd.read_csv(tmp_path / "solute_solvent.csv")
        assert df["NSolv"].tolist() == [1, 1]

    @pytest.mark.parametrize(
        ("options", "transfers"),
        [((), []), (("--traj-memory", "--in-memory-step", "2"), [2])],
    )
    def test_in_memory_options_are_forwarded(
        self,
        cli_env: RunCli,
        fake_universe: FakeUniverse,
        config_file: Path,
        monkeypatch: pytest.MonkeyPatch,
        options: tuple[str, ...],
        transfers: list[int],
    ):
        steps = []
        monkeypatch.setattr(
            fake_universe["uni"], "transfer_to_memory", lambda step: steps.append(step)
        )

        cli_env("traj.xtc", "top.tpr", str(config_file), *options)

        assert fake_universe["kwargs"] == {"format": None}
        assert steps == transfers

    def test_lammps_timestep_is_forwarded_for_a_dump(
        self, cli_env: RunCli, fake_universe: FakeUniverse, tmp_path: Path
    ):
        config = write_config(tmp_path, "solvent: [SOL]\nlammps_timestep: 2 fs\n")

        cli_env("traj.lammpsdump", "top.data", str(config))

        assert fake_universe["kwargs"]["format"] == "LAMMPSDUMP"
        assert fake_universe["kwargs"]["dt"] == pytest.approx(0.002)

    def test_in_memory_step_requires_traj_memory(
        self, cli_env: RunCli, config_file: Path, capsys: pytest.CaptureFixture
    ):
        with pytest.raises(SystemExit) as exit_info:
            cli_env("traj.xtc", "top.tpr", str(config_file), "--in-memory-step", "2")

        assert exit_info.value.code == 2
        assert "--traj-memory" in capsys.readouterr().err

    def test_log_level_is_forwarded(
        self,
        cli_env: RunCli,
        fake_universe: FakeUniverse,
        config_file: Path,
        monkeypatch: pytest.MonkeyPatch,
    ):
        levels = []
        monkeypatch.setattr(
            cli, "start_logging", lambda level, **_: levels.append(level)
        )

        cli_env("traj.xtc", "top.tpr", str(config_file), "--log-level", "debug")

        assert levels == ["DEBUG"]

    def test_errors_are_logged_and_exit(
        self,
        cli_env: RunCli,
        fake_universe: FakeUniverse,
        config_file: Path,
        captured_logs: list[str],
    ):
        fake_universe["uni"].atoms[1].element = "Qq"

        with pytest.raises(SystemExit) as exit_info:
            cli_env("traj.xtc", "top.tpr", str(config_file))

        assert exit_info.value.code == 1
        # the traceback is logged on its own, then a readable summary
        traceback, summary = captured_logs[-2:]
        assert traceback.startswith("Traceback:")
        assert summary.startswith("KeyError: No van der Waals radius")
        assert "Full traceback in" in summary

    def test_interrupt_is_logged(
        self,
        cli_env: RunCli,
        config_file: Path,
        monkeypatch: pytest.MonkeyPatch,
        captured_logs: list[str],
    ):
        def interrupt(*_: object, **__: object) -> None:
            raise KeyboardInterrupt

        monkeypatch.setattr(cli.mda, "Universe", interrupt)

        with pytest.raises(SystemExit) as exit_info:
            cli_env("traj.xtc", "top.tpr", str(config_file))

        assert exit_info.value.code == 130
        assert captured_logs[-1].startswith("Interrupted")

    def test_command_line_is_logged(
        self,
        cli_env: RunCli,
        fake_universe: FakeUniverse,
        config_file: Path,
        captured_logs: list[str],
        tmp_path: Path,
    ):
        cli_env("traj.xtc", "top.tpr", str(config_file), "--traj-memory")

        (logged,) = [m for m in captured_logs if m.startswith("Running in")]
        assert logged.rstrip() == (
            f"Running in {tmp_path}: molclusters traj.xtc top.tpr "
            f"{shlex.quote(str(config_file))} --traj-memory"
        )

    def test_runs_on_a_tpr_topology(self, cli_env: RunCli, tmp_path: Path):
        # TPR files have no elements, and atom names such as 'CMB' or 'HM3' can't
        # be turned into one by just dropping their digits
        config = tmp_path / "input.yaml"
        config.write_text("rules:\n  MOL:\n    MOL: cm 5.0\n    MAL: cm 5.0\n")

        cli_env(
            str(DATA_DIR / "met-mal" / "start.pdb"),
            str(DATA_DIR / "met-mal" / "met-mal.tpr"),
            str(config),
        )

        assert (tmp_path / "molclusters.json").exists()

    @pytest.mark.parametrize("ext", ["lammpsdump", "lammpstrj", "dump"])
    def test_runs_on_a_lammps_data_topology_with_a_dump(
        self, cli_env: RunCli, tmp_path: Path, ext: str
    ):
        # DATA files have neither elements nor atom names: they come from the dump
        dump = tmp_path / f"traj.{ext}"
        dump.write_text(Path(LAMMPS_DUMP).read_text())
        config = tmp_path / "input.yaml"
        config.write_text(
            "rules:\n  SOL:\n    SOL: cm 5.0\nlammps_resnames:\n  SOL: 1-2\n"
        )

        cli_env(str(dump), LAMMPS_DATA, str(config))

        assert (tmp_path / "molclusters.json").exists()

    def test_lammps_dump_without_elements_is_reported(
        self, cli_env: RunCli, tmp_path: Path, captured_logs: list[str]
    ):
        dump = write_dump(
            tmp_path, "id mol type x y z",
            ["1 1 1 0 0 0", "2 1 1 1 0 0", "3 2 1 10 0 0", "4 2 1 11 0 0"],
        )  # fmt: skip
        config = tmp_path / "input.yaml"
        config.write_text(
            "rules:\n  SOL:\n    SOL: cm 5.0\nlammps_resnames:\n  SOL: 1-2\n"
        )

        with pytest.raises(SystemExit):
            cli_env(dump, LAMMPS_DATA, str(config))

        assert "'element' column" in captured_logs[-1]


class TestDescribeError:
    def test_names_the_error_type(self):
        assert _describe_error(ValueError("bad value")) == "ValueError: bad value"

    def test_key_error_message_is_not_quoted(self):
        assert _describe_error(KeyError("no such key")) == "KeyError: no such key"

    def test_notes_follow_the_message(self):
        err = ValueError("boom")
        err.add_note("Raised by Broken.analyse() on frame 1")

        assert _describe_error(err) == (
            "ValueError: boom\nRaised by Broken.analyse() on frame 1"
        )

    @pytest.mark.parametrize(
        ("config", "expected"),
        [
            (
                {"rules": {"MOL": {"MOL": "cm 3"}}, "distance_backend": "GPU"},
                "\n  distance_backend: Input should be 'serial' or 'OpenMP'",
            ),
            (
                {"rules": {"MOL": {"MOL": "cm abc"}}},
                "\n  Invalid rule MOL:MOL ('cm abc'): Invalid number for 'cm'",
            ),
        ],
    )
    def test_lists_config_validation_errors_plainly(
        self, config: dict[str, Any], expected: str
    ):
        with pytest.raises(ValidationError) as err_info:
            MolClsConfig(**config)

        described = _describe_error(err_info.value)

        assert described.startswith("Invalid configuration:\n")
        assert expected in described
        assert "pydantic.dev" not in described

    def test_lists_each_error_of_a_group(self):
        err = ExceptionGroup(
            "Errors in input file", [ValueError("first"), ValueError("second")]
        )

        assert _describe_error(err) == (
            "Errors in input file:\n  ValueError: first\n  ValueError: second"
        )


class TestAssignRadii:
    def test_takes_empty_elements_from_names(
        self, make_universe: UniverseFactory, captured_logs: list[str]
    ):
        uni = make_universe([[]], 2, blank_elements=[0, 2])

        _assign_radii(uni)

        assert uni.atoms.elements.tolist() == ["C", "O", "C", "O"]
        assert uni.atoms.radii.tolist() == [1.7, 1.52, 1.7, 1.52]
        assert captured_logs[-1].startswith("2 atom(s) had an empty element")

    def test_lists_every_unknown_element(self, make_universe: UniverseFactory):
        uni = make_universe([[]], 2)
        uni.atoms[1].element = "Qq"
        uni.atoms[3].element = "Xx"

        with pytest.raises(KeyError) as err_info:
            _assign_radii(uni)

        message = err_info.value.args[0]
        assert "element(s) ['Qq', 'Xx']" in message
        assert "'Qq', e.g. <Atom 2: O1" in message

    def test_reports_atoms_left_without_element(self, make_universe: UniverseFactory):
        uni = make_universe([[]], 1, blank_elements=[0])
        uni.atoms[0].name = "1"  # no letters to take an element from

        with pytest.raises(KeyError, match="no element, e.g. <Atom 1: 1"):
            _assign_radii(uni)


class TestLogSystem:
    def test_logs_composition_times_and_box(
        self, make_universe: UniverseFactory, captured_logs: list[str]
    ):
        uni = make_universe([[], [], []], 3, ["MOL", "SOL", "SOL"])

        _log_system(uni)

        assert [m.rstrip() for m in captured_logs] == [
            "System: 6 atoms in 3 residues (1 MOL, 2 SOL)",
            "Trajectory: 3 frame(s) from 0 to 2 ps, every 1 ps",
            "Box (first frame): 2000.00 x 2000.00 x 2000.00 angstrom",
        ]

    def test_warns_without_box(
        self, make_universe: UniverseFactory, captured_logs: list[str]
    ):
        uni = make_universe([[]], 1)
        uni.dimensions = None

        _log_system(uni)

        assert captured_logs[-1].startswith("The trajectory has no box")


def test_version(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture):
    monkeypatch.setattr(sys, "argv", ["molclusters", "--version"])

    with pytest.raises(SystemExit) as exit_info:
        runpy.run_module("molclusters", run_name="__main__")

    assert exit_info.value.code == 0
    assert capsys.readouterr().out.strip() == __version__
