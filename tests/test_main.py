# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

import runpy
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from MDAnalysis import Universe

import molclusters.main as cli
from molclusters.config import MolClsConfig
from molclusters.main import _apply_lammps_resnames
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
        assert fake_universe["kwargs"] == {"in_memory": False, "in_memory_step": 1}

        uni = fake_universe["uni"]
        assert uni.atoms[0].element == "C"
        assert uni.atoms.radii[:2].tolist() == [1.7, 1.52]

        df = pd.read_csv(tmp_path / "solute_solvent.csv")
        assert df["NSolv"].tolist() == [1, 1]
        assert (tmp_path / "molclusters.json").exists()

    def test_solvent_defaults_to_every_non_solute_resname(
        self, cli_env: RunCli, fake_universe: FakeUniverse, tmp_path: Path
    ):
        cli_env("traj.xtc", "top.tpr", str(write_config(tmp_path)))

        df = pd.read_csv(tmp_path / "solute_solvent.csv")
        assert df["NSolv"].tolist() == [1, 1]

    def test_in_memory_options_are_forwarded(
        self, cli_env: RunCli, fake_universe: FakeUniverse, config_file: Path
    ):
        cli_env(
            "traj.xtc", "top.tpr", str(config_file),
            "--traj-memory", "--in-memory-step", "2",
        )  # fmt: skip

        assert fake_universe["kwargs"] == {"in_memory": True, "in_memory_step": 2}

    def test_in_memory_step_requires_traj_memory(
        self, cli_env: RunCli, config_file: Path, capsys: pytest.CaptureFixture
    ):
        with pytest.raises(SystemExit) as exit_info:
            cli_env("traj.xtc", "top.tpr", str(config_file), "--in-memory-step", "2")

        assert exit_info.value.code == 2
        assert "--traj-memory" in capsys.readouterr().err

    def test_unknown_element_is_reported(
        self, cli_env: RunCli, fake_universe: FakeUniverse, config_file: Path
    ):
        fake_universe["uni"].atoms[1].element = "Qq"

        with pytest.raises(KeyError, match="does not have an element"):
            cli_env("traj.xtc", "top.tpr", str(config_file))

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


def test_version(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture):
    monkeypatch.setattr(sys, "argv", ["molclusters", "--version"])

    with pytest.raises(SystemExit) as exit_info:
        runpy.run_module("molclusters", run_name="__main__")

    assert exit_info.value.code == 0
    assert capsys.readouterr().out.strip() == __version__
