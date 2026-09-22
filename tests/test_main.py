# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

from pathlib import Path

import pytest
from MDAnalysis import Universe

from molclusters.config import MolClsConfig
from molclusters.main import _apply_lammps_resnames

DATA_DIR = (Path(__file__).parent / "data").resolve()


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

    def test_raises_when_a_molecule_id_is_not_covered(self):
        uni = Universe(str(DATA_DIR / "lammps_mini.data"))
        config = MolClsConfig(
            rules={"SOL": {"SOL": "cm 5.0"}},
            lammps_resnames={"SOL": 1},
        )

        with pytest.raises(ValueError, match="does not cover"):
            _apply_lammps_resnames(uni, config)
