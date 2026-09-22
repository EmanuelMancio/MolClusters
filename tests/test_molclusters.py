# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

from pathlib import Path

import pytest
from MDAnalysis import Universe

from molclusters.config import MolClsConfig
from molclusters.molclusters import MolClusters

DATA_DIR = (Path(__file__).parent / "data" / "met-mal").resolve()


class TestWriteCoordinates:
    def test_writes_both_size_and_id_grouped_files(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ):
        uni = Universe(str(DATA_DIR / "met-mal.tpr"), str(DATA_DIR / "start.pdb"))
        config = MolClsConfig(rules={"MOL": {"MOL": "cm 15.0"}}, solute=["MOL"])
        molcls = MolClusters(uni, config)

        solute_clusters = [
            cls
            for cls in molcls.clusters.values()
            if set(config.solute).intersection(cls.resnames)
        ]
        assert solute_clusters, "fixture/rule setup should yield a solute cluster"

        monkeypatch.chdir(tmp_path)
        molcls._MolClusters__write_coordinates()

        size_files = list(tmp_path.glob("cls-n*.gro"))
        id_files = list(tmp_path.glob("cls-id*.gro"))

        assert size_files, "expected a size-grouped .gro output"
        assert len(id_files) == len(solute_clusters), (
            "expected one id-grouped .gro output per solute cluster"
        )
