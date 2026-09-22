# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

import numpy as np
import pytest
from MDAnalysis import Universe

from molclusters.cluster import Cluster
from molclusters.config import MolClsConfig
from molclusters.conntable import ConnectionTable


@pytest.fixture
def populated_cluster() -> Cluster:
    """Build a real two-molecule Cluster via the same path production code uses.

    Returns
    -------
    Cluster
        A cluster containing the two connected MOL residues from the fixture.
    """
    uni = Universe(
        "tests/data/met-mal/met-mal.tpr",
        "tests/data/met-mal/start.pdb",
    )
    config = MolClsConfig(rules={"MOL": {"MOL": "cm 15.0"}})
    sels = {res: uni.select_atoms(f"resname {res}") for res in config._rules.all_keys()}
    conntab = ConnectionTable(uni, config._rules, sels)

    return Cluster(uni, next(conntab.subconntables()))


class TestCluster:
    cls = Cluster(
        Universe(
            "tests/data/met-mal/met-mal.tpr",
            "tests/data/met-mal/start.pdb",
        )
    )

    def test_cluster_creation(self):
        assert self.cls is not None

    def test_cluster_length(self):
        assert len(self.cls) == 0

    def test_id(self):
        assert self.cls.id == 1

    def test_add_con_exception(self):
        with pytest.raises(ValueError):
            self.cls.add_con(0, 2, 2.0)

    def test_set_dist_exception_moli_not_in_cluster(self):
        with pytest.raises(ValueError):
            self.cls.set_dist(0, 1, 2.0)

    def test_add_unsupported_operand(self):
        with pytest.raises(TypeError):
            self.cls + "not a residue group"

    def test_print_cluster_np(self):
        print(np.array(self.cls.cluster))


class TestPopulatedCluster:
    def test_distance_between_mols(self, populated_cluster: Cluster):
        moli = next(iter(populated_cluster))
        molj = next(iter(populated_cluster.neighbors(moli)))

        expected = populated_cluster.cluster[moli][molj]["distance"]
        assert populated_cluster.get_dist(moli, molj) == expected

    def test_contains(self, populated_cluster: Cluster):
        mol = next(iter(populated_cluster))

        assert mol in populated_cluster
        assert -1 not in populated_cluster
