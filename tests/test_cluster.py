# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

import numpy as np
import pytest
from MDAnalysis import Universe

from molclusters.cluster import Cluster


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

    @pytest.mark.skip()
    def test_distance_between_mols(self):
        assert self.cls.get_dist(0, 1) == 1.0

    def test_id(self):
        assert self.cls.id == 1

    @pytest.mark.skip()
    def test_contains(self):
        assert 0 in self.cls

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
