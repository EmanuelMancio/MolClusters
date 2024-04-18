from MolClusters import Cluster
from MDAnalysis import Universe

import pytest
import numpy as np


class TestCluster:
    cls = Cluster(
        Universe(
            "../data/met-mal/met-mal.tpr",
            "../data/met-mal/start.pdb",
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

    def test_print_cluster_np(self):
        print(np.array(self.cls.cluster))