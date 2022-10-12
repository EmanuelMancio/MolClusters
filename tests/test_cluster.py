from pathlib import Path
import sys
import pytest

directory = Path(__file__).absolute()

sys.path.append(str(directory.parent.parent))

from molclusters import Cluster


class TestCluster:
    cls = Cluster(0, 1, "A", "B", 1.0, 0.0)

    def test_cluster_creation(self):
        assert self.cls is not None

    def test_cluster_length(self):
        assert len(self.cls) == 2

    def test_distance_between_mols(self):
        assert self.cls.get_dist(0, 1) == 1.0

    def test_id(self):
        assert self.cls.id == 1

    def test_contains(self):
        assert 0 in self.cls

    def test_add_con_exception(self):
        with pytest.raises(ValueError):
            self.cls.add_con(0, 2, 2.0)
