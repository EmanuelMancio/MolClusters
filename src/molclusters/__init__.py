import importlib.metadata

__version__ = importlib.metadata.version("MolClusters")

from .molclusters import MolClusters
from . import cluster as cluster

__all__ = ["MolClusters"]
