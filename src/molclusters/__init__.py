import importlib.metadata

from .molclusters import MolClusters
from . import cluster as cluster

__all__ = ["MolClusters"]

__version__ = importlib.metadata.version("MolClusters")
