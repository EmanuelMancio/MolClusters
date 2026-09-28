# API reference

Generated from the docstrings of the `molclusters` package. Most users need only a
few names, all importable from the package itself:

```python
from molclusters import (
    MolClusters,        # run the analyses over a trajectory
    ClusterTracker,     # follow the clusters, without analyses
    Transition,         # how a frame's clusters came from the previous frame's
    FrameAnalysis,      # base class of an analysis
    Run, Frame,         # what an analysis sees of the run and of a frame
    OutputFile,         # declaration of a file an analysis writes
    Cluster, MolGroup,  # a tracked cluster, and any group of molecules
    start_logging,      # log like the command line does
)
from molclusters.config import MolClsConfig, read_config
from molclusters.report import read_report
```

| Module | Contents |
| --- | --- |
| [`molclusters`](molclusters.md) | `MolClusters`, the runner |
| [`molclusters.config`](config.md) | `MolClsConfig`, `read_config`, the rule classes |
| [`molclusters.tracker`](tracker.md) | `ClusterTracker`, `Transition` |
| [`molclusters.cluster`](cluster.md) | `MolGroup`, `Cluster` and their properties |
| [`molclusters.analysis`](analysis.md) | `FrameAnalysis`, `Run`, `Frame`, the built-in analyses and their registry |
| [`molclusters.output`](output.md) | `OutputFile`, `RunOutput` |
| [`molclusters.report`](report.md) | the JSON report's format and reader |
| [`molclusters.conntable`](conntable.md) | `ConnectionTable`, the connections of a frame |
| [`molclusters.log`](log.md) | `start_logging` |
