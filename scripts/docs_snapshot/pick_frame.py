# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Suggest frames and views of the met-mal test run for the home page's picture.

Usage: ``python pick_frame.py REPORT [RADIUS] [STEP]``, with REPORT the
``molclusters.jsonl`` of a run on ``tests/data/met-mal``. The picture should show
clusters that are visibly apart, so that each colour is one group the reader can
see: for every STEP-th frame (default 2) and every cluster of 5 or more molecules,
the view of frame.py centred on it (RADIUS, default 40 Å) is scored by its
smallest gap between two clusters and between a cluster and a free molecule
(closest atoms, in Å). Views with 3 to 6 clusters (the palette's six colours)
holding at least 25 molecules and no gap under 5 Å are kept, and those with the
most molecules in clusters are printed as
``frame time seed sizes free gap free_gap``, for ``frame.py REPORT FRAME SEED``.
"""

import sys
from itertools import combinations

import numpy as np
from frame import RADIUS, near, read_report, units, universe
from MDAnalysis import AtomGroup
from MDAnalysis.lib.distances import distance_array


def gap(a: AtomGroup, b: AtomGroup, box: np.ndarray) -> float:
    """Measure how far apart two groups are.

    Returns
    -------
    float
        The distance between their closest atoms, minimum image (Å).
    """
    return float(distance_array(a.positions, b.positions, box=box).min())


def main() -> None:
    """Print the best views."""
    report = sys.argv[1]
    radius = float(sys.argv[2]) if len(sys.argv) > 2 else RADIUS
    step = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    frames = read_report(report)
    u = universe()
    found = []
    for ts in u.trajectory[::step]:
        clusters = frames[ts.frame]
        groups = units(u, clusters)
        box = u.dimensions
        for seed in [i for i, c in enumerate(clusters) if c["size"] >= 5]:
            kept = near(groups, groups[seed].center_of_mass(), radius)
            cls = [i for i in kept if i < len(clusters)]
            free = [i for i in kept if i >= len(clusters)]
            sizes = [clusters[i]["size"] for i in cls]
            if not 3 <= len(cls) <= 6 or sum(sizes) < 25:
                continue
            g = min(gap(groups[i], groups[j], box) for i, j in combinations(cls, 2))
            fg = min(
                (gap(groups[i], groups[j], box) for i in cls for j in free),
                default=np.inf,
            )
            if min(g, fg) >= 5:
                found.append(
                    (sum(sizes), g, fg, ts.frame, ts.time, seed + 1, sizes, len(free))
                )
    found.sort(key=lambda f: -f[0])
    print("frame  time/ps  seed  sizes  free  gap  free_gap")
    for _, g, fg, frame, time, seed, sizes, nfree in found[:15]:
        print(f"{frame} {time:.0f} {seed} {sizes} {nfree} {g:.1f} {fg:.1f}")


if __name__ == "__main__":
    main()
