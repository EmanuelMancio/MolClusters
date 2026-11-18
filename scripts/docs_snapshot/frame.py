# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Write a neighbourhood of one frame of the met-mal test run as a PDB for VMD.

Usage: ``python frame.py REPORT FRAME SEED OUT.pdb [RADIUS]``, with REPORT the
``molclusters.jsonl`` of a run on ``tests/data/met-mal`` (see
docs/developer-guide/docs.md), SEED the size rank (1 the largest) of the cluster
the view is centred on, and RADIUS (Å, default 40) how far around it molecules are
kept; pick_frame.py suggests FRAME and SEED. Every cluster within RADIUS is written
whole, at its nearest image to the seed, with the free molecules around it, and
the whole is turned so its clusters spread across the x-y plane. The B-factor
column holds each molecule's cluster rank by size within the view (1 the largest)
or 0 if free; scene.tcl colours by it.
"""

import json
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np

DATA = Path(__file__).resolve().parents[2] / "tests" / "data" / "met-mal"
RADIUS = 40.0


def universe() -> mda.Universe:
    """Load the test system.

    Returns
    -------
    MDAnalysis.Universe
        The met-mal topology and trajectory.
    """
    return mda.Universe(str(DATA / "met-mal.tpr"), str(DATA / "met-traj.trr"))


def read_report(report: str) -> dict[int, list[dict]]:
    """Read a run's report.

    Returns
    -------
    dict
        Each frame's clusters, by frame index, largest first.
    """
    with Path(report).open(encoding="utf-8") as f:
        frames = {d["frame"]: d["clusters"] for d in map(json.loads, f) if "frame" in d}
    for clusters in frames.values():
        clusters.sort(key=lambda c: -c["size"])
    return frames


def units(u: mda.Universe, clusters: list[dict]) -> list[mda.AtomGroup]:
    """Make each cluster whole; the current frame's positions are changed.

    Returns
    -------
    list of MDAnalysis.AtomGroup
        Each cluster's atoms, in the order given, then each free molecule's.
    """
    box = u.dimensions[:3]
    u.atoms.unwrap(compound="fragments")  # each molecule whole
    res = {r.resid: r for r in u.residues}
    out = []
    for c in clusters:
        ref = res[c["resids"][0]].atoms.center_of_mass()
        for rid in c["resids"]:  # each molecule at its nearest image to the first one
            ag = res.pop(rid).atoms
            ag.positions -= np.round((ag.center_of_mass() - ref) / box) * box
        out.append(u.residues[np.isin(u.residues.resids, c["resids"])].atoms)
    return out + [r.atoms for r in res.values()]


def near(groups: list[mda.AtomGroup], center: np.ndarray, radius: float) -> list[int]:
    """Move each group to its nearest image to center.

    Returns
    -------
    list of int
        The indices of the groups whose center of mass is within radius of center.
    """
    box = groups[0].dimensions[:3]
    kept = []
    for i, ag in enumerate(groups):
        ag.positions -= np.round((ag.center_of_mass() - center) / box) * box
        if np.linalg.norm(ag.center_of_mass() - center) <= radius:
            kept.append(i)
    return kept


def main() -> None:
    """Write the PDB."""
    report, frame, seed, out = (
        sys.argv[1],
        int(sys.argv[2]),
        int(sys.argv[3]),
        sys.argv[4],
    )
    radius = float(sys.argv[5]) if len(sys.argv) > 5 else RADIUS
    clusters = read_report(report)[frame]
    u = universe()
    u.trajectory[frame]
    groups = units(u, clusters)
    center = groups[seed - 1].center_of_mass()
    kept = near(groups, center, radius)
    shown = [i for i in kept if i < len(clusters)]  # largest first

    u.add_TopologyAttr("tempfactors")
    u.atoms.tempfactors = -1
    for rank, i in enumerate(shown, start=1):
        groups[i].tempfactors = rank
    for i in kept[len(shown) :]:
        groups[i].tempfactors = 0
    sel = u.atoms[u.atoms.tempfactors >= 0]
    # Turn the clusters' principal axes onto x, y, z, so that VMD, looking down z,
    # sees them spread as far apart as they are.
    xyz = sel.positions - sel.center_of_mass()
    axes = np.linalg.eigh(np.cov(xyz[sel.tempfactors > 0].T))[1][:, ::-1]
    if np.linalg.det(axes) < 0:
        axes[:, 2] *= -1
    sel.positions = xyz @ axes + u.dimensions[:3] / 2
    sel.write(out, bonds="all")
    sizes = [clusters[i]["size"] for i in shown]
    print(
        f"{len(shown)} clusters (sizes {sizes}), {len(kept) - len(shown)} free molecules"
    )


if __name__ == "__main__":
    main()
