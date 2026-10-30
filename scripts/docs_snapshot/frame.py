# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Write one frame of the met-mal test run as a PDB for VMD, for the home page's picture.

Usage: ``python frame.py REPORT FRAME OUT.pdb``, with REPORT the ``molclusters.jsonl``
of a run on ``tests/data/met-mal`` (see docs/developer-guide/docs.md). The droplet
is put at the box center, every cluster whole and next to the largest one, and only
its clusters and the free molecules within 12 Å of them are written. The B-factor
column holds each molecule's cluster rank by size (1 the largest) or 0 if free;
scene.tcl colours by it.
"""

import json
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np

DATA = Path(__file__).resolve().parents[2] / "tests" / "data" / "met-mal"


def main() -> None:
    """Write the PDB."""
    report, frame, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    u = mda.Universe(str(DATA / "met-mal.tpr"), str(DATA / "met-traj.trr"))
    u.trajectory[frame]
    box = u.dimensions[:3]
    with Path(report).open(encoding="utf-8") as f:
        clusters = next(
            d["clusters"] for d in map(json.loads, f) if d.get("frame") == frame
        )
    clusters.sort(key=lambda c: -c["size"])

    u.atoms.unwrap(compound="fragments")  # each molecule whole
    res = {r.resid: r for r in u.residues}
    rank = dict.fromkeys(res, 0)
    center = None
    for k, c in enumerate(clusters, start=1):
        ref = res[c["resids"][0]].atoms.center_of_mass()
        for rid in c["resids"]:  # each molecule at its nearest image to the first one
            ag = res[rid].atoms
            ag.positions -= np.round((ag.center_of_mass() - ref) / box) * box
            rank[rid] = k
        atoms = sum((res[r].atoms for r in c["resids"][1:]), res[c["resids"][0]].atoms)
        if center is None:
            center = atoms.center_of_mass()
        atoms.positions -= np.round((atoms.center_of_mass() - center) / box) * box
    for rid, r in res.items():
        if not rank[rid]:
            r.atoms.positions -= (
                np.round((r.atoms.center_of_mass() - center) / box) * box
            )
    u.atoms.positions += box / 2 - center

    u.add_TopologyAttr("tempfactors")
    for rid, r in res.items():
        r.atoms.tempfactors = rank[rid]
    in_clusters = u.atoms[u.atoms.tempfactors > 0]
    near = u.select_atoms(
        "tempfactor 0 and same residue as around 12 group cls", cls=in_clusters
    )
    (in_clusters | near).write(out, bonds="all")


if __name__ == "__main__":
    main()
