# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Put chosen frames of a ``cls-id<id>.gro`` file side by side, in one PDB for VMD.

Usage: ``python cluster_row.py IN.gro OUT.pdb SPACING TIME [TIME ...]``, times in ps
as the title lines give them (``Cluster-<id> - Time = <t>``), SPACING in Å between
the frames' centers. The frames of these files have different atom counts, which
VMD can't load as one trajectory, so they become one structure, left to right in
the order given, each centered on its own center of geometry.
"""

import io
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np


def main() -> None:
    """Write the PDB."""
    src, out, spacing, *times = sys.argv[1:]
    wanted = [float(t) for t in times]
    lines = Path(src).read_text(encoding="utf-8").splitlines(keepends=True)
    frames = {}
    i = 0
    while i < len(lines):
        n = int(lines[i + 1])
        time = float(lines[i].rsplit("=", 1)[1])
        if time in wanted:
            frames[time] = "".join(lines[i : i + n + 3])
        i += n + 3
    missing = [t for t in wanted if t not in frames]
    if missing:
        sys.exit(f"no frame at {missing} ps")
    parts = []
    for k, t in enumerate(wanted):
        u = mda.Universe(io.StringIO(frames[t]), format="GRO")
        u.atoms.positions -= u.atoms.positions.mean(axis=0)
        u.atoms.positions += np.array([k * float(spacing), 0, 0])
        parts.append(u.atoms)
    mda.Merge(*parts).atoms.write(out)


if __name__ == "__main__":
    main()
