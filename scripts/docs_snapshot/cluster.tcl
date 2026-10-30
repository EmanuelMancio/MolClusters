# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# One cluster at a few times, cluster_row.py's PDB, as a Tachyon scene file: the
# Output files page's picture of a cls-id<id>.gro file.
# Usage: vmd -dispdev text -e cluster.tcl -args IN.pdb OUT.dat
# Malic acid (MAL, the solute) blue, 2-methylerythritol (MOL) grey.

lassign $argv pdb out
mol new $pdb waitfor all
mol delrep 0 top
foreach {sel cid rgb} {
  "resname MAL" 23 {0.165 0.471 0.839}
  "resname MOL" 24 {0.537 0.529 0.506}
} {
  color change rgb $cid {*}$rgb
  mol representation VDW 1.0 24
  mol color ColorID $cid
  mol selection $sel
  mol material AOChalky
  mol addrep top
}
color Display Background white
display projection orthographic
display depthcue off
display shadows on
display ambientocclusion on
display aoambient 0.85
display aodirect 0.45
axes location off
display resetview
scale by 1.75
render Tachyon $out
quit
