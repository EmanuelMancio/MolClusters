# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# The home page's picture: a VMD scene of frame.py's PDB, written as a Tachyon file.
# Usage: vmd -dispdev text -e scene.tcl -args IN.pdb OUT.dat
# (no `display resize`: VMD 1.9.3's text mode exits on it; render.sh sets the size)

lassign $argv pdb out
mol new $pdb waitfor all
mol delrep 0 top

# Each cluster in the view (by size rank, at most six) in its own colour, free
# molecules smaller and grey.
set palette {
  {0.165 0.471 0.839}
  {0.922 0.408 0.204}
  {0.106 0.686 0.478}
  {0.929 0.631 0.000}
  {0.910 0.482 0.643}
  {0.384 0.314 0.839}
}
proc addrep {cid rgb sel scale} {
  color change rgb $cid {*}$rgb
  mol representation VDW $scale 24
  mol color ColorID $cid
  mol selection $sel
  mol material AOChalky
  mol addrep top
}
foreach rgb $palette cid {23 24 25 26 27 28} rank {1 2 3 4 5 6} {
  addrep $cid $rgb "beta $rank" 1.0
}
addrep 29 {0.62 0.62 0.60} "beta 0" 0.6

color Display Background white
display projection orthographic
display depthcue off
display shadows on
display ambientocclusion on
display aoambient 0.85
display aodirect 0.45
axes location off
display resetview
rotate x by -20
rotate y by 15
scale by 1.6
render Tachyon $out
quit
