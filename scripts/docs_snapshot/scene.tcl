# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# The home page's picture: a VMD scene of frame.py's PDB, written as a Tachyon file.
# Usage: vmd -dispdev text -e scene.tcl -args IN.pdb OUT.dat
# (no `display resize`: VMD 1.9.3's text mode exits on it; render.sh sets the size)

lassign $argv pdb out
mol new $pdb waitfor all
mol delrep 0 top

# The six clusters (by size rank) most visible from this view, found by rendering
# each cluster in red against the others in grey; yellow is kept off orange, the
# palette's least distinct pair. Other clusters are grey, free molecules faint.
set shown {2 6 3 14 1 15}
set palette {
  {0.165 0.471 0.839}
  {0.922 0.408 0.204}
  {0.106 0.686 0.478}
  {0.929 0.631 0.000}
  {0.910 0.482 0.643}
  {0.384 0.314 0.839}
}
proc addrep {cid rgb sel material} {
  color change rgb $cid {*}$rgb
  mol representation VDW 1.0 24
  mol color ColorID $cid
  mol selection $sel
  mol material $material
  mol addrep top
}
foreach b $shown rgb $palette cid {23 24 25 26 27 28} {
  addrep $cid $rgb "beta $b" AOChalky
}
addrep 29 {0.537 0.529 0.506} "beta > 0 and not beta $shown" AOChalky
addrep 30 {0.80 0.80 0.78} "beta 0" Transparent
material change opacity Transparent 0.3

color Display Background white
display projection orthographic
display depthcue off
display shadows on
display ambientocclusion on
display aoambient 0.85
display aodirect 0.45
axes location off
display resetview
rotate x by -15
rotate y by 20
scale by 1.15
render Tachyon $out
quit
