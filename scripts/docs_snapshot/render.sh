#!/bin/sh
# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# Ray-trace the home page's picture on white and on black, for matte.py.
# Usage (Linux or WSL, VMD 1.9.3 on the PATH): render.sh IN.pdb WORKDIR
# Tachyon writes no alpha channel, so the scene is rendered twice, identically but
# for the background, and matte.py recovers the transparency from the two.
set -e
here=$(cd "$(dirname "$0")" && pwd)
pdb=$1
work=$2
tachyon=${TACHYON:-$(dirname "$(readlink -f "$(command -v vmd)")")/../lib/vmd/tachyon_LINUXAMD64}
[ -x "$tachyon" ] || tachyon=/usr/local/lib/vmd/tachyon_LINUXAMD64
vmd -dispdev text -e "$here/scene.tcl" -args "$pdb" "$work/scene.dat" < /dev/null > "$work/vmd.log" 2>&1
sed 's/^Background 1 1 1/Background 0 0 0/' "$work/scene.dat" > "$work/scene-black.dat"
for bg in white black; do
    dat="$work/scene.dat"
    [ "$bg" = black ] && dat="$work/scene-black.dat"
    "$tachyon" "$dat" -res 3000 3000 -aasamples 16 -format TARGA -o "$work/$bg.tga" > /dev/null
done
