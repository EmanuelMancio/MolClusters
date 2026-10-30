#!/bin/sh
# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# Smoke test of a built container image: that it doesn't run as root, `molclusters
# --version`, then a run on the met-mal test system as the calling user, as a user
# would run it. Run from the
# repository's root:
#
#   scripts/check_container.sh IMAGE [EXPECTED_VERSION]
#
# The engine is `docker`, or $ENGINE (e.g. ENGINE=podman).

set -eu

image=$1
expected=${2:-}
engine=${ENGINE:-docker}

user=$("$engine" image inspect --format '{{.Config.User}}' "$image")
case "$user" in
    "" | root | 0 | 0:* | root:*)
        echo "error: the image runs as root (USER '$user')" >&2; exit 1 ;;
esac

version=$("$engine" run --rm "$image" --version)
echo "$version"
if [ -n "$expected" ]; then
    case "$version" in
        *"$expected"*) ;;
        *) echo "error: expected version $expected" >&2; exit 1 ;;
    esac
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
cp tests/data/met-mal/start.pdb tests/data/met-mal/met-mal.tpr "$work"
printf 'rules:\n  MOL:\n    MOL: cm 5.0\n    MAL: cm 5.0\n' > "$work/input.yaml"

case "$engine" in
    podman) as_caller=--userns=keep-id ;;  # rootless Podman maps --user to a subuid
    *) as_caller="--user=$(id -u):$(id -g)" ;;
esac
"$engine" run --rm "$as_caller" -v "$work:/data" "$image" \
    start.pdb met-mal.tpr input.yaml

if [ ! -s "$work/molclusters.jsonl.zst" ]; then
    echo "error: the run wrote no molclusters.jsonl.zst" >&2
    ls -la "$work" >&2
    exit 1
fi
echo "Container image $image works."
