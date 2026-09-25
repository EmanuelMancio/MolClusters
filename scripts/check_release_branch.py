# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Guard `cz bump` so it only runs on the release branch."""

import shutil
import subprocess  # noqa: S404
import sys

RELEASE_BRANCH = "main"

git = shutil.which("git")
branch = subprocess.run(  # noqa: S603
    [git, "branch", "--show-current"], capture_output=True, text=True, check=True
).stdout.strip()

if branch != RELEASE_BRANCH:
    print(f"cz bump must run on '{RELEASE_BRANCH}', not '{branch}'.", file=sys.stderr)
    sys.exit(1)
