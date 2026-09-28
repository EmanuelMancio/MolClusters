# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Set CITATION.cff's version and release date to the ones `cz bump` is making.

Run as a commitizen pre-bump hook, which passes the new version in
CZ_PRE_NEW_VERSION; `cz bump` commits with `-a`, so the change lands in the
release commit. GitHub's "Cite this repository" reads the file, and Zenodo
takes its authors, title, abstract, keywords and license for each release's DOI.
"""

import datetime
import os
import re
import sys
from pathlib import Path

CITATION = Path(__file__).resolve().parent.parent / "CITATION.cff"

version = os.environ.get("CZ_PRE_NEW_VERSION")
if not version:
    print("CZ_PRE_NEW_VERSION is not set; run this through `cz bump`.", file=sys.stderr)
    sys.exit(1)

text = CITATION.read_text(encoding="utf-8")
for key, value in (
    ("version", version),
    ("date-released", f'"{datetime.date.today().isoformat()}"'),
):
    text, n = re.subn(rf"^{key}:.*$", f"{key}: {value}", text, flags=re.MULTILINE)
    if n != 1:
        print(f"CITATION.cff has {n} '{key}:' lines, expected 1.", file=sys.stderr)
        sys.exit(1)
CITATION.write_text(text, encoding="utf-8", newline="\n")
