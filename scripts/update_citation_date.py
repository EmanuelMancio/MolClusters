# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Set CITATION.cff's release date to today, for the release `cz bump` is making.

Run as a commitizen pre-bump hook; commitizen itself sets the file's `version`
(`version_files`), but has no way to write a date. `cz bump` commits with `-a`,
so the change lands in the release commit. GitHub's "Cite this repository"
reads the file, and Zenodo takes its authors, title, abstract, keywords and
license for each release's DOI.
"""

import datetime
import re
import sys
from pathlib import Path

CITATION = Path(__file__).resolve().parent.parent / "CITATION.cff"

text = CITATION.read_text(encoding="utf-8")
text, n = re.subn(
    r"^date-released:.*$",
    f'date-released: "{datetime.date.today().isoformat()}"',
    text,
    flags=re.MULTILINE,
)
if n != 1:
    print(f"CITATION.cff has {n} 'date-released:' lines, expected 1.", file=sys.stderr)
    sys.exit(1)
CITATION.write_text(text, encoding="utf-8", newline="\n")
