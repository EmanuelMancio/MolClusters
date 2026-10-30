# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Label a row of renders (matte.py's output) and write it as light and dark SVG files.

Usage: ``python label_row.py IN.png OUT_PREFIX "title|subtitle" ...``, one label per
object in the image, left to right; objects are found as columns of opaque pixels
separated by empty ones. Writes ``OUT_PREFIX-light.svg`` and ``-dark.svg``, the
image embedded as WebP, the labels in the docs' grey ink.
"""

import base64
import io
import pathlib
import sys
from xml.sax.saxutils import escape

import numpy as np
from PIL import Image

WIDTH = 760  # display width, px
PAD = 48  # room at the sides for the labels
INK = {"light": ("#0b0b0b", "#52514e"), "dark": ("#f0efec", "#c3c2b7")}
FONT = "Roboto, 'Segoe UI', system-ui, -apple-system, sans-serif"


def main() -> None:
    """Write the two SVG files."""
    src, prefix, *labels = sys.argv[1:]
    img = Image.open(src)
    filled = (np.asarray(img)[..., 3] > 8).any(axis=0)
    edges = np.flatnonzero(np.diff(filled.astype(int)))
    starts, ends = edges[::2] + 1, edges[1::2]
    if filled[0]:
        starts = np.r_[0, starts]
    centers = (starts[: len(ends)] + ends) / 2
    if len(centers) != len(labels):
        sys.exit(f"found {len(centers)} objects for {len(labels)} labels")
    scale = (WIDTH - 2 * PAD) / img.width
    h = img.height * scale
    buf = io.BytesIO()
    img.save(buf, "WEBP", quality=90, method=6)
    data = base64.b64encode(buf.getvalue()).decode()
    for theme, (ink, ink2) in INK.items():
        items = [
            f'<image href="data:image/webp;base64,{data}" x="{PAD}" y="0" '
            f'width="{WIDTH - 2 * PAD}" height="{h:.1f}"/>'
        ]
        for x, label in zip(PAD + centers * scale, labels, strict=True):
            title, _, sub = label.partition("|")
            items.append(
                f'<text x="{x:.1f}" y="{h + 18:.1f}" font-size="13" font-weight="bold" '
                f'fill="{ink}" text-anchor="middle">{escape(title)}</text>'
            )
            items.append(
                f'<text x="{x:.1f}" y="{h + 35:.1f}" font-size="12" fill="{ink2}" '
                f'text-anchor="middle">{escape(sub)}</text>'
            )
        body = "\n  ".join(items)
        with pathlib.Path(f"{prefix}-{theme}.svg").open(
            "w", encoding="utf-8", newline="\n"
        ) as f:
            f.write(
                f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" '
                f'height="{h + 44:.0f}" viewBox="0 0 {WIDTH} {h + 44:.0f}" '
                f'font-family="{FONT}">\n  {body}\n</svg>\n'
            )


if __name__ == "__main__":
    main()
