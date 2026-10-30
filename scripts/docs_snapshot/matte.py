# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Recover a transparent image from two renders of one scene, on black and on white.

Usage: ``python matte.py BLACK WHITE OUT [SCALE]``. On black a pixel is
c = a F, on white c' = a F + (1 - a), so a = 1 - (c' - c) and F = c / a: exact,
anti-aliased edges included. The result is downsampled SCALE times (2) in
premultiplied colour, cropped to what's drawn, and saved in the format OUT names
(e.g. WebP, with its alpha channel).
"""

import sys

import numpy as np
from PIL import Image


def main() -> None:
    """Write the transparent image."""
    black, white, out = sys.argv[1:4]
    scale = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    b = np.asarray(Image.open(black).convert("RGB"), dtype=np.float64) / 255
    w = np.asarray(Image.open(white).convert("RGB"), dtype=np.float64) / 255
    alpha = np.clip(1 - (w - b).mean(axis=2), 0, 1)[..., None]
    pre = np.dstack([b, alpha])  # premultiplied: colour times alpha is the black render
    h, wd = pre.shape[0] // scale, pre.shape[1] // scale
    pre = (
        pre[: h * scale, : wd * scale].reshape(h, scale, wd, scale, 4).mean(axis=(1, 3))
    )
    a = pre[..., 3:]
    rgb = np.where(a > 1e-4, pre[..., :3] / np.maximum(a, 1e-4), 0)
    rgba = np.dstack([np.clip(rgb, 0, 1), a])
    img = Image.fromarray((rgba * 255 + 0.5).astype(np.uint8), "RGBA")
    ys, xs = np.nonzero(rgba[..., 3] > 0.01)
    m = 8
    img = img.crop(
        (max(xs.min() - m, 0), max(ys.min() - m, 0), xs.max() + m + 1, ys.max() + m + 1)
    )
    img.save(out, quality=90, method=6)


if __name__ == "__main__":
    main()
