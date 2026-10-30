# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Draw the explanatory figures of the docs, as light and dark SVG files.

Run from anywhere: ``uv run python scripts/docs_figures.py``. Writes into
``docs/assets/concepts/``, one ``<name>-light.svg`` and ``<name>-dark.svg`` per
figure; the pages show each with Material's ``#only-light``/``#only-dark``. The
figures are schematic (hand-placed molecules, made-up frames), so they depend on
nothing but this file. Colours are the palette of the quick start's plot: the
first categorical slots for cluster ids, grey ink for text and hairlines, on a
transparent background.
"""

import math
from collections.abc import Callable
from pathlib import Path
from xml.sax.saxutils import escape

OUT = Path(__file__).resolve().parent.parent / "docs" / "assets" / "concepts"

THEMES = {
    "light": {
        "ink": "#0b0b0b",
        "ink2": "#52514e",
        "muted": "#898781",
        "hair": "#e1e0d9",
        "axis": "#c3c2b7",
        "ids": ["#2a78d6", "#eb6834", "#1baf7a"],
    },
    "dark": {
        "ink": "#f0efec",
        "ink2": "#c3c2b7",
        "muted": "#898781",
        "hair": "#2c2c2a",
        "axis": "#383835",
        "ids": ["#3987e5", "#d95926", "#199e70"],
    },
}

FONT = "Roboto, 'Segoe UI', system-ui, -apple-system, sans-serif"
R = 7  # molecule radius in the events figure

Molecule = tuple[float, float, int]  # x, y, cluster at frame t (0: free)
Frame = tuple[list[Molecule], list[tuple[list[int], int, str]]]
Flow = tuple[str, str, int]  # from, to, molecules
FlowStyle = Callable[[str, str, int], tuple[str | None, float, str | None]]


class Svg:
    """A list of SVG elements, written out with a fixed size."""

    def __init__(self, width: float, height: float, theme: str) -> None:
        self.width, self.height = width, height
        self.t = THEMES[theme]
        self.items: list[str] = []

    def add(self, item: str) -> None:
        """Append an SVG element."""
        self.items.append(item)

    def text(
        self,
        x: float,
        y: float,
        s: str,
        size: float = 13,
        color: str = "ink2",
        anchor: str = "start",
        weight: str = "normal",
    ) -> None:
        """Write text at (x, y); colours are theme keys or hex values."""
        fill = self.t.get(color, color)
        self.add(
            f'<text x="{x:g}" y="{y:g}" font-size="{size}" fill="{fill}" '
            f'text-anchor="{anchor}" font-weight="{weight}">{escape(s)}</text>'
        )

    def line(
        self,
        x1: float,
        y1: float,
        x2: float,
        y2: float,
        color: str = "axis",
        width: float = 1.5,
        dash: str | None = None,
        opacity: float = 1,
    ) -> None:
        """Draw a straight line with round caps."""
        stroke = self.t.get(color, color)
        extra = f' stroke-dasharray="{dash}"' if dash else ""
        if opacity != 1:
            extra += f' stroke-opacity="{opacity}"'
        self.add(
            f'<line x1="{x1:g}" y1="{y1:g}" x2="{x2:g}" y2="{y2:g}" stroke="{stroke}" '
            f'stroke-width="{width}" stroke-linecap="round"{extra}/>'
        )

    def circle(
        self,
        x: float,
        y: float,
        r: float,
        fill: str | None = None,
        stroke: str | None = None,
        width: float = 1.5,
    ) -> None:
        """Draw a circle, filled, outlined or both."""
        f = self.t.get(fill, fill) if fill else "none"
        s = (
            f' stroke="{self.t.get(stroke, stroke)}" stroke-width="{width}"'
            if stroke
            else ""
        )
        self.add(f'<circle cx="{x:g}" cy="{y:g}" r="{r:g}" fill="{f}"{s}/>')

    def rect(
        self,
        x: float,
        y: float,
        w: float,
        h: float,
        rx: float = 6,
        stroke: str = "hair",
        fill: str | None = None,
        fill_opacity: float = 1,
        width: float = 1,
        dash: str | None = None,
    ) -> None:
        """Draw a rounded rectangle."""
        f = self.t.get(fill, fill) if fill else "none"
        extra = f' fill-opacity="{fill_opacity}"' if fill and fill_opacity != 1 else ""
        if dash:
            extra += f' stroke-dasharray="{dash}"'
        self.add(
            f'<rect x="{x:g}" y="{y:g}" width="{w:g}" height="{h:g}" rx="{rx:g}" fill="{f}" '
            f'stroke="{self.t.get(stroke, stroke)}" stroke-width="{width}"{extra}/>'
        )

    def path(
        self,
        d: str,
        color: str,
        width: float,
        opacity: float = 1,
        dash: str | None = None,
    ) -> None:
        """Draw an unfilled path."""
        extra = f' stroke-opacity="{opacity}"' if opacity != 1 else ""
        if dash:
            extra += f' stroke-dasharray="{dash}"'
        self.add(
            f'<path d="{d}" fill="none" stroke="{self.t.get(color, color)}" '
            f'stroke-width="{width:g}"{extra}/>'
        )

    def atom(self, x: float, y: float, r: float, label: str = "") -> None:
        """Draw an atom: a translucent grey disc, with an optional letter."""
        c = self.t["muted"]
        self.add(
            f'<circle cx="{x:g}" cy="{y:g}" r="{r:g}" fill="{c}" fill-opacity="0.3" '
            f'stroke="{c}" stroke-width="1.2"/>'
        )
        if label:
            self.text(x, y + 4.5, label, color="ink", anchor="middle", weight="bold")

    def ellipse(self, x: float, y: float, rx: float, ry: float, color: str) -> None:
        """Draw an outlined, lightly filled ellipse."""
        c = self.t.get(color, color)
        self.add(
            f'<ellipse cx="{x:g}" cy="{y:g}" rx="{rx:g}" ry="{ry:g}" fill="{c}" '
            f'fill-opacity="0.15" stroke="{c}" stroke-width="1.5"/>'
        )

    def arrow(self, x1: float, y: float, x2: float, color: str = "muted") -> None:
        """Draw a horizontal arrow from x1 to x2."""
        self.line(x1, y, x2 - 2, y, color, 1.5)
        c = self.t[color]
        self.add(f'<path d="M{x2:g},{y:g} l-7,-4 v8 z" fill="{c}"/>')

    def write(self, path: Path) -> None:
        """Write the SVG file."""
        body = "\n  ".join(self.items)
        path.write_text(
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.width}" '
            f'height="{self.height}" viewBox="0 0 {self.width} {self.height}" '
            f'font-family="{FONT}">\n  {body}\n</svg>\n',
            encoding="utf-8",
            newline="\n",
        )


def idcolor(svg: Svg, cid: int) -> str | None:
    """Give the colour of a cluster id.

    Returns
    -------
    str or None
        The colour of id 1, 2 or 3; None for 0 (free).
    """
    return svg.t["ids"][cid - 1] if cid else None


# --- Events ----------------------------------------------------------------------

# Each event: title, caption, then frame t and frame t+1. A frame is a list of
# molecules (x, y, cluster at t; 0 = free) and of clusters (molecule indices, id, label).
EVENTS: list[tuple[str, str, Frame, Frame]] = [
    (
        "formation",
        "Free molecules come together: a new id.",
        ([(30, 30, 0), (72, 48, 0), (42, 82, 0), (118, 76, 0)], []),
        (
            [(50, 42, 0), (70, 42, 0), (60, 59, 0), (118, 80, 0)],
            [([0, 1, 2], 3, "3 new")],
        ),
    ),
    (
        "growth",
        "Molecules join: the cluster keeps its id.",
        (
            [(38, 50, 1), (58, 50, 1), (48, 67, 1), (112, 34, 0), (116, 86, 0)],
            [([0, 1, 2], 1, "1")],
        ),
        (
            [(38, 50, 1), (58, 50, 1), (48, 67, 1), (68, 67, 0), (78, 50, 0)],
            [([0, 1, 2, 3, 4], 1, "1")],
        ),
    ),
    (
        "split",
        "The piece with most of its molecules keeps the id.",
        (
            [
                (30, 50, 1),
                (50, 50, 1),
                (40, 67, 1),
                (60, 67, 1),
                (80, 67, 1),
                (100, 67, 1),
            ],
            [([0, 1, 2, 3, 4, 5], 1, "1")],
        ),
        (
            [
                (26, 50, 1),
                (46, 50, 1),
                (36, 67, 1),
                (56, 67, 1),
                (100, 72, 1),
                (120, 72, 1),
            ],
            [([0, 1, 2, 3], 1, "1"), ([4, 5], 3, "3 new")],
        ),
    ),
    (
        "merge",
        "1 gave more molecules and keeps its id; 2 ends.",
        (
            [
                (26, 44, 1),
                (46, 44, 1),
                (36, 61, 1),
                (56, 61, 1),
                (102, 78, 2),
                (122, 78, 2),
            ],
            [([0, 1, 2, 3], 1, "1"), ([4, 5], 2, "2")],
        ),
        (
            [
                (26, 56, 1),
                (46, 56, 1),
                (36, 73, 1),
                (56, 73, 1),
                (76, 73, 2),
                (96, 73, 2),
            ],
            [([0, 1, 2, 3, 4, 5], 1, "1")],
        ),
    ),
    (
        "dissolution",
        "1 leaves single molecules only: it ends.",
        (
            [
                (30, 40, 1),
                (50, 40, 1),
                (40, 57, 1),
                (92, 76, 2),
                (112, 76, 2),
                (102, 93, 2),
            ],
            [([0, 1, 2], 1, "1"), ([3, 4, 5], 2, "2")],
        ),
        (
            [
                (24, 30, 1),
                (40, 62, 1),
                (72, 76, 1),
                (92, 76, 2),
                (112, 76, 2),
                (102, 93, 2),
            ],
            [([2, 3, 4, 5], 2, "2")],
        ),
    ),
    (
        "mixed dimer",
        "One molecule from each of two clusters: always new.",
        (
            [
                (26, 40, 1),
                (46, 40, 1),
                (36, 57, 1),
                (98, 74, 2),
                (118, 74, 2),
                (108, 91, 2),
            ],
            [([0, 1, 2], 1, "1"), ([3, 4, 5], 2, "2")],
        ),
        (
            [
                (22, 36, 1),
                (42, 36, 1),
                (66, 72, 1),
                (86, 72, 2),
                (110, 40, 2),
                (130, 40, 2),
            ],
            [([0, 1], 1, "1"), ([4, 5], 2, "2"), ([2, 3], 3, "3 new")],
        ),
    ),
]

BOX_W, BOX_H = 150, 112


def draw_frame(svg: Svg, ox: float, oy: float, frame: Frame) -> None:
    """Draw one frame: hairline box, then clusters' hulls, connections and molecules."""
    mols, clusters = frame
    svg.rect(ox, oy, BOX_W, BOX_H, rx=8)
    for members, cid, label in clusters:
        xs = [mols[i][0] for i in members]
        ys = [mols[i][1] for i in members]
        pad = R + 5
        x0, y0 = ox + min(xs) - pad, oy + min(ys) - pad
        w, h = max(xs) - min(xs) + 2 * pad, max(ys) - min(ys) + 2 * pad
        color = idcolor(svg, cid)
        svg.rect(
            x0, y0, w, h, rx=pad, stroke=color, fill=color, fill_opacity=0.1, width=1.5
        )
        svg.text(x0 + 2, y0 - 4, label, size=12, color="ink", weight="bold")
        # connections: each molecule to its nearest one among those before it
        for k in range(1, len(members)):
            a = mols[members[k]]
            b = min(
                (mols[i] for i in members[:k]),
                key=lambda m: (m[0] - a[0]) ** 2 + (m[1] - a[1]) ** 2,
            )
            # stop at the circles' edges: free molecules are hollow
            d = ((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2) ** 0.5
            ux, uy = (b[0] - a[0]) / d * R, (b[1] - a[1]) / d * R
            x1, y1, x2, y2 = (
                ox + a[0] + ux,
                oy + a[1] + uy,
                ox + b[0] - ux,
                oy + b[1] - uy,
            )
            svg.line(x1, y1, x2, y2, "muted", 1.5)
    for x, y, cid in mols:
        if cid:
            svg.circle(ox + x, oy + y, R, fill=idcolor(svg, cid))
        else:
            svg.circle(ox + x, oy + y, R - 0.75, stroke="muted")


def events(theme: str) -> Svg:
    """Draw six events, each as frame t and frame t + 1.

    Returns
    -------
    Svg
        The figure.
    """
    pw, ph, gx, gy, top = 350, 184, 30, 18, 34
    svg = Svg(2 * pw + gx, top + 3 * ph + 2 * gy, theme)
    # key
    x = 0
    svg.circle(x + 7, 12, R, fill=idcolor(svg, 1))
    svg.circle(x + 25, 12, R, fill=idcolor(svg, 2))
    svg.text(x + 40, 16, "in cluster 1 or 2 at t (molecules keep this colour)", size=12)
    x = 330
    svg.circle(x + 7, 12, R - 0.75, stroke="muted")
    svg.text(x + 22, 16, "free at t", size=12)
    x = 420
    svg.rect(
        x,
        3,
        30,
        18,
        rx=9,
        stroke=idcolor(svg, 3),
        fill=idcolor(svg, 3),
        fill_opacity=0.1,
        width=1.5,
    )
    svg.text(x + 38, 16, "a cluster in that frame, with its id", size=12)
    for n, (title, caption, before, after) in enumerate(EVENTS):
        ox = (n % 2) * (pw + gx)
        oy = top + (n // 2) * (ph + gy)
        svg.text(ox, oy + 14, title, size=14, color="ink", weight="bold")
        draw_frame(svg, ox, oy + 24, before)
        svg.arrow(ox + BOX_W + 8, oy + 24 + BOX_H / 2, ox + pw - BOX_W - 8)
        draw_frame(svg, ox + pw - BOX_W, oy + 24, after)
        svg.text(
            ox + BOX_W / 2,
            oy + 24 + BOX_H + 15,
            "t",
            size=12,
            color="muted",
            anchor="middle",
        )
        svg.text(
            ox + pw - BOX_W / 2,
            oy + 24 + BOX_H + 15,
            "t + 1",
            size=12,
            color="muted",
            anchor="middle",
        )
        svg.text(ox, oy + 24 + BOX_H + 35, caption, size=13)
    return svg


# --- Id assignment ---------------------------------------------------------------

# Frame t: clusters 7, 9, 12 and the free molecules; frame t+1: groups A, B, C and
# the free molecules. Flows: (from, to, molecules). Group A is tracking.md's
# Counter({7: 4, 9: 2, 0: 1}).
LEFT = [("7", 1, 6), ("9", 2, 4), ("12", 3, 2), ("free", 0, 2)]
RIGHT = ["A", "B", "C", "free"]
FLOWS: list[Flow] = [
    ("7", "A", 4),
    ("7", "B", 2),
    ("9", "A", 2),
    ("9", "C", 1),
    ("9", "free", 1),
    ("12", "C", 1),
    ("12", "free", 1),
    ("free", "A", 1),
    ("free", "C", 1),
]
PICKS = {("7", "A"), ("9", "A")}  # each cluster's best group (12 picks nothing)
COLOR_OF = {name: c for name, c, _ in LEFT}

LX, RX, IBW, IBH, UNIT = 150, 450, 124, 30, 3.2  # box positions, size, px per molecule
ROW0, ROWSTEP = 62, 44


def dots(svg: Svg, x: float, y: float, origins: list[int]) -> None:
    """Draw a row of molecules inside a box, coloured by their cluster at t."""
    for k, c in enumerate(origins):
        cx = x + 14 + 16 * k
        if c:
            svg.circle(cx, y, 6, fill=idcolor(svg, c))
        else:
            svg.circle(cx, y, 5.25, stroke="muted")


def flow_ports(
    flows: list[Flow],
) -> tuple[dict[tuple[str, str], float], dict[tuple[str, str], float]]:
    """Find where each flow leaves its left box and enters its right box.

    Returns
    -------
    tuple of dict
        The y offsets, by (from, to), on the left and on the right boxes.
    """
    rows_l = {name: i for i, (name, _, _) in enumerate(LEFT)}
    rows_r = {name: i for i, name in enumerate(RIGHT)}
    out: dict[tuple[str, str], float] = {}
    inn: dict[tuple[str, str], float] = {}
    for side, key, order, rows in ((out, 0, 1, rows_l), (inn, 1, 0, rows_r)):
        for name in rows:
            mine = sorted(
                (f for f in flows if f[key] == name),
                key=lambda f: rows_l[f[0]] if order == 0 else rows_r[f[1]],
            )
            total = sum(f[2] for f in mine) * UNIT
            y = ROW0 + rows[name] * ROWSTEP - total / 2
            for f in mine:
                side[f[:2]] = y + f[2] * UNIT / 2
                y += f[2] * UNIT
    return out, inn


def draw_flows(svg: Svg, oy: float, flows: list[Flow], style: FlowStyle) -> None:
    """Draw the flows as curves as wide as the molecules they carry."""
    out, inn = flow_ports(flows)
    for a, b, n in flows:
        y1, y2 = oy + out[(a, b)], oy + inn[(a, b)]
        x1, x2 = LX + IBW + 2, RX - 2
        xm = (x1 + x2) / 2
        d = f"M{x1:g},{y1:g} C{xm:g},{y1:g} {xm:g},{y2:g} {x2:g},{y2:g}"
        color, opacity, dash = style(a, b, n)
        svg.path(d, color or svg.t["muted"], n * UNIT - 0.6, opacity, dash)


def id_assignment(theme: str) -> Svg:
    """Draw the two steps of id assignment on one pair of frames.

    Returns
    -------
    Svg
        The figure.
    """
    width, ph = 760, 236
    svg = Svg(width, 2 * ph + 20, theme)
    content = {
        "A": [1, 1, 1, 1, 2, 2, 0],
        "B": [1, 1],
        "C": [2, 3, 0],
        "free": [2, 3],
    }
    for panel in (0, 1):
        oy = panel * (ph + 20)
        title = (
            "1. Each cluster of frame t picks the group holding most of its molecules"
            if panel == 0
            else "2. Each group continues the cluster that gave it most molecules, among those that picked it"
        )
        svg.text(0, oy + 16, title, size=14, color="ink", weight="bold")
        svg.text(
            LX + IBW / 2, oy + 38, "frame t", size=12, color="muted", anchor="middle"
        )
        svg.text(
            RX + IBW / 2,
            oy + 38,
            "frame t + 1",
            size=12,
            color="muted",
            anchor="middle",
        )
        if panel == 0:

            def style(a: str, b: str, n: int) -> tuple[str | None, float, str | None]:
                color = idcolor(svg, COLOR_OF[a])
                if n == 1:
                    return color, 0.45, "3 3"  # a single molecule never carries an id
                return color, 1 if (a, b) in PICKS else 0.4, None

            draw_flows(svg, oy, FLOWS, style)
            notes_l = ["picks A (4 of 6)", "picks A (2 of 4)", "picks nothing", ""]
            notes_r = ["picked by 7 and 9", "picked by nobody", "picked by nobody", ""]
        else:

            def style(a: str, b: str, n: int) -> tuple[str | None, float, str | None]:  # noqa: ARG001
                return idcolor(svg, COLOR_OF[a]), 1 if a == "7" else 0.4, None

            draw_flows(svg, oy, [f for f in FLOWS if f[:2] in PICKS], style)
            notes_l = ["continues in A", "ends: merged into 7", "ends: dissolved", ""]
            notes_r = ["continues 7 (4 beat 2)", "new id 13", "new id 14", ""]
        for i, (name, cid, n) in enumerate(LEFT):
            y = oy + ROW0 + i * ROWSTEP
            svg.rect(LX, y - IBH / 2, IBW, IBH, rx=8, stroke="axis", width=1)
            dots(svg, LX, y, [cid] * n)
            label = name if name == "free" else f"cluster {name}"
            svg.text(
                LX - 10,
                y - 2 if notes_l[i] else y + 4,
                label,
                size=13,
                color="ink",
                anchor="end",
                weight="bold",
            )
            if notes_l[i]:
                svg.text(LX - 10, y + 13, notes_l[i], size=12, anchor="end")
        for i, name in enumerate(RIGHT):
            y = oy + ROW0 + i * ROWSTEP
            stroke = "axis"
            if panel == 1 and name == "A":
                stroke = idcolor(svg, 1)
            svg.rect(
                RX,
                y - IBH / 2,
                IBW,
                IBH,
                rx=8,
                stroke=stroke,
                width=1.5 if stroke != "axis" else 1,
            )
            dots(svg, RX, y, content[name])
            label = name if name == "free" else f"group {name}"
            svg.text(
                RX + IBW + 10,
                y - 2 if notes_r[i] else y + 4,
                label,
                size=13,
                color="ink",
                weight="bold",
            )
            if notes_r[i]:
                svg.text(RX + IBW + 10, y + 13, notes_r[i], size=12)
    return svg


# --- Lifetimes -------------------------------------------------------------------

# Clusters as (id, first frame seen, last frame seen) over frames 0-9.
LIVES = [
    (1, 0, 3, "born at start: ≥ 4Δt"),
    (2, 2, 6, "lifetime 5Δt"),
    (3, 5, 5, "one frame: lifetime Δt"),
    (4, 6, 9, "alive at end: ≥ 3Δt"),
]
NFRAMES = 10


def lifetimes(theme: str) -> Svg:
    """Draw four clusters' sightings, birth and death times and lifetimes.

    Returns
    -------
    Svg
        The figure.
    """
    width, x0, dx, top, rowh = 760, 150, 46, 62, 46
    xf = [x0 + k * dx for k in range(NFRAMES)]
    last = xf[-1]
    height = top + len(LIVES) * rowh + 52
    svg = Svg(width, height, theme)
    blue, ink2 = idcolor(svg, 1), svg.t["ink2"]
    # key
    svg.circle(7, 12, 5.5, fill=blue)
    svg.text(20, 16, "seen in the frame", size=12)
    svg.circle(145, 12, 5, stroke="ink2")
    svg.text(158, 16, "first frame gone (DeathTime)", size=12)
    svg.line(340, 12, 370, 12, blue, 3, dash="4 4")
    svg.text(378, 16, "formed or ended somewhere in here", size=12)
    svg.line(590, 8, 590, 16, "ink2", 1.2)
    svg.line(590, 12, 616, 12, "ink2", 1.2)
    svg.line(616, 8, 616, 16, "ink2", 1.2)
    svg.text(624, 16, "recorded Lifetime", size=12)
    # frame gridlines and axis
    axis_y = top + len(LIVES) * rowh - 10
    for k, x in enumerate(xf):
        svg.line(x, top - 22, x, axis_y, "hair", 1)
        svg.text(
            x,
            axis_y + 18,
            "0" if k == 0 else ("Δt" if k == 1 else f"{k}Δt"),
            size=12,
            color="muted",
            anchor="middle",
        )
    svg.line(xf[0] - 30, axis_y, last + 30, axis_y, "axis", 1)
    svg.text(xf[0], axis_y + 36, "first frame", size=12, color="muted", anchor="middle")
    svg.text(last, axis_y + 36, "last frame", size=12, color="muted", anchor="middle")
    for i, (cid, first, lastseen, note) in enumerate(LIVES):
        y = top + i * rowh
        svg.text(
            x0 - 40,
            y + 4,
            f"cluster {cid}",
            size=13,
            color="ink",
            anchor="end",
            weight="bold",
        )
        born_at_start, alive_at_end = first == 0, lastseen == NFRAMES - 1
        # when it formed: before frame 0 (unknown), or between the frame before and the first seen
        lead = xf[0] - 34 if born_at_start else xf[first - 1]
        svg.line(lead, y, xf[first], y, blue, 3, dash="4 4")
        if lastseen > first:
            svg.line(xf[first], y, xf[lastseen], y, blue, 3)
        if alive_at_end:
            svg.line(last, y, last + 34, y, blue, 3, dash="4 4")
            end = last
        else:
            end = xf[lastseen + 1]
            svg.line(xf[lastseen], y, end, y, blue, 3, dash="4 4")
            svg.circle(end, y, 5, stroke="ink2")
        for k in range(first, lastseen + 1):
            svg.circle(xf[k], y, 5.5, fill=blue)
        # the recorded lifetime: BirthTime to DeathTime, or to the last frame
        b = y + 14
        svg.line(xf[first], b - 3, xf[first], b + 3, ink2, 1.2)
        svg.line(xf[first], b, end, b, ink2, 1.2)
        svg.line(end, b - 3, end, b + 3, ink2, 1.2)
        svg.text(last + 50, y + 4, note, size=13)
    return svg


# --- Connection rules ------------------------------------------------------------

MOLECULE = [(0, 0, 8), (15, -7, 7), (11, 13, 6), (-13, 9, 6)]  # atoms: dx, dy, r


def molecule(svg: Svg, x: float, y: float, label: str) -> None:
    """Draw a small molecule with its center of mass marked by a cross."""
    for dx, dy, r in MOLECULE:
        svg.atom(x + dx, y + dy, r)
    svg.line(x - 5, y - 5, x + 5, y + 5, "ink", 1.6)
    svg.line(x - 5, y + 5, x + 5, y - 5, "ink", 1.6)
    svg.text(x, y + 38, label, size=12, color="muted", anchor="middle")


def rule_cm(theme: str) -> Svg:
    """Draw a cm rule: centers of mass within the cutoff connect.

    Returns
    -------
    Svg
        The figure.
    """
    svg = Svg(600, 250, theme)
    blue = idcolor(svg, 1)
    a, b, c, cut = (150, 125), (228, 100), (355, 175), 112
    svg.add(
        f'<circle cx="{a[0]}" cy="{a[1]}" r="{cut}" fill="{blue}" fill-opacity="0.06" '
        f'stroke="{blue}" stroke-width="1.5" stroke-dasharray="5 4"/>'
    )
    svg.line(*a, *b, blue, 2.5)
    svg.line(*a, *c, "muted", 1.5, dash="2 4")
    # the cutoff, as a radius pointing down-left
    ang = math.radians(215)
    ex, ey = a[0] + cut * math.cos(ang), a[1] - cut * math.sin(ang)
    svg.line(*a, ex, ey, "ink2", 1.2)
    svg.text((a[0] + ex) / 2 - 6, (a[1] + ey) / 2 - 8, "cutoff", anchor="end")
    for (x, y), label in ((a, "molecule 1"), (b, "molecule 2"), (c, "molecule 3")):
        molecule(svg, x, y, label)
    svg.text(400, 70, "1 and 2: centers of mass at", color="ink")
    svg.text(400, 88, "most the cutoff apart:", color="ink")
    svg.text(400, 106, "connected", color="ink", weight="bold")
    svg.text(400, 200, "1 and 3: farther apart:")
    svg.text(400, 218, "not connected", weight="bold")
    svg.line(400, 238, 410, 248, "ink", 1.6)
    svg.line(400, 248, 410, 238, "ink", 1.6)
    svg.text(418, 247, "center of mass", size=12)
    return svg


def arc(cx: float, cy: float, r: float, a1: float, a2: float) -> str:
    """Give the SVG path of the short arc from angle a1 to a2 (radians, screen axes).

    Returns
    -------
    str
        The path's ``d`` attribute.
    """
    delta = (a2 - a1 + math.pi) % (2 * math.pi) - math.pi
    x1, y1 = cx + r * math.cos(a1), cy + r * math.sin(a1)
    x2, y2 = cx + r * math.cos(a2), cy + r * math.sin(a2)
    return f"M{x1:g},{y1:g} A{r:g},{r:g} 0 0 {1 if delta > 0 else 0} {x2:g},{y2:g}"


def rule_hb(theme: str) -> Svg:
    """Draw an hb rule's geometry: the D-A distance and the D-H-A angle.

    Returns
    -------
    Svg
        The figure.
    """
    svg = Svg(600, 280, theme)
    blue, orange = idcolor(svg, 1), idcolor(svg, 2)
    d, h, a = (190, 140), (232, 122), (360, 114)
    d2, a2 = (148, 168), (400, 142)  # the atoms D and A are bonded to
    svg.line(*d2, *d, "muted", 2)
    svg.line(*d, *h, "muted", 2)
    svg.line(*a, *a2, "muted", 2)
    svg.line(*h, *a, "ink2", 1.5, dash="2 4")
    # D-A distance, drawn below the atoms
    off = 82
    svg.line(d[0], d[1] + 16, d[0], d[1] + off + 4, "axis", 1)
    svg.line(a[0], a[1] + 16, a[0], a[1] + off + 12, "axis", 1)
    y1, y2 = d[1] + off, a[1] + off + 8
    svg.line(d[0], y1, a[0], y2, blue, 2)
    svg.text(
        (d[0] + a[0]) / 2,
        (y1 + y2) / 2 + 24,
        "donor-acceptor distance: at most d",
        color="ink",
        anchor="middle",
    )
    # D-H-A angle at H
    ad = math.atan2(d[1] - h[1], d[0] - h[0])
    aa = math.atan2(a[1] - h[1], a[0] - h[0])
    svg.path(arc(*h, 24, ad, aa), orange, 2)
    svg.line(h[0] + 4, h[1] + 25, h[0] + 20, h[1] + 44, "axis", 1)
    svg.text(
        h[0] + 24, h[1] + 54, "D-H-A angle: at least a (180° is straight)", color="ink"
    )
    svg.atom(*d2, 8)
    svg.atom(*a2, 8)
    svg.atom(*d, 13, "D")
    svg.atom(*h, 9, "H")
    svg.atom(*a, 13, "A")
    svg.text(d[0] - 22, d[1] - 30, "donor: the atom", size=12, anchor="end")
    svg.text(d[0] - 22, d[1] - 15, "H is bonded to", size=12, anchor="end")
    svg.text(a[0] + 22, a[1] - 18, "acceptor: charge < aq", size=12)
    svg.line(h[0], h[1] - 11, h[0], h[1] - 34, "axis", 1)
    svg.text(h[0], h[1] - 54, "hydrogen: mass between hmin", size=12, anchor="middle")
    svg.text(h[0], h[1] - 39, "and hmax, charge > hq", size=12, anchor="middle")
    svg.text(110, 200, "molecule 1", size=12, color="muted", anchor="middle")
    svg.text(416, 150, "molecule 2", size=12, color="muted")
    return svg


# --- Size ------------------------------------------------------------------------


def size(theme: str) -> Svg:
    """Draw the equivalent sphere: R_g, the sqrt(5/3) R_g sphere, and R.

    Returns
    -------
    Svg
        The figure.
    """
    svg = Svg(640, 300, theme)
    blue = idcolor(svg, 1)
    cx, cy, rs, ra, n = 150, 150, 100, 13, 46
    rg = rs * math.sqrt(3 / 5)
    big = rs + ra / 2
    # atoms spread evenly over a disk (Vogel's spiral), the outer ones' centers at rs
    for k in range(n):
        r = rs * math.sqrt((k + 0.5) / (n - 0.5))
        t = k * math.pi * (3 - math.sqrt(5))
        svg.atom(cx + r * math.cos(t), cy + r * math.sin(t), ra)
    svg.add(
        f'<circle cx="{cx}" cy="{cy}" r="{rg:g}" fill="none" stroke="{svg.t["ink2"]}" '
        f'stroke-width="1.5" stroke-dasharray="5 4"/>'
    )
    svg.circle(cx, cy, rs, stroke="ink2", width=1.5)
    svg.circle(cx, cy, big, stroke=blue, width=2.5)
    svg.circle(cx, cy, 3, fill="ink")
    labels = [
        (-42, big, blue, "R = √(5/3) Rg + ½ r̄vdW: the equivalent sphere"),
        (-14, rs, "ink2", "√(5/3) Rg: as far as the outer atoms' centers"),
        (14, rg, "ink2", "Rg: radius of gyration (dashed)"),
        (40, rs + ra, "ink2", "atoms, at their van der Waals radius"),
    ]
    for deg, r, color, text in labels:
        t = math.radians(deg)
        px, py = cx + r * math.cos(t), cy + r * math.sin(t)
        ty = cy + 108 * math.sin(t)
        svg.line(px, py, 300, ty, color, 1)
        svg.circle(px, py, 2.5, fill=color)
        svg.text(306, ty + 4, text, color="ink" if color == blue else "ink2")
    return svg


# --- Shape -----------------------------------------------------------------------


def descriptors(lam: tuple[float, float, float]) -> tuple[float, float]:
    """Give the shape parameter and sphericity of gyration-tensor eigenvalues.

    Returns
    -------
    tuple of float
        S and 1 - asphericity, as defined on the Cluster properties page.
    """
    total = sum(lam)
    dev = [x - total / 3 for x in lam]
    asph = 1.5 * sum(x * x for x in dev) / total**2
    return 27 * dev[0] * dev[1] * dev[2] / total**3, 1 - asph


def shape(theme: str) -> Svg:
    """Draw where shapes fall by shape parameter and sphericity.

    Returns
    -------
    Svg
        The figure.
    """
    svg = Svg(620, 340, theme)
    blue = idcolor(svg, 1)
    x0, x1, y0, y1 = 80, 590, 290, 30  # plot area: S -0.25..2, sphericity 0..1

    def px(sv: float) -> float:
        return x0 + (sv + 0.25) / 2.25 * (x1 - x0)

    def py(v: float) -> float:
        return y0 + v * (y1 - y0)

    for v in (0, 0.25, 0.5, 0.75, 1):
        svg.line(x0, py(v), x1, py(v), "hair", 1)
        svg.text(x0 - 8, py(v) + 4, f"{v:g}", size=12, color="muted", anchor="end")
    for sv in (-0.25, 0, 0.5, 1, 1.5, 2):
        svg.line(px(sv), y0, px(sv), y0 + 5, "axis", 1)
        label = f"{sv:g}".replace("-", "−")
        svg.text(px(sv), y0 + 20, label, size=12, color="muted", anchor="middle")
    svg.line(x0, y0, x1, y0, "axis", 1)
    svg.text((x0 + x1) / 2, y0 + 42, "shape parameter S", anchor="middle")
    ym = (y0 + y1) / 2
    svg.add(
        f'<text x="24" y="{ym:g}" font-size="13" fill="{svg.t["ink2"]}" '
        f'text-anchor="middle" transform="rotate(-90 24 {ym:g})">sphericity</text>'
    )
    # every set of eigenvalues falls between the spheroids and the flat shapes
    steps = [k / 100 for k in range(101)]
    prolate = [descriptors((1, 1 - b, 1 - b)) for b in steps]  # sphere -> rod
    flat = [descriptors((1, b, 0)) for b in steps]  # rod -> disk
    oblate = [descriptors((1, 1, b)) for b in steps]  # disk -> sphere
    pts = " ".join(f"{px(sv):.1f},{py(v):.1f}" for sv, v in prolate + flat + oblate)
    svg.add(
        f'<polygon points="{pts}" fill="{blue}" fill-opacity="0.1" stroke="{blue}" '
        f'stroke-width="1.5" stroke-linejoin="round"/>'
    )
    for lam, label, dx, dy, anchor in (
        ((1, 0.35, 0.35), "prolate", 10, -8, "start"),
        ((1, 0.25, 0), "flat", -4, 22, "end"),
        ((1, 1, 0.4), "oblate", -10, 4, "end"),
    ):
        sv, v = descriptors(lam)
        svg.text(px(sv) + dx, py(v) + dy, label, size=12, anchor=anchor)
    # reference shapes: eigenvalues, label offset, sketch offset and (rx, ry)
    marks = [
        ((1, 1, 1), "sphere, cube", (12, 4), (112, 0), (10, 10)),
        ((1, 1, 0), "disk, ring", (-8, 74), (24, 94), (16, 4)),
        ((1, 0, 0), "rod", (-12, -14), (-24, -38), (20, 2.5)),
        ((1, 0.15, 0.15), "cigar", (12, 4), (76, 0), (16, 6)),
        ((1, 0.5, 0), "A4 sheet", (14, 22), (36, 46), (14, 10)),
    ]
    for lam, label, (dx, dy), (ix, iy), (rx, ry) in marks:
        sv, v = descriptors(lam)
        x, y = px(sv), py(v)
        svg.circle(x, y, 4, fill="ink")
        anchor = "end" if dx < 0 and dy < 40 else "start"
        if dy > 40:  # a label far below its point: draw a leader
            svg.line(x, y + 7, x, y + dy - 16, "axis", 1)
        svg.text(x + dx, y + dy, label, color="ink", weight="bold", anchor=anchor)
        svg.ellipse(x + ix, y + iy, rx, ry, "ink2")
    return svg


# --- Nucleus ---------------------------------------------------------------------

# A cluster of solutes (S) and solvent molecules (V): positions, then connections.
NUC_MOLS = {
    "S1": (70, 130),
    "S2": (102, 106),
    "S3": (102, 154),
    "V1": (142, 130),
    "V2": (182, 130),
    "V5": (184, 172),
    "S4": (222, 130),
    "V3": (262, 106),
    "S5": (302, 106),
    "S6": (334, 130),
    "V4": (298, 158),
}
NUC_EDGES = [
    ("S1", "S2"),
    ("S1", "S3"),
    ("S2", "S3"),
    ("S2", "V1"),
    ("S3", "V1"),
    ("V1", "V2"),
    ("V2", "S4"),
    ("V2", "V5"),
    ("S4", "V3"),
    ("V3", "S5"),
    ("S5", "S6"),
    ("S6", "V4"),
]
NUCLEI = [
    (["S1", "S2", "S3"], "nucleus 1", -1),
    (["S4"], "nucleus 2", 1),
    (["S5", "S6"], "nucleus 3", -1),
]


def nucleus(theme: str) -> Svg:
    """Draw a cluster's nuclei: its graph restricted to the solutes, in pieces.

    Returns
    -------
    Svg
        The figure.
    """
    svg = Svg(560, 230, theme)
    blue, orange = idcolor(svg, 1), idcolor(svg, 2)
    r = 10
    xs = [x for x, _ in NUC_MOLS.values()]
    ys = [y for _, y in NUC_MOLS.values()]
    pad = r + 16
    x0, y0 = min(xs) - pad, min(ys) - pad - 14
    svg.rect(
        x0,
        y0,
        max(xs) - min(xs) + 2 * pad,
        max(ys) - min(ys) + 2 * pad + 14,
        rx=24,
        stroke="axis",
        width=1.5,
    )
    svg.text(x0 + 4, y0 - 8, "cluster", size=12, color="ink", weight="bold")
    for members, label, side in NUCLEI:
        mx = [NUC_MOLS[m][0] for m in members]
        my = [NUC_MOLS[m][1] for m in members]
        p = r + 5
        svg.rect(
            min(mx) - p,
            min(my) - p,
            max(mx) - min(mx) + 2 * p,
            max(my) - min(my) + 2 * p,
            rx=p,
            stroke=orange,
            width=1.5,
            dash="4 3",
        )
        ty = max(my) + p + 15 if side > 0 else min(my) - p - 6
        svg.text(
            (min(mx) + max(mx)) / 2, ty, label, size=12, color="ink", anchor="middle"
        )
    for u, v in NUC_EDGES:
        (x1, y1), (x2, y2) = NUC_MOLS[u], NUC_MOLS[v]
        dist = math.hypot(x2 - x1, y2 - y1)
        ux, uy = (x2 - x1) / dist * r, (y2 - y1) / dist * r
        svg.line(x1 + ux, y1 + uy, x2 - ux, y2 - uy, "muted", 1.5)
    for name, (x, y) in NUC_MOLS.items():
        if name[0] == "S":
            svg.circle(x, y, r, fill=blue)
        else:
            svg.circle(x, y, r - 0.75, stroke="muted")
    kx = 410
    svg.circle(kx + 8, 82, r, fill=blue)
    svg.text(kx + 26, 87, "solute", size=12)
    svg.circle(kx + 8, 112, r - 0.75, stroke="muted")
    svg.text(kx + 26, 117, "solvent", size=12)
    svg.rect(kx - 4, 133, 24, 18, rx=8, stroke=orange, width=1.5, dash="4 3")
    svg.text(kx + 26, 146, "nucleus", size=12)
    svg.line(kx - 2, 172, kx + 18, 172, "muted", 1.5)
    svg.text(kx + 26, 176, "connection", size=12)
    return svg


def main() -> None:
    """Write every figure in both themes."""
    OUT.mkdir(parents=True, exist_ok=True)
    for name, draw in (
        ("events", events),
        ("id-assignment", id_assignment),
        ("lifetimes", lifetimes),
        ("rule-cm", rule_cm),
        ("rule-hb", rule_hb),
        ("size", size),
        ("shape", shape),
        ("nucleus", nucleus),
    ):
        for theme in THEMES:
            draw(theme).write(OUT / f"{name}-{theme}.svg")


if __name__ == "__main__":
    main()
