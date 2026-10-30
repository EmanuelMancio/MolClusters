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


def main() -> None:
    """Write every figure in both themes."""
    OUT.mkdir(parents=True, exist_ok=True)
    for name, draw in (
        ("events", events),
        ("id-assignment", id_assignment),
        ("lifetimes", lifetimes),
    ):
        for theme in THEMES:
            draw(theme).write(OUT / f"{name}-{theme}.svg")


if __name__ == "__main__":
    main()
