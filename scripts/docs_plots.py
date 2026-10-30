# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Plot the docs' figures of a run's results, as light and dark SVG files.

Usage: ``uv run --with matplotlib python scripts/docs_plots.py RESULTS``, with
RESULTS the output directory of the quick start's run on ``tests/data/met-mal``.
Writes ``<name>-light.svg`` and ``-dark.svg`` into ``docs/assets/quickstart/`` (the
quick start's ``evo`` plot) and ``docs/assets/outputs/`` (the Output files page's).
The pages show the plain matplotlib code that makes each plot; this only adds the
docs' palette and theme, on transparent backgrounds.
"""

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.use("svg")
ASSETS = Path(__file__).resolve().parent.parent / "docs" / "assets"

THEMES = {
    "light": {
        "series": ["#2a78d6", "#eb6834"],
        "title": "#0b0b0b",
        "text": "#52514e",
        "grid": "#e1e0d9",
        "axis": "#c3c2b7",
    },
    "dark": {
        "series": ["#3987e5", "#d95926"],
        "title": "#ffffff",
        "text": "#c3c2b7",
        "grid": "#2c2c2a",
        "axis": "#383835",
    },
}


def style(theme: str) -> None:
    """Set matplotlib's defaults to a theme's colours."""
    t = THEMES[theme]
    mpl.rcParams.update(
        {
            "font.family": ["Segoe UI", "DejaVu Sans"],
            "font.size": 10,
            "text.color": t["text"],
            "axes.labelcolor": t["text"],
            "axes.edgecolor": t["axis"],
            "axes.facecolor": "none",
            "axes.grid": True,
            "axes.grid.axis": "y",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.spines.left": False,
            "axes.titlecolor": t["text"],
            "axes.prop_cycle": mpl.cycler(color=t["series"]),
            "grid.color": t["grid"],
            "xtick.color": t["text"],
            "ytick.color": t["text"],
            "ytick.left": False,
            "figure.facecolor": "none",
            "savefig.facecolor": "none",
            "svg.hashsalt": "molclusters",  # stable ids, so unchanged plots stay unchanged
        }
    )


def lifetimes(results: Path, theme: str) -> plt.Figure:
    """Plot how many clusters lived how long, uncensored ones only.

    Returns
    -------
    matplotlib.figure.Figure
        The figure.
    """
    lives = pd.read_csv(results / "cluster_lifetimes.csv")
    whole = lives[~lives.BornAtStart & ~lives.AliveAtEnd]
    counts = whole.Lifetime.value_counts().sort_index()
    dt = counts.index.min()
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    ax.bar(counts.index / 1000, counts, width=0.8 * dt / 1000)
    ax.set_yscale("log")
    ax.set_xlabel("Lifetime (ns)")
    ax.set_title("Clusters born and ended within the run, by lifetime", loc="left")
    ax.set_xlim(0, None)
    fig.tight_layout()
    return fig


def solute_solvent(results: Path, theme: str) -> plt.Figure:
    """Plot the mixed clusters' average solute and solvent counts over time.

    Returns
    -------
    matplotlib.figure.Figure
        The figure.
    """
    data = pd.read_csv(results / "solute_solvent.csv")
    mean = data.rolling(21, center=True).mean()  # 1 ns at 48 ps per frame
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    for column, label in (
        ("NSolt", "malic acid (solute)"),
        ("NSolv", "2-methylerythritol"),
    ):
        (line,) = ax.plot(data.Time / 1000, mean[column], linewidth=2)
        last = mean[column].dropna().iloc[-1]
        ax.annotate(
            label,
            (data.Time.iloc[-1] / 1000, last),
            xytext=(6, 0),
            textcoords="offset points",
            va="center",
            color=THEMES[theme]["text"],
        )
    ax.set_xlabel("Time (ns)")
    ax.set_title("Molecules per solute-solvent cluster (1 ns running mean)", loc="left")
    ax.set_ylim(0, None)
    fig.tight_layout()
    return fig


def evo(results: Path, theme: str) -> plt.Figure:
    """Plot the number of clusters and the largest one's size over time.

    Returns
    -------
    matplotlib.figure.Figure
        The figure.
    """
    t = THEMES[theme]
    time, n_clusters, _, _, max_size = np.loadtxt(results / "evo.txt", unpack=True)
    fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(7.2, 4.6))
    # the layout of the first version of this plot, in points of its 518.4 x 331.2
    fig.subplots_adjust(
        left=25.82 / 518.4,
        right=507.6 / 518.4,
        top=1 - 24.808301 / 331.2,
        bottom=1 - 294.130469 / 331.2,
        hspace=(177.34585 - 141.59292) / 116.784619,
    )
    for ax, y, color, title in (
        (top, n_clusters, t["series"][0], "Number of clusters"),
        (bottom, max_size, t["series"][1], "Largest cluster (molecules)"),
    ):
        ax.plot(time / 1000, y, color=color, lw=2, solid_capstyle="round")
        ax.set_title(title, loc="left", fontsize=11, color=t["title"])
        ax.margins(x=0)
        ax.set_ylim(0, y.max() * 1.05)
        ax.set_axisbelow(True)
        ax.tick_params(length=0)
        ax.spines["bottom"].set_linewidth(0.8)
    bottom.set_xlabel("Time (ns)")
    return fig


PLOTS = (
    ("quickstart", "evo", evo),
    ("outputs", "lifetimes", lifetimes),
    ("outputs", "solute-solvent", solute_solvent),
)


def main() -> None:
    """Write every plot in both themes."""
    results = Path(sys.argv[1])
    for folder, name, plot in PLOTS:
        (ASSETS / folder).mkdir(parents=True, exist_ok=True)
        for theme in THEMES:
            style(theme)
            fig = plot(results, theme)
            fig.savefig(
                ASSETS / folder / f"{name}-{theme}.svg", metadata={"Date": None}
            )
            plt.close(fig)


if __name__ == "__main__":
    main()
