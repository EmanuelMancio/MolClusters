# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Shared fixtures: small synthetic systems whose clustering is known by construction."""

from collections.abc import Callable, Sequence

import MDAnalysis as mda
import numpy as np
import pytest
from loguru import logger
from MDAnalysis.coordinates.memory import MemoryReader

# Residues in the same group sit on a line, BOND_STEP apart (< CUTOFF, so each one
# connects to its neighbours); groups and loose residues sit GROUP_STEP apart, far
# beyond CUTOFF, so nothing connects across groups. The box is large enough that
# periodic images never come into play.
CUTOFF = 3.0
BOND_STEP = 2.0
GROUP_STEP = 50.0
BOX = 2000.0

MOL_RULES = {"MOL": {"MOL": f"cm {CUTOFF}"}}
ALL_PAIRS_RULES = {
    "MOL": {"MOL": f"cm {CUTOFF}", "SOL": f"cm {CUTOFF}"},
    "SOL": {"SOL": f"cm {CUTOFF}"},
}

type Groups = Sequence[Sequence[int]]
type UniverseFactory = Callable[..., mda.Universe]


@pytest.fixture
def captured_logs():
    # The package disables its logger by default (see `molclusters/__init__.py`);
    # `start_logging()` re-enables it for a real run. A dedicated sink is used
    # instead of capsys/capfd, since loguru's default handler binds its own
    # stderr reference ahead of pytest's output capture.
    messages: list[str] = []
    logger.enable("molclusters")
    handler_id = logger.add(messages.append, format="{message}")
    yield messages
    logger.remove(handler_id)
    logger.disable("molclusters")


def _frame_positions(groups: Groups, n_res: int) -> np.ndarray:
    """Place residue centers so that exactly `groups` (1-based resids) are connected.

    Returns
    -------
    np.ndarray
        Array of shape ``(n_res, 3)`` with each residue's reference position.
    """
    positions = np.zeros((n_res, 3))
    placed = set()
    slot = 0

    for group in groups:
        anchor = np.array([100.0 + slot * GROUP_STEP, 100.0, 100.0])
        for k, resid in enumerate(group):
            positions[resid - 1] = anchor + [0.0, k * BOND_STEP, 0.0]
            placed.add(resid)
        slot += 1

    for resid in range(1, n_res + 1):
        if resid not in placed:
            positions[resid - 1] = [100.0 + slot * GROUP_STEP, 100.0, 100.0]
            slot += 1

    return positions


@pytest.fixture
def make_universe() -> UniverseFactory:
    """Build an in-memory Universe of two-atom (C-O) residues from frame layouts.

    Each frame is given as a list of groups of 1-based resids that must form one
    connected cluster under a ``cm CUTOFF`` rule in that frame; any resid not listed
    is left isolated.

    Returns
    -------
    UniverseFactory
        ``make_universe(frames, n_res, resnames=None, blank_elements=())``.
    """

    def factory(
        frames: Sequence[Groups],
        n_res: int,
        resnames: Sequence[str] | None = None,
        blank_elements: Sequence[int] = (),
    ) -> mda.Universe:
        n_atoms = 2 * n_res
        uni = mda.Universe.empty(
            n_atoms,
            n_residues=n_res,
            atom_resindex=np.repeat(np.arange(n_res), 2),
            trajectory=True,
        )

        elements = ["C", "O"] * n_res
        for idx in blank_elements:
            elements[idx] = ""

        uni.add_TopologyAttr("names", ["C1", "O1"] * n_res)
        uni.add_TopologyAttr("types", ["C", "O"] * n_res)
        uni.add_TopologyAttr("elements", elements)
        uni.add_TopologyAttr("resnames", list(resnames or ["MOL"] * n_res))
        uni.add_TopologyAttr("resids", np.arange(1, n_res + 1))
        uni.add_TopologyAttr("masses", [12.011, 15.999] * n_res)
        uni.add_TopologyAttr("charges", [0.3, -0.3] * n_res)
        uni.add_TopologyAttr("bonds", [(2 * i, 2 * i + 1) for i in range(n_res)])

        coords = np.empty((len(frames), n_atoms, 3), dtype=np.float32)
        for f, groups in enumerate(frames):
            centers = _frame_positions(groups, n_res)
            coords[f, 0::2] = centers
            coords[f, 1::2] = centers + [0.0, 0.0, 1.2]

        uni.load_new(
            coords,
            format=MemoryReader,
            dimensions=[BOX, BOX, BOX, 90.0, 90.0, 90.0],
            dt=1.0,
        )
        # JsonReport records both paths
        uni.filename = "synthetic.top"
        uni.trajectory.filename = "synthetic.traj"
        return uni

    return factory
