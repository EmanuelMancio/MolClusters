# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides `BUILTINS`, the built-in analyses, in the order they run.

Classes:
--------
- Builtin: A built-in analysis, how the config builds it and what it needs.

Functions:
----------
- build_builtins: Build the built-in analyses a config runs.

A built-in runs unless the config's `analyses` turns it off (by class name),
and only when the config options it needs are set. Adding one takes a line in
`BUILTINS`; the config checks `analyses` against it.
"""

from collections.abc import Callable
from dataclasses import dataclass

from ..config import MolClsConfig
from .base import FrameAnalysis
from .coordinates import ClusterCoordinates
from .nucleus import Nucleus
from .report import JsonReport
from .size import SizeEvolution
from .solute import SoluteSolvent


@dataclass(frozen=True, slots=True)
class Builtin:
    """A built-in analysis, how the config builds it and what it needs.

    Attributes
    ----------
    kind : type[FrameAnalysis]
        The analysis class; its name turns it on or off in the config's `analyses`.
    build : Callable[[MolClsConfig], FrameAnalysis]
        Builds the analysis from the config, once the options it needs are set.
    needs : tuple[str, ...]
        The config options it needs, without which it doesn't run.
    """

    kind: type[FrameAnalysis]
    build: Callable[[MolClsConfig], FrameAnalysis]
    needs: tuple[str, ...] = ()

    @property
    def name(self) -> str:
        """The analysis' class name, as the config's `analyses` and the log name it.

        Returns
        -------
        str
            E.g. ``"JsonReport"``.
        """
        return self.kind.__name__

    def runs(self, config: MolClsConfig) -> bool:
        """Whether the config runs this analysis.

        Parameters
        ----------
        config : MolClsConfig
            The analysis configuration.

        Returns
        -------
        bool
            True unless `analyses` turns it off or an option it needs isn't set.
        """
        return config.analyses.get(self.name, True) and all(
            getattr(config, option) is not None for option in self.needs
        )


# `JsonReport` records the nuclei of `Nucleus` when that one runs, so it comes after it
BUILTINS: tuple[Builtin, ...] = (
    Builtin(SizeEvolution, lambda c: SizeEvolution()),
    Builtin(
        SoluteSolvent, lambda c: SoluteSolvent(c.solute, c.solvent), needs=("solute",)
    ),
    Builtin(
        ClusterCoordinates,
        lambda c: ClusterCoordinates(c.solute, follow=c._follow_solute),
        needs=("solute",),
    ),
    Builtin(Nucleus, lambda c: Nucleus(c.nucleus), needs=("nucleus",)),
    Builtin(JsonReport, lambda c: JsonReport()),
)


def build_builtins(config: MolClsConfig) -> list[FrameAnalysis]:
    """Build the built-in analyses a config runs, in the order they run.

    Parameters
    ----------
    config : MolClsConfig
        The analysis configuration, with the topology-dependent defaults filled in.

    Returns
    -------
    list[FrameAnalysis]
        The analyses of `BUILTINS` the config runs (see `Builtin.runs`).
    """
    return [builtin.build(config) for builtin in BUILTINS if builtin.runs(config)]
