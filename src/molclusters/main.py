# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Serves as the main entry point for analyzing molecular clusters using the MolClusters library.

Functions:
----------
- main():
    Sets up the command-line interface, parses arguments, initializes the MDAnalysis Universe,
    fills in missing van der Waals radii, reads the analysis config, and runs the molecular
    cluster analysis.

Dependencies:
-------------
- argparse: For parsing command-line arguments.
- datetime: For timestamping the log file.
- MDAnalysis: For molecular dynamics trajectory and structure analysis.
- loguru: For logging.
- MolClusters: The core library for molecular cluster analysis.

Usage:
------
Run via the installed console script with the required arguments:
    molclusters <trajectory_file> <topology_file> <input_file>

Arguments:
----------
- traj: The trajectory file for the molecular dynamics simulation.
- top: The topology file for the molecular dynamics simulation.
- inp: The input YAML/JSON/TOML file containing analysis settings.

Optional Arguments:
-------------------
- --traj-memory: Load the trajectory into memory.
- --in-memory-step: Step for in-memory trajectory loading (requires --traj-memory).
- --version: Displays the version of the MolClusters library.

Example:
--------
    molclusters trajectory.dcd topology.pdb input.yaml
"""

import argparse as arg
import warnings
from datetime import datetime

import MDAnalysis as mda
import numpy as np
from loguru import logger
from MDAnalysis.exceptions import NoDataError
from MDAnalysis.guesser.tables import vdwradii
from MDAnalysis.lib.util import guess_format
from MDAnalysis.topology.LAMMPSParser import LammpsDumpParser

from .config import MolClsConfig, read_config
from .log import start_logging
from .molclusters import MolClusters
from .version import __version__


def _apply_lammps_resnames(uni: mda.Universe, config: MolClsConfig) -> None:
    """Fill in residue names on a LAMMPS topology from the `lammps_resnames` config.

    LAMMPS files carry no residue names, only a numeric molecule id per atom
    (exposed by MDAnalysis as `resid`). This assigns each residue's `resname` by
    looking up its `resid` in `config.lammps_resnames`, so the rest of the
    pipeline can treat a LAMMPS system exactly like a GROMACS one. A no-op when
    `lammps_resnames` isn't set.

    Parameters
    ----------
    uni : mda.Universe
        The Universe to update in place.
    config : MolClsConfig
        The analysis configuration.

    Raises
    ------
    ValueError
        If some residue's molecule id isn't covered by `lammps_resnames`.
    """
    if config.lammps_resnames is None:
        return

    resids = uni.residues.resids
    names = [config.resname_for_resid(int(resid)) for resid in resids]

    missing = sorted({int(r) for r, n in zip(resids, names, strict=True) if n is None})
    if missing:
        raise ValueError(
            "'lammps_resnames' does not cover LAMMPS molecule id(s) "
            f"{missing}. Add an entry for every molecule id present in the topology."
        )

    if hasattr(uni.residues, "resnames"):
        uni.residues.resnames = names
    else:
        uni.add_TopologyAttr("resnames", values=names)

    logger.debug(f"Assigned resnames from 'lammps_resnames' to {len(names)} residues.")


def _apply_lammps_dump_elements(uni: mda.Universe, traj: str) -> None:
    """Fill in elements on a topology from the `element` column of a LAMMPS dump.

    MDAnalysis' LAMMPS dump *trajectory* reader only reads coordinates, so when a
    DATA topology (which has no elements or atom names) is paired with a dump
    trajectory, the dump's `element` column (``dump_modify ... element ...``) is
    lost and elements can't be guessed. This reads that column from the dump's
    first frame and assigns it by atom id. A no-op when the topology already has
    elements, when `traj` isn't a LAMMPS dump, or when the dump has no `element`
    column.

    Parameters
    ----------
    uni : mda.Universe
        The Universe to update in place.
    traj : str
        The trajectory file the Universe was built from.

    Raises
    ------
    ValueError
        If the dump's atom ids don't match the topology's.
    """
    if hasattr(uni.atoms, "elements"):
        return
    try:
        if guess_format(traj) != "LAMMPSDUMP":
            return
    except ValueError:
        return

    with warnings.catch_warnings():
        # only the elements are used, so the parser's mass/type fallbacks don't apply
        warnings.filterwarnings(
            "ignore", message="No mass column|Guessed all Masses|Set all atom types"
        )
        top = LammpsDumpParser(traj).parse()
    if not hasattr(top, "elements"):
        return

    dump_ids = top.ids.values  # sorted by the parser
    idx = np.searchsorted(dump_ids, uni.atoms.ids)
    idx[idx == len(dump_ids)] = 0
    if len(dump_ids) != len(uni.atoms) or np.any(dump_ids[idx] != uni.atoms.ids):
        raise ValueError(
            f"The atom ids in the LAMMPS dump {traj!r} don't match the topology's."
        )

    uni.add_TopologyAttr("elements", values=top.elements.values[idx])
    logger.debug(f"Assigned elements from the LAMMPS dump {traj!r}.")


def main() -> None:
    """Main entry point for the MolClusters analysis script.

    This function sets up the command-line interface, parses arguments, initializes the MDAnalysis Universe,
    and runs the molecular cluster analysis using the MolClusters library.

    Command-line Arguments:
    -----------------------
    - traj: The trajectory file for the molecular dynamics simulation.
    - top: The topology file for the molecular dynamics simulation.
    - inp: The input YAML/JSON/TOML file containing analysis settings.
    - --traj-memory: Load the trajectory into memory.
    - --in-memory-step: Step for in-memory trajectory loading (requires --traj-memory).
    - --version: Displays the version of the MolClusters library.

    Raises
    ------
    KeyError
        If an atom in the topology does not have an associated element.
    ValueError
        If the topology has neither elements nor atom names to guess them from.
    """  # noqa: D401
    parser = arg.ArgumentParser()

    parser.add_argument("traj", type=str, help="Trajectory File.")
    parser.add_argument("top", type=str, help="Topology file.")
    parser.add_argument("inp", type=str, help="Input file with analysis settings.")
    parser.add_argument(
        "--traj-memory",
        action="store_true",
        help="Load trajectory in memory. Use with caution.",
    )
    parser.add_argument(
        "--in-memory-step",
        type=int,
        default=1,
        help="Step for in-memory trajectory loading. Default is 1.",
    )
    parser.add_argument("--version", action="version", version=__version__)

    args = parser.parse_args()

    if not args.traj_memory and args.in_memory_step != 1:
        parser.error("--in-memory-step can only be used when --traj-memory is enabled.")

    start_logging(filename=f"molclusters_{datetime.now().strftime('%Y%m%d_%H%M')}.log")

    logger.info("Reading configuration")
    cls_args = read_config(args.inp)

    logger.info("Starting Cluster Analysis")
    uni = mda.Universe(
        args.top,
        args.traj,
        in_memory=args.traj_memory,
        in_memory_step=args.in_memory_step,
    )

    _apply_lammps_resnames(uni, cls_args)

    if cls_args.solvent is None and cls_args.solute is not None:
        cls_args.solvent = sorted(set(uni.residues.resnames) - set(cls_args.solute))
        logger.debug(f"Setting solvent to {cls_args.solvent}")

    _apply_lammps_dump_elements(uni, args.traj)

    if not hasattr(uni.atoms, "elements"):
        if not hasattr(uni.atoms, "names"):
            # e.g. a LAMMPS DATA topology with a trajectory that has no element column
            raise ValueError(
                "The topology has neither elements nor atom names to guess them from. "
                "For LAMMPS, write an 'element' column to the dump trajectory "
                "(dump_modify ... element ...)."
            )
        # e.g. TPR topologies carry types and masses but no elements
        logger.debug("Topology has no elements, guessing them from atom names.")
        uni.guess_TopologyAttrs(to_guess=["elements"])

    try:
        radiis = []
        for at in uni.atoms:
            if at.element == "":
                at.element = "".join(
                    filter(str.isalpha, at.name)
                )  # Extract only the alphabetic part of the name
            radiis.append(vdwradii[at.element.upper()])
    except (KeyError, NoDataError) as err:
        raise KeyError(f"Atom: {str(at)} does not have an element.") from err

    uni.add_TopologyAttr("radii", values=radiis)

    molclusters = MolClusters(uni, cls_args)
    molclusters.run()
