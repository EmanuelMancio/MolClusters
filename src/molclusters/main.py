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
from datetime import datetime

import MDAnalysis as mda
from loguru import logger
from MDAnalysis.guesser.tables import vdwradii

from .config import read_config
from .log import start_logging
from .molclusters import MolClusters
from .version import __version__


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

    if cls_args.solvent is None and cls_args.solute is not None:
        cls_args.solvent = set(uni.residues.resnames) - set(cls_args.solute)
        logger.debug(f"Setting solvent to {cls_args.solvent}")

    try:
        radiis = []
        for at in uni.atoms:
            if at.element == "":
                at.element = "".join(
                    filter(str.isalpha, at.name)
                )  # Extract only the alphabetic part of the name
            radiis.append(vdwradii[at.element.upper()])
    except KeyError as err:
        raise KeyError(f"Atom: {str(at)} does not have an element.") from err

    uni.add_TopologyAttr("radii", values=radiis)

    molclusters = MolClusters(uni, cls_args)
    molclusters.run()
