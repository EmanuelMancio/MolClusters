# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Serves as the main entry point for analyzing molecular clusters using the MolClusters library.

Functions:
----------
- parse_input_file(in_file: Type[arg.FileType]) -> Dict[str, Dict[str, Tuple[float, str]]]:
    Parses the input YAML configuration file and processes the rules for molecular cluster analysis.

- main():
    The main function that sets up the command-line interface, parses arguments, initializes the MDAnalysis Universe,
    and runs the molecular cluster analysis.

Dependencies:
-------------
- argparse: For parsing command-line arguments.
- copy: For deep copying configuration data.
- pathlib: For handling file paths.
- typing: For type annotations.
- MDAnalysis: For molecular dynamics trajectory and structure analysis.
- yaml: For parsing YAML configuration files.
- MolClusters: The core library for molecular cluster analysis.

Usage:
------
Run this script from the command line with the required arguments:
    python main.py <trajectory_file> <topology_file> <input_file>

Arguments:
----------
- traj: The trajectory file for the molecular dynamics simulation.
- top: The topology file for the molecular dynamics simulation.
- inp: The input YAML file containing analysis settings.

Optional Arguments:
-------------------
- --version: Displays the version of the MolClusters library.

Example:
--------
    python main.py trajectory.dcd topology.pdb input.yaml
"""

import argparse as arg
import copy
import pathlib as path
from typing import Dict, List, Tuple

import MDAnalysis as mda
import yaml
from MDAnalysis.guesser.tables import vdwradii

from . import __version__
from .molclusters import MolClusters


# TODO: add Config class for better config capability
def parse_input_file(
    in_file: arg.FileType,
) -> Dict[str, Dict[str, Tuple[str, Dict[str, float]]] | List[str]]:
    """Parse the input YAML configuration file for molecular cluster analysis.

    This function reads the input YAML file, processes the rules for molecular cluster analysis,
    and returns a structured configuration dictionary.

    Parameter
    ----------
    in_file : arg.FileType
        The input YAML file containing analysis settings.

    Returns
    -------
    Dict[str, Dict[str, Tuple[str, Dict[str, float]]] | List[str]]:
        A dictionary containing the parsed configuration, including rules and other settings.

    Raises
    ------
    KeyError:
        If the input file contains invalid or missing keys.
    """
    config = yaml.safe_load(in_file)
    config["filename"] = str(path.Path(in_file.name).absolute())

    for mi, val in copy.deepcopy(config["rules"]).items():
        for mj, rule in val.items():
            rule = rule.split()
            op = rule[0].lower()
            if op == "cm":
                dist = float(rule[1])
                if mj in config["rules"]:
                    config["rules"][mj][mi] = (op, dist)
                else:
                    config["rules"][mj] = {mi: (op, dist)}
                config["rules"][mi][mj] = (op, dist)
            else:
                dist = float(rule[rule.index("d") + 1]) if "d" in rule else 3.5
                ang = float(rule[rule.index("a") + 1]) if "a" in rule else 150.0

                if mj in config["rules"]:
                    config["rules"][mj][mi] = (op, {"d": dist, "a": ang})
                else:
                    config["rules"][mj] = {mi: (op, {"d": dist, "a": ang})}
                config["rules"][mi][mj] = (op, {"d": dist, "a": ang})

    # TODO: Implement input correctness analysis

    if "nucleus" in config:
        if "solute" in config["nucleus"]:
            if "solute" in config:
                config["nucleus"].extend(config["solute"])
            else:
                print(
                    "'solute' in nucleus being desconsidered because solute was not defined!"
                )

            config["nucleus"].remove("solute")

    return config


def main() -> None:
    """Main entry point for the MolClusters analysis script.

    This function sets up the command-line interface, parses arguments, initializes the MDAnalysis Universe,
    and runs the molecular cluster analysis using the MolClusters library.

    Command-line Arguments:
    -----------------------
    - traj: The trajectory file for the molecular dynamics simulation.
    - top: The topology file for the molecular dynamics simulation.
    - inp: The input YAML file containing analysis settings.
    - --version: Displays the version of the MolClusters library.

    Raises
    ------
    KeyError:
        If an atom in the topology does not have an associated element.
    """  # noqa: D401
    parser = arg.ArgumentParser()

    parser.add_argument("traj", type=str, help="Trajectory File")
    parser.add_argument("top", type=str, help="Topology file")
    parser.add_argument(
        "inp", type=arg.FileType("r"), help="Input file with analysis settings"
    )
    parser.add_argument("--version", action="version", version=__version__)

    args = parser.parse_args()

    cls_args = parse_input_file(args.inp)

    uni = mda.Universe(
        args.top, args.traj, in_memory=True
    )  # TODO: add in_memory_step as option on cmdline

    try:
        radiis = []
        for at in uni.atoms:
            if at.element == "":
                at.element = at.name[
                    0
                ]  # FIXME: will stop working for atoms with names longer than 1 char
            radiis.append(vdwradii[at.element])
    except KeyError as err:
        raise KeyError(f"Atom: {str(at)} does not have a element.") from err

    uni.add_TopologyAttr("radii", values=radiis)

    molclusters = MolClusters(uni, cls_args)
    molclusters.run()
