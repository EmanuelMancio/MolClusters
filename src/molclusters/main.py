import argparse as arg
import copy
import MDAnalysis as mda
import pathlib as path
import yaml

from typing import Type, Dict, Tuple

from .molclusters import MolClusters
from . import __version__


def parse_input_file(
    in_file: Type[arg.FileType],
) -> Dict[str, Dict[str, Tuple[float, str]]]:
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

    return config


def main():
    parser = arg.ArgumentParser()

    parser.add_argument("traj", type=str, help="Trajectory File")
    parser.add_argument("top", type=str, help="Topology file")
    parser.add_argument(
        "inp", type=arg.FileType("r"), help="Input file with analysis settings"
    )
    parser.add_argument("--version",action="version",version=__version__)

    args = parser.parse_args()

    cls_args = parse_input_file(args.inp)

    uni = mda.Universe(
        args.top, args.traj, in_memory=True
    )  # TODO: add in_memory_step as option on cmdline

    molclusters = MolClusters(uni, cls_args)
    molclusters.run()
