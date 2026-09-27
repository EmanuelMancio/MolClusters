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
- --log-level: Minimum level of the messages logged (default INFO).
- --output-dir: Directory for the results and the log (default: current directory).
- --version: Displays the version of the MolClusters library.

Example:
--------
    molclusters trajectory.dcd topology.pdb input.yaml
"""

import argparse as arg
import shlex
import sys
import warnings
from collections import Counter
from datetime import datetime
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from loguru import logger
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.coordinates.timestep import Timestep
from MDAnalysis.guesser.tables import vdwradii
from MDAnalysis.lib.util import guess_format
from MDAnalysis.topology.LAMMPSParser import LammpsDumpParser
from pydantic import ValidationError

from .config import MolClsConfig, read_config
from .log import FILE_ONLY, start_logging
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

    logger.info(
        f"Assigned residue names from 'lammps_resnames' to {len(names)} residues"
    )


# File extensions (compressed or not) read as LAMMPS dump trajectories. MDAnalysis
# only recognizes '.lammpsdump' on its own.
_LAMMPS_DUMP_EXTS = {"LAMMPSDUMP", "LAMMPSTRJ", "DUMP"}


def _traj_format(traj: str) -> str | None:
    """Pick the MDAnalysis format to read a trajectory file with.

    Parameters
    ----------
    traj : str
        The trajectory file.

    Returns
    -------
    str | None
        ``"LAMMPSDUMP"`` for a LAMMPS dump (see `_LAMMPS_DUMP_EXTS`), or None to
        let MDAnalysis guess the format from the extension.
    """
    try:
        ext = guess_format(traj)
    except ValueError:  # no extension
        return None
    return "LAMMPSDUMP" if ext in _LAMMPS_DUMP_EXTS else None


def _lammps_dump_timestep(config: MolClsConfig, traj: str) -> dict[str, float]:
    """Give the trajectory reader the LAMMPS timestep size, for a LAMMPS dump.

    A LAMMPS dump records step numbers only, and MDAnalysis turns them into times
    as step x dt with dt = 1 unless told otherwise, so times are step numbers
    unless `lammps_timestep` is set. Warns when that's the case, and when
    `lammps_timestep` is set but `traj` isn't a LAMMPS dump.

    Parameters
    ----------
    config : MolClsConfig
        The configuration, with the optional `lammps_timestep`.
    traj : str
        The trajectory file.

    Returns
    -------
    dict[str, float]
        ``{"dt": timestep in ps}`` to pass to the reader, or an empty dict.
    """
    dt = config.lammps_timestep_ps
    if _traj_format(traj) != "LAMMPSDUMP":
        if dt is not None:
            logger.warning(
                f"'lammps_timestep' has no effect: {traj!r} is not a LAMMPS dump."
            )
        return {}

    if dt is None:
        logger.warning(
            "LAMMPS dumps record step numbers, not times: every time in the results "
            "(Time, cluster ages) will be a step number, not ps. Set 'lammps_timestep' "
            "in the config (e.g. '2 fs') to get them in ps."
        )
        return {}

    return {"dt": dt}


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
    if hasattr(uni.atoms, "elements") or _traj_format(traj) != "LAMMPSDUMP":
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
    logger.info(f"Assigned elements from the LAMMPS dump {traj!r}")


def _log_system(uni: mda.Universe) -> None:
    """Log a summary of the loaded system, for the user to check against expectations.

    Parameters
    ----------
    uni : mda.Universe
        The loaded Universe.
    """
    if hasattr(uni.residues, "resnames"):
        counts = Counter(uni.residues.resnames)
        composition = ", ".join(f"{n} {name}" for name, n in counts.items())
    else:
        composition = "no residue names"

    logger.info(
        f"System: {len(uni.atoms)} atoms in {len(uni.residues)} residues "
        f"({composition})"
    )

    # times as the frames carry them: the dt is inferred by MDAnalysis for some
    # formats (it warns when it can't), so they're worth a check. Not the reader's
    # dt and totaltime, which for a LAMMPS dump count MD steps, not frames
    traj = uni.trajectory
    first = traj[0].time
    every = traj[1].time - first if len(traj) > 1 else 0.0
    last = traj[-1].time
    traj[0]
    logger.info(
        f"Trajectory: {len(traj)} frame(s) from {first:g} to {last:g} ps, "
        f"every {every:g} ps"
    )

    box = uni.dimensions
    if box is None or not np.any(box[:3]):
        logger.warning(
            "The trajectory has no box: distances are computed without periodic "
            "boundary conditions."
        )
    else:
        lengths = " x ".join(f"{x:.2f}" for x in box[:3])
        angles = (
            ""
            if np.allclose(box[3:], 90)
            else " (angles " + ", ".join(f"{x:.1f}" for x in box[3:]) + ")"
        )
        logger.info(f"Box (first frame): {lengths} angstrom{angles}")


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
    - --log-level: Minimum level of the messages logged (default INFO).
    - --output-dir: Directory for the results and the log (default: current directory).
    - --version: Displays the version of the MolClusters library.

    Errors during the analysis are logged (with their traceback in the log file
    only) and end the program with exit code 1, or 130 when interrupted.
    """  # noqa: D401
    parser = arg.ArgumentParser(
        prog="molclusters",
        description=(
            f"MolClusters {__version__}: analyze the formation and lifetime of "
            "molecular clusters in a molecular dynamics trajectory."
        ),
    )

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
    parser.add_argument(
        "--log-level",
        type=str.upper,
        choices=["TRACE", "DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Minimum level of the messages logged. Default is INFO.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(),
        help=(
            "Directory to write the results and the log to, created if missing. "
            "Default is the current directory."
        ),
    )
    parser.add_argument("--version", action="version", version=__version__)

    args = parser.parse_args()

    if not args.traj_memory and args.in_memory_step != 1:
        parser.error("--in-memory-step can only be used when --traj-memory is enabled.")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    log_file = (
        args.output_dir / f"molclusters_{datetime.now().strftime('%Y%m%d_%H%M')}.log"
    )
    start_logging(level=args.log_level, filename=log_file)
    logger.info(
        f"Running in {Path.cwd()}: {shlex.join(['molclusters', *sys.argv[1:]])}"
    )

    try:
        _analyse(args)
    except KeyboardInterrupt:
        logger.warning("Interrupted: the analysis stopped before writing its results")
        sys.exit(130)
    except Exception as err:
        # the traceback is for bug reports; the terminal only needs what went wrong
        logger.bind(**FILE_ONLY).opt(exception=err).error("Traceback:")
        logger.error(
            f"{_describe_error(err)}\n(Full traceback in {log_file.resolve()})"
        )
        sys.exit(1)


def _describe_error(err: BaseException) -> str:
    """Describe an error in one message, listing each error of an exception group.

    Returns
    -------
    str
        ``<type>: <message>``, one line per error for a group or a config
        validation error, followed by the error's notes (e.g. which analysis
        raised it).
    """
    if isinstance(err, BaseExceptionGroup):
        lines = [f"{err.message}:"]
        for sub in err.exceptions:
            lines += ["  " + line for line in _describe_error(sub).splitlines()]
        return "\n".join(lines)

    if isinstance(err, ValidationError):
        # pydantic's own text adds error types, the raw input and a docs link
        lines = ["Invalid configuration:"]
        for error in err.errors():
            message = error["msg"].removeprefix("Value error, ")
            loc = ".".join(str(part) for part in error["loc"])
            lines.append(f"  {loc}: {message}" if loc else f"  {message}")
        return "\n".join(lines)

    # str(KeyError) quotes its message
    message = err.args[0] if isinstance(err, KeyError) and err.args else err
    # notes say where the error came from, e.g. which analysis raised it
    notes = "".join(f"\n{note}" for note in getattr(err, "__notes__", []))
    return f"{type(err).__name__}: {message}{notes}"


def _assign_radii(uni: mda.Universe) -> None:
    """Assign van der Waals radii from elements, filling in empty ones from names.

    Parameters
    ----------
    uni : mda.Universe
        The Universe to update in place; it must have elements and atom names.

    Raises
    ------
    KeyError
        If some element has no van der Waals radius in MDAnalysis' table, or an
        atom has no element at all.
    """
    radii = []
    n_from_name = 0
    unknown: dict[str, str] = {}  # element -> first atom that has it

    for at in uni.atoms:
        if at.element == "":
            # Extract only the alphabetic part of the name
            at.element = "".join(filter(str.isalpha, at.name))
            n_from_name += 1
        radius = vdwradii.get(at.element.upper())
        if radius is None:
            unknown.setdefault(at.element, str(at))
        radii.append(radius)

    if unknown:
        details = "; ".join(
            f"{element!r}, e.g. {atom}" if element else f"no element, e.g. {atom}"
            for element, atom in unknown.items()
        )
        raise KeyError(
            "No van der Waals radius is known for element(s) "
            f"{sorted(unknown)} ({details}). Check the elements in the topology "
            "(or, when it has none, the atom names they are guessed from)."
        )

    if n_from_name:
        logger.info(
            f"{n_from_name} atom(s) had an empty element, took it from the letters "
            "of their atom names"
        )

    uni.add_TopologyAttr("radii", values=radii)


def _load_into_memory(uni: mda.Universe, step: int | None) -> None:
    """Load the trajectory into memory, keeping every frame's time.

    Not MDAnalysis' own ``Universe(in_memory=...)`` or ``transfer_to_memory``:
    the first also hands the reader's options to the in-memory reader, which
    then gets a LAMMPS dump's dt twice, and both time frame i as i x dt from 0,
    which loses a trajectory's start time (e.g. a continuation run), uneven
    spacing, and for a LAMMPS dump even the spacing itself (its reader's dt is
    the MD timestep, not the time between frames). So the frames are copied here,
    their times with them, into a `_TimedMemoryReader`.

    Parameters
    ----------
    uni : mda.Universe
        The Universe, changed in place.
    step : int | None
        Keep every `step`-th frame (every frame when None).
    """
    traj = uni.trajectory
    if isinstance(traj, MemoryReader):  # nothing to copy
        uni.transfer_to_memory(step=step)
        return

    frames = traj[:: step or 1]
    n_frames = len(frames)
    first = traj[0]
    has_box = first.dimensions is not None
    coordinates = np.empty((n_frames, traj.n_atoms, 3), dtype=np.float32)
    dimensions = np.empty((n_frames, 6), dtype=np.float32) if has_box else None
    velocities = np.empty_like(coordinates) if first.has_velocities else None
    forces = np.empty_like(coordinates) if first.has_forces else None
    times = np.empty(n_frames)
    for i, ts in enumerate(frames):
        coordinates[i] = ts.positions
        times[i] = ts.time
        if dimensions is not None:
            dimensions[i] = ts.dimensions
        if velocities is not None:
            velocities[i] = ts.velocities
        if forces is not None:
            forces[i] = ts.forces

    uni.trajectory = _TimedMemoryReader(
        coordinates,
        dimensions=dimensions,
        # for anything reading the reader's dt: the time between frames
        dt=times[1] - times[0] if n_frames > 1 else traj.dt,
        filename=traj.filename,
        velocities=velocities,
        forces=forces,
        times=times,
    )


class _TimedMemoryReader(MemoryReader):
    """An in-memory trajectory whose frames keep their own times.

    MDAnalysis' `MemoryReader` times frame i as i x dt; this one gives each frame
    the time it had in the trajectory it was copied from (see
    `_load_into_memory`).

    Attributes
    ----------
    times : np.ndarray | None
        Each frame's time, in ps (None: i x dt, as MemoryReader).
    """

    def __init__(
        self, *args: object, times: np.ndarray | None = None, **kwargs: object
    ) -> None:
        # before MemoryReader.__init__, which reads the first frame
        self.times = None if times is None else np.asarray(times, dtype=float)
        super().__init__(*args, **kwargs)

    def _read_next_timestep(self, ts: Timestep | None = None) -> Timestep:
        ts = super()._read_next_timestep(ts)
        if self.times is not None:
            ts.time = float(self.times[ts.frame])
        return ts

    def copy(self) -> "_TimedMemoryReader":
        new = super().copy()
        new.times = None if self.times is None else self.times.copy()
        new[self.ts.frame]  # re-read, with its time
        return new


def _analyse(args: arg.Namespace) -> None:
    """Load the system described by the command-line arguments and analyse it.

    Parameters
    ----------
    args : arg.Namespace
        The parsed command-line arguments.

    Raises
    ------
    ValueError
        If the topology has neither elements nor atom names to guess them from.
    """
    cls_args = read_config(args.inp)

    logger.info(f"Loading topology {args.top!r} and trajectory {args.traj!r}")
    uni = mda.Universe(
        args.top,
        args.traj,
        format=_traj_format(args.traj),
        **_lammps_dump_timestep(cls_args, args.traj),
    )
    if args.traj_memory:
        _load_into_memory(uni, args.in_memory_step)

    _apply_lammps_resnames(uni, cls_args)
    _log_system(uni)

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
        logger.info("Topology has no elements, guessing them from atom names")
        uni.guess_TopologyAttrs(to_guess=["elements"])

    if not hasattr(uni.atoms, "names"):
        # e.g. LAMMPS topologies, which carry no atom names of their own. Without
        # this, MDAnalysis' GRO writer would warn and write every atom as "X".
        logger.info("Topology has no atom names, using elements as names")
        uni.add_TopologyAttr("names", values=uni.atoms.elements)

    _assign_radii(uni)

    molclusters = MolClusters(uni, cls_args)
    molclusters.run(output_dir=args.output_dir)
