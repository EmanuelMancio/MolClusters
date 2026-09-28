# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Provides the `MolClsConfig` class that stores the configuration information for the analysis."""

import bisect
import json
import tomllib
from collections.abc import Iterable
from enum import StrEnum, auto
from pathlib import Path
from typing import Annotated, Literal, Self

import yaml
from loguru import logger
from pydantic import Field, PositiveInt, PrivateAttr, model_validator
from pydantic.dataclasses import dataclass
from pydantic_settings import BaseSettings, SettingsConfigDict

from .symdict import SymmetricDict


class RuleType(StrEnum):
    """Enum class for types of rules."""

    CM = auto()
    HB = auto()


# A distance in Angstrom: positive and finite (rules are written as plain-text
# numbers, so `inf`/`nan` must be rejected explicitly rather than relying on a
# hand-rolled number regex).
_Distance = Annotated[float, Field(gt=0.0, allow_inf_nan=False)]

# The atomic mass (in u) and partial charge (in e) thresholds of an 'hb' rule's
# hydrogen/acceptor criteria.
_Mass = Annotated[float, Field(ge=0.0, allow_inf_nan=False)]
_Charge = Annotated[float, Field(allow_inf_nan=False)]

# The MDAnalysis acceleration backends usable by this project's pinned MDAnalysis
# version without extra optional dependencies (excludes "distopia", which needs the
# separate `distopia` package). Restricting to a Literal here validates the value
# eagerly at config-read time; left to MDAnalysis itself, an invalid backend can go
# unnoticed indefinitely, since it's only checked when a distance call actually falls
# back to its "bruteforce"/"pkdtree" methods (see `ConnectionTable`'s docstring).
DistanceBackend = Literal["serial", "OpenMP"]

# How the JSON report (molclusters.jsonl) is compressed: zstd (.zst), gzip (.gz),
# or not at all.
ReportCompression = Literal["zstd", "gzip", "none"]


@dataclass
class Rule:
    """Base rule class."""

    type: RuleType


@dataclass(kw_only=True)
class CMRule(Rule):
    """Rule class for `cm` rule."""

    dist: _Distance
    type: RuleType = Field(RuleType.CM, frozen=True)

    def __str__(self) -> str:
        """Format the rule in the config file's own syntax.

        Returns
        -------
        str
            E.g. ``cm 5.0``.
        """
        return f"cm {self.dist}"


@dataclass(kw_only=True)
class HBRule(Rule):
    """Rule class for `hb` rule.

    Besides the donor-acceptor distance and D-H-A angle cutoffs, it holds the
    criteria that pick out the rule's hydrogens (mass strictly between
    `h_mass_min` and `h_mass_max`, charge above `h_charge_min`) and acceptors
    (charge below `a_charge_max`). The defaults are MDAnalysis'
    `HydrogenBondAnalysis.guess_hydrogens()`/`guess_acceptors()` ones.
    """

    dist: _Distance = 3.5
    ang: Annotated[float, Field(ge=0.0, le=180.0)] = 150.0
    h_mass_min: _Mass = 0.9
    h_mass_max: _Mass = 1.1
    h_charge_min: _Charge = 0.3
    a_charge_max: _Charge = -0.5
    type: RuleType = Field(RuleType.HB, frozen=True)

    @model_validator(mode="after")
    def _check_mass_bounds(self) -> Self:
        if self.h_mass_min >= self.h_mass_max:
            raise ValueError(
                f"Hydrogen mass bounds must satisfy hmin < hmax, got "
                f"hmin {self.h_mass_min} and hmax {self.h_mass_max}."
            )
        return self

    def __str__(self) -> str:
        """Format the rule in the config file's own syntax, defaults filled in.

        Returns
        -------
        str
            E.g. ``hb d 3.5 a 150.0 hmin 0.9 hmax 1.1 hq 0.3 aq -0.5``.
        """
        return " ".join(
            [
                "hb",
                *(f"{flag} {getattr(self, attr)}" for flag, attr in _HB_FLAGS.items()),
            ]
        )


# each 'hb' rule flag, in the order `HBRule.__str__` writes them, and the field it sets
_HB_FLAGS = {
    "d": "dist",
    "a": "ang",
    "hmin": "h_mass_min",
    "hmax": "h_mass_max",
    "hq": "h_charge_min",
    "aq": "a_charge_max",
}


def _to_float(raw: str, label: str) -> float:
    try:
        return float(raw)
    except ValueError:
        raise ValueError(f"Invalid number for '{label}': {raw!r}.") from None


def _parse_rule(spec: str) -> Rule:
    """
    Parse a rule string such as ``"cm 5.0"`` or ``"hb d 3.5 a 150"`` into a `Rule`.

    An 'hb' rule takes any of these ``<flag> <number>`` pairs, in any order, each
    optional: ``d`` (donor-acceptor distance cutoff), ``a`` (D-H-A angle cutoff),
    ``hmin``/``hmax`` (exclusive hydrogen mass bounds), ``hq`` (exclusive minimum
    hydrogen charge) and ``aq`` (exclusive maximum acceptor charge).

    Parameters
    ----------
    spec : str
        The raw rule string as written in the configuration file.

    Returns
    -------
    Rule
        The parsed `CMRule` or `HBRule`.

    Raises
    ------
    ValueError
        If `spec` is empty, uses an unknown command, or has an invalid/malformed
        argument for its command.
    """
    parts = spec.split()

    if not parts:
        raise ValueError("Empty rule is not accepted.")

    try:
        rule_type = RuleType(parts[0])
    except ValueError:
        raise ValueError("Command must be either 'cm' or 'hb'.") from None

    args = parts[1:]

    if rule_type == RuleType.CM:
        if len(args) != 1:
            raise ValueError("Format for 'cm' rule must be: cm <number>")
        return CMRule(dist=_to_float(args[0], "cm"))

    if len(args) % 2:
        raise ValueError("'hb' flags must be given as '<flag> <number>' pairs.")

    flags = dict(zip(args[0::2], args[1::2], strict=True))
    if len(flags) != len(args) // 2:
        raise ValueError("Each 'hb' flag can only be used once.")

    unknown = flags.keys() - _HB_FLAGS.keys()
    if unknown:
        raise ValueError(
            f"'hb' rule only supports flags {list(_HB_FLAGS)}, got {sorted(unknown)}."
        )

    return HBRule(
        **{_HB_FLAGS[flag]: _to_float(value, flag) for flag, value in flags.items()}
    )


# One `lammps_resnames` entry: a single LAMMPS molecule id (``501``), or an
# inclusive id range written ``"first-last"`` (``"1-500"``).
_LammpsResidSpec = int | str


def _parse_lammps_resid_range(spec: _LammpsResidSpec, name: str) -> tuple[int, int]:
    """
    Parse one `lammps_resnames` entry into an inclusive (start, end) id range.

    Kept as a `(start, end)` pair rather than the individual ids it covers, since a
    range like ``"1-500000"`` is meant to describe a whole block of solvent
    molecules cheaply, not to be expanded into half a million entries.

    Parameters
    ----------
    spec : int | str
        A single molecule id, or an inclusive range written ``"first-last"``.
    name : str
        The resname this entry is being assigned to, used in error messages.

    Returns
    -------
    tuple[int, int]
        The inclusive `(start, end)` range covered by `spec` (`start == end` for a
        single id).

    Raises
    ------
    ValueError
        If `spec` is not a plain integer or a well-formed ``"first-last"`` range.
    """
    if isinstance(spec, int):
        return spec, spec

    text = spec.strip()

    if "-" in text:
        start_s, _, end_s = text.partition("-")
        try:
            start, end = int(start_s), int(end_s)
        except ValueError:
            raise ValueError(
                f"Invalid LAMMPS molecule id range {spec!r} for {name!r}: expected "
                "'first-last', e.g. '1-500'."
            ) from None

        if start > end:
            raise ValueError(
                f"Invalid LAMMPS molecule id range {spec!r} for {name!r}: the first "
                "id must not be greater than the last."
            )

        return start, end

    try:
        value = int(text)
    except ValueError:
        raise ValueError(
            f"Invalid LAMMPS molecule id {spec!r} for {name!r}: expected a whole "
            "number or a range like '1-500'."
        ) from None

    return value, value


# `lammps_timestep` units, and how many ps each is
_TIME_UNITS_PS = {"fs": 1e-3, "ps": 1.0, "ns": 1e3}


def _parse_lammps_timestep(spec: str | float) -> float:
    """
    Parse `lammps_timestep`, a LAMMPS timestep size with its unit, into ps.

    The unit is required: LAMMPS' own depends on the run's unit style (fs for
    ``real``, ps for ``metal``), so a bare number would be ambiguous.

    Parameters
    ----------
    spec : str | float
        The timestep size and its unit, e.g. ``"2 fs"`` or ``"0.001 ps"``.

    Returns
    -------
    float
        The timestep size, in ps.

    Raises
    ------
    ValueError
        If `spec` isn't a positive number followed by fs, ps or ns.
    """
    units = ", ".join(_TIME_UNITS_PS)
    if not isinstance(spec, str):
        raise ValueError(
            f"lammps_timestep {spec!r} needs its unit ({units}), e.g. '2 fs' for "
            "LAMMPS' real units or '0.001 ps' for metal units."
        )

    text = spec.strip()
    unit = text[-2:].lower()
    try:
        size = float(text[:-2])
    except ValueError:
        size = None
    if unit not in _TIME_UNITS_PS or size is None:
        raise ValueError(
            f"Invalid lammps_timestep {spec!r}: expected a number and its unit "
            f"({units}), e.g. '2 fs'."
        )
    if size <= 0:
        raise ValueError(f"Invalid lammps_timestep {spec!r}: it must be positive.")

    return size * _TIME_UNITS_PS[unit]


class MolClsConfig(BaseSettings):
    """Settings class for input parameters."""

    model_config = SettingsConfigDict(
        env_prefix="MOLCLS_", use_enum_values=True, extra="ignore"
    )

    rules: dict[str, dict[str, str]]
    _rules: SymmetricDict[str, Rule] = PrivateAttr(default_factory=SymmetricDict)

    solute: list[str] | None = None
    solvent: list[str] | None = None

    nucleus: list[str] | None = None

    follow: list[str | PositiveInt] | None = None
    _follow_solute: bool = PrivateAttr(default=False)

    ignore_composition: list[list[str]] | None = None

    _ignore_composition: set[frozenset[str]] = PrivateAttr(default_factory=set)

    # LAMMPS topologies have no residue names, only numeric molecule ids: this maps
    # each name used elsewhere in this file (rules, solute, ...) to the LAMMPS
    # molecule id(s) it stands for, e.g. {"SOL": "1-500", "NA": 501}.
    lammps_resnames: dict[str, _LammpsResidSpec | list[_LammpsResidSpec]] | None = None

    # The LAMMPS run's timestep size with its unit, e.g. "2 fs". LAMMPS dumps only
    # record step numbers, so without it their times are step numbers, not ps.
    # Optional; used when MolClusters loads a LAMMPS dump (the CLI does).
    lammps_timestep: str | float | None = None
    _lammps_timestep_ps: float | None = PrivateAttr(default=None)

    # Acceleration backend for the "cm" rule's distance calculations (see
    # `ConnectionTable`); "hb" rules always run serially, since HydrogenBondAnalysis
    # doesn't expose a backend option.
    distance_backend: DistanceBackend = "serial"

    # Turns built-in analyses on or off by class name, e.g. {"JsonReport": False}.
    # One left out runs when the options it needs are set (see `analysis.builtins`).
    analyses: dict[str, bool] = Field(default_factory=dict)

    # How `JsonReport` compresses its report: zstd is the fastest and smallest,
    # gzip the most widely readable (see `report` for the file's format).
    report_compression: ReportCompression = "zstd"

    # Sorted, non-overlapping (start, end, name) ranges built from `lammps_resnames`,
    # kept as ranges (not one dict entry per id) so a config spanning millions of
    # molecule ids costs only as much memory as the handful of lines the user wrote.
    _lammps_resid_ranges: list[tuple[int, int, str]] = PrivateAttr(default_factory=list)
    _lammps_resid_starts: list[int] = PrivateAttr(default_factory=list)

    def resname_for_resid(self, resid: int) -> str | None:
        """
        Look up the resname configured for a LAMMPS molecule id.

        Parameters
        ----------
        resid : int
            A LAMMPS molecule id (the ``mol`` column), which MDAnalysis exposes as
            `resid` for a topology parsed from a LAMMPS DATA file.

        Returns
        -------
        str | None
            The name assigned to `resid` via `lammps_resnames`, or None if it isn't
            covered.
        """
        idx = bisect.bisect_right(self._lammps_resid_starts, resid) - 1
        if idx < 0:
            return None

        start, end, name = self._lammps_resid_ranges[idx]
        return name if start <= resid <= end else None

    @property
    def lammps_timestep_ps(self) -> float | None:
        """
        The LAMMPS timestep size, in ps.

        Returns
        -------
        float | None
            `lammps_timestep` converted to ps, or None if it isn't set.
        """
        return self._lammps_timestep_ps

    def is_ignored_composition(self, resnames: Iterable[str]) -> bool:
        """
        Check whether a cluster composition is configured to be ignored.

        Parameters
        ----------
        resnames : Iterable[str]
            The residue names making up a candidate cluster.

        Returns
        -------
        bool
            True if this composition matches an `ignore_composition` entry.
        """
        return frozenset(resnames) in self._ignore_composition

    def describe(self) -> str:
        """
        Describe the configuration the analysis actually runs with.

        Unlike the raw config file, this shows the parsed rules with their
        defaults filled in, the `solute` keyword expanded, and the LAMMPS molecule
        id ranges merged, so a user can check what was understood.

        Returns
        -------
        str
            A multi-line, human-readable description.
        """

        def names(values: Iterable[str] | None) -> str:
            return ", ".join(values) if values else "none"

        lines = ["rules (distances in angstrom, angles in degrees):"]
        lines += [
            f"  {mi} - {mj}: {self._rules[mi, mj]}" for mi, mj in sorted(self._rules)
        ]
        lines.append(f"solute: {names(self.solute)}")
        lines.append(f"solvent: {names(self.solvent)}")
        lines.append(f"nucleus: {names(self.nucleus)}")
        lines.append(
            "follow: "
            + (
                "solute (one solute-<resid>.gro per solute)"
                if self._follow_solute
                else "none"
            )
        )

        comps = self.ignore_composition or []
        lines.append(
            f"ignore_composition: {'; '.join(' + '.join(sorted(c)) for c in comps) or 'none'}"
        )
        lines.append(f"distance_backend: {self.distance_backend} (cm rules only)")

        # both import this module
        from .analysis.builtins import BUILTINS
        from .report import report_name

        running = [b.name for b in BUILTINS if b.runs(self)]
        off = [b.name for b in BUILTINS if not self.analyses.get(b.name, True)]
        lines.append(
            f"analyses: {', '.join(running) or 'none'}"
            + (f" (turned off: {', '.join(off)})" if off else "")
        )
        lines.append(
            f"report_compression: {self.report_compression} "
            f"({report_name(self.report_compression)})"
        )

        ranges = [
            f"{name} = {start}" if start == end else f"{name} = {start}-{end}"
            for start, end, name in self._lammps_resid_ranges
        ]
        lines.append(f"lammps_resnames: {', '.join(ranges) or 'none'}")
        lines.append(
            "lammps_timestep: none (LAMMPS dump times are step numbers)"
            if self._lammps_timestep_ps is None
            else f"lammps_timestep: {self._lammps_timestep_ps:g} ps"
        )

        return "\n".join(lines)

    @model_validator(mode="after")
    def _build_lammps_timestep(self) -> Self:
        if self.lammps_timestep is not None:
            self._lammps_timestep_ps = _parse_lammps_timestep(self.lammps_timestep)

        return self

    @model_validator(mode="after")
    def _build_rules(self) -> Self:
        # each pair's rule and how it was written ("A:B ('cm 5.0')"): a rule
        # applies both ways, so A:B and B:A are the same pair, and a conflict
        # names both
        given: SymmetricDict[str, tuple[Rule, str]] = SymmetricDict()
        for mi, neighbors in self.rules.items():
            for mj, spec in neighbors.items():
                if "solute" in (mi, mj):
                    raise ValueError(
                        "'solute' is not supported in 'rules' configuration."
                    )

                try:
                    rule = _parse_rule(spec)
                except ValueError as e:
                    # pydantic builds its message from str(e) alone and drops
                    # notes, so the rule has to be named in the message itself
                    raise ValueError(f"Invalid rule {mi}:{mj} ({spec!r}): {e}") from e

                if (mi, mj) in given and given[mi, mj][0] != rule:
                    raise ValueError(
                        f"Conflicting rules {given[mi, mj][1]} and {mi}:{mj} "
                        f"({spec!r}): a rule applies both ways, so give each "
                        "pair one rule."
                    )

                given[mi, mj] = (rule, f"{mi}:{mj} ({spec!r})")
                logger.trace(f"Using rule between {mi} and {mj}: {rule}")
                self._rules[mi, mj] = rule

        return self

    @model_validator(mode="after")
    def _parse_solute(self) -> Self:
        if self.solute is not None and "solute" in self.solute:
            raise ValueError("'solute' configuration cannot be solute.")

        errors = []

        if self.follow is not None:
            self._follow_solute = "solute" in self.follow

        def expand(values: list, attr_nm: str) -> list:
            if "solute" not in values:
                return values

            if self.solute is None:
                errors.append(
                    ValueError(
                        f"'solute' keyword used in {attr_nm}, but 'solute' not defined."
                    )
                )
                return values

            return list({*values, *self.solute} - {"solute"})

        for attr_nm in ("nucleus", "follow"):
            values = getattr(self, attr_nm)
            if values is not None:
                setattr(self, attr_nm, expand(values, attr_nm))

        if self.ignore_composition is not None:
            self.ignore_composition = [
                expand(values, "ignore_composition")
                for values in self.ignore_composition
            ]

        if errors:
            raise ExceptionGroup("Errors in input file", errors)

        if self.follow is not None:
            # only the `solute` keyword drives any output so far; after expansion
            # the solute resnames are indistinguishable from ones listed by hand
            unused = [
                f
                for f in self.follow
                if not (self._follow_solute and f in (self.solute or []))
            ]
            if unused:
                logger.warning(
                    f"'follow' entries {unused} have no effect: only the 'solute' "
                    "keyword is supported so far."
                )

        return self

    @model_validator(mode="after")
    def _default_solvent(self) -> Self:
        # only residues named in the rules can join a cluster, so any other
        # residue is no solvent of theirs
        if self.solvent is None and self.solute is not None:
            self.solvent = sorted(self._rules.all_keys() - set(self.solute))

        return self

    @model_validator(mode="after")
    def _check_analyses(self) -> Self:
        from .analysis.builtins import BUILTINS  # the analysis package imports this one

        known = [b.name for b in BUILTINS]
        unknown = [name for name in self.analyses if name not in known]
        if unknown:
            # a misspelled name would otherwise silently leave an analysis on
            raise ValueError(
                f"Unknown analyses {unknown} in 'analyses': the built-ins are {known}."
            )

        missing = [
            f"{b.name} needs '{option}'"
            for b in BUILTINS
            if self.analyses.get(b.name)
            for option in b.needs
            if getattr(self, option) is None
        ]
        if missing:
            raise ValueError(f"Analyses turned on without what they need: {missing}.")

        return self

    @model_validator(mode="after")
    def _build_ignore_composition(self) -> Self:
        if self.ignore_composition is not None:
            self._ignore_composition = {frozenset(i) for i in self.ignore_composition}

        return self

    @model_validator(mode="after")
    def _build_lammps_resnames(self) -> Self:
        if self.lammps_resnames is None:
            return self

        ranges: list[tuple[int, int, str]] = []
        for name, spec in self.lammps_resnames.items():
            specs = spec if isinstance(spec, list) else [spec]
            ranges.extend(
                (*_parse_lammps_resid_range(one, name), name) for one in specs
            )

        ranges.sort(key=lambda r: r[0])

        # Ranges are user-authored (one line per species, not per molecule), so
        # this stays a handful of entries even for a huge system: an O(R^2) check
        # is simpler than a sweep and is negligible at that size.
        for i, (start_i, end_i, name_i) in enumerate(ranges):
            for start_j, end_j, name_j in ranges[i + 1 :]:
                if start_j > end_i:
                    break  # sorted by start: nothing further can overlap `i`
                if name_j != name_i:
                    lo, hi = max(start_i, start_j), min(end_i, end_j)
                    raise ValueError(
                        f"LAMMPS molecule id(s) {lo}-{hi} are assigned to both "
                        f"{name_i!r} and {name_j!r} in 'lammps_resnames'."
                    )

        merged: list[tuple[int, int, str]] = []
        for start, end, name in ranges:
            if merged and merged[-1][2] == name and start <= merged[-1][1] + 1:
                prev_start, prev_end, _ = merged[-1]
                new_end = max(prev_end, end)

                if start <= prev_end:
                    logger.warning(
                        f"'lammps_resnames' for {name!r} has overlapping entries: "
                        f"molecule id(s) {start}-{min(end, prev_end)} are listed "
                        f"more than once. They will be treated as a single range, "
                        f"{prev_start}-{new_end}, covering every id in either entry."
                    )

                merged[-1] = (prev_start, new_end, name)
            else:
                merged.append((start, end, name))

        self._lammps_resid_ranges = merged
        self._lammps_resid_starts = [r[0] for r in merged]

        return self


def read_config(path: Path | str) -> MolClsConfig:
    """
    Read a configuration file and return an MolClsConfig object.

    Support JSON, YAML, and TOML file formats based on the file extension.

    Parameters
    ----------
    path : Path
        The path to the configuration file.

    Returns
    -------
    MolClsConfig
        An instance of MolClsConfig initialized with the data from the configuration file.

    Raises
    ------
    ValueError
        If the file extension is not supported.

    Notes
    -----
    - Supported file extensions: '.json', '.yaml', '.yml', '.toml'.
    - Requires the corresponding libraries for each file type (e.g., `json`, `yaml`, `tomllib`).
    """
    path = Path(path)

    logger.info(f"Reading config file: {path.resolve()}")

    match path.suffix:
        case ".json":
            with path.open() as f:
                data = json.load(f)
        case ".yaml" | ".yml":
            with path.open() as f:
                data = yaml.safe_load(f)
        case ".toml":
            with path.open("rb") as f:
                data = tomllib.load(f)
        case _:
            raise ValueError(
                f"Unsupported config type: {path}. Only JSON, YAML, and TOML supported."
            )

    logger.debug(f"Config file contents: {data}")

    # `extra="ignore"` would otherwise silently drop a misspelled key
    unknown = sorted(data.keys() - MolClsConfig.model_fields.keys())
    if unknown:
        logger.warning(
            f"Ignoring unknown config key(s) {unknown}; check them for typos. "
            f"Known keys: {sorted(MolClsConfig.model_fields)}."
        )

    return MolClsConfig(**data)
