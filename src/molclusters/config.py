# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides the `MolClsConfig` class that stores the configuration information for the analysis."""

import json
import tomllib
from collections.abc import Iterable
from enum import StrEnum, auto
from pathlib import Path
from typing import Annotated, Self

import yaml
from loguru import logger
from pydantic import Field, PositiveInt, PrivateAttr, model_validator
from pydantic.dataclasses import dataclass
from pydantic_settings import BaseSettings, SettingsConfigDict

from .log import _logger_wraps
from .symdict import SymmetricDict


class RuleType(StrEnum):
    """Enum class for types of rules."""

    CM = auto()
    HB = auto()


# A distance in Angstrom: positive and finite (rules are written as plain-text
# numbers, so `inf`/`nan` must be rejected explicitly rather than relying on a
# hand-rolled number regex).
_Distance = Annotated[float, Field(gt=0.0, allow_inf_nan=False)]


@dataclass
class Rule:
    """Base rule class."""

    type: RuleType


@dataclass(kw_only=True)
class CMRule(Rule):
    """Rule class for `cm` rule."""

    dist: _Distance
    type: RuleType = Field(RuleType.CM, frozen=True)


@dataclass(kw_only=True)
class HBRule(Rule):
    """Rule class for `hb` rule."""

    dist: _Distance = 3.5
    ang: Annotated[float, Field(ge=0.0, le=180.0)] = 150.0
    type: RuleType = Field(RuleType.HB, frozen=True)


_HB_FLAGS = {"d": "dist", "a": "ang"}


def _to_float(raw: str, label: str) -> float:
    try:
        return float(raw)
    except ValueError:
        raise ValueError(f"Invalid number for '{label}': {raw!r}.") from None


def _parse_rule(spec: str) -> Rule:
    """
    Parse a rule string such as ``"cm 5.0"`` or ``"hb d 3.5 a 150"`` into a `Rule`.

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
            f"'hb' rule only supports flags 'd' and 'a', got {sorted(unknown)}."
        )

    return HBRule(
        **{_HB_FLAGS[flag]: _to_float(value, flag) for flag, value in flags.items()}
    )


# One `lammps_resnames` entry: a single LAMMPS molecule id (``501``), or an
# inclusive id range written ``"first-last"`` (``"1-500"``).
_LammpsResidSpec = int | str


def _expand_lammps_resid_spec(spec: _LammpsResidSpec, name: str) -> list[int]:
    """
    Expand one `lammps_resnames` entry into the LAMMPS molecule ids it covers.

    Parameters
    ----------
    spec : int | str
        A single molecule id, or an inclusive range written ``"first-last"``.
    name : str
        The resname this entry is being assigned to, used in error messages.

    Returns
    -------
    list[int]
        The molecule ids covered by `spec`.

    Raises
    ------
    ValueError
        If `spec` is not a plain integer or a well-formed ``"first-last"`` range.
    """
    if isinstance(spec, int):
        return [spec]

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

        return list(range(start, end + 1))

    try:
        return [int(text)]
    except ValueError:
        raise ValueError(
            f"Invalid LAMMPS molecule id {spec!r} for {name!r}: expected a whole "
            "number or a range like '1-500'."
        ) from None


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

    _lammps_resid_to_name: dict[int, str] = PrivateAttr(default_factory=dict)

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
        return self._lammps_resid_to_name.get(resid)

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

    @model_validator(mode="after")
    def _build_rules(self) -> Self:
        for mi, neighbors in self.rules.items():
            for mj, spec in neighbors.items():
                if "solute" in (mi, mj):
                    raise ValueError(
                        "'solute' is not supported in 'rules' configuration."
                    )

                try:
                    rule = _parse_rule(spec)
                except ValueError as e:
                    e.add_note(f"From rule {mi}:{mj}")
                    raise

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

        for name, spec in self.lammps_resnames.items():
            specs = spec if isinstance(spec, list) else [spec]
            for one in specs:
                for resid in _expand_lammps_resid_spec(one, name):
                    other = self._lammps_resid_to_name.get(resid)
                    if other is not None and other != name:
                        raise ValueError(
                            f"LAMMPS molecule id {resid} is assigned to both "
                            f"{other!r} and {name!r} in 'lammps_resnames'."
                        )
                    self._lammps_resid_to_name[resid] = name

        return self


@_logger_wraps(entry=False, exit=False)
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

    config = MolClsConfig(**data)

    logger.info(f"Configuration read:\n{config}")

    return config
