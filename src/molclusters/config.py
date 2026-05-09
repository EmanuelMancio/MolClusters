# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides the `MolClsConfig` class that stores the configuration information for the analysis."""

import json
import re
import tomllib
from enum import StrEnum, auto
from pathlib import Path
from typing import Annotated, Literal, Self

import yaml
from loguru import logger
from pydantic import (
    BaseModel,
    Field,
    PositiveFloat,
    PositiveInt,
    PrivateAttr,
    field_validator,
    model_validator,
)
from pydantic.dataclasses import dataclass
from pydantic_settings import BaseSettings, SettingsConfigDict

from .log import _logger_wraps
from .symdict import SymmetricDict


class RuleType(StrEnum):
    """Enum class for types of rules."""

    CM = auto()
    HB = auto()


NUMBER_RE = re.compile(r"(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")


class _RuleInput(BaseModel):
    rules: dict[
        str,
        dict[str, str],
    ]

    @field_validator("rules")
    @classmethod
    def validate_rules(
        cls,
        value: dict[
            str,
            dict[str, str],
        ],
    ) -> dict[
        str,
        dict[str, str],
    ]:
        for k1 in value:
            for k2, rl in value[k1].items():
                try:
                    cls._validate_rule(rl)
                except ValueError as e:
                    e.add_note(f"From rule {k1}:{k2}")
                    raise
        return value

    @staticmethod
    def _validate_rule(value: str) -> None:
        parts = value.split()

        if not parts:
            raise ValueError("Empty rule is not accepted")

        cmd = parts[0]

        if cmd == "cm":
            if len(parts) != 2:
                raise ValueError("Formar for 'cm' rule must be: cm <number>")
            if not NUMBER_RE.fullmatch(parts[1]):
                raise ValueError(
                    "cm requires a valid positive number (float or scientific notation allowed)."
                )

        elif cmd == "hb":
            used_flags = set()
            i = 1

            while i < len(parts):
                flag = parts[i]

                if flag not in {"d", "a"}:
                    raise ValueError("'hb' rule only supporrs flags 'd' and 'a'.")

                if flag in used_flags:
                    raise ValueError(f"Flag '{flag}' cannot be used twice.")

                if i + 1 >= len(parts):
                    raise ValueError(f"Flag '{flag}' must be followed by a number.")

                number = parts[i + 1]

                if not NUMBER_RE.fullmatch(number):
                    raise ValueError(
                        f"Invalid number for flag '{flag}'. Must be positive (float or scientific notation allowed)."
                    )

                used_flags.add(flag)
                i += 2
        else:
            raise ValueError("Command must be either 'cm' or 'hb'.")


@dataclass
class Rule:
    """Base rule class."""

    type: RuleType


@dataclass(kw_only=True)
class CMRule(Rule):
    """Rule class for `cm` rule."""

    dist: PositiveFloat
    type: RuleType = Field(RuleType.CM, frozen=True)


@dataclass(kw_only=True)
class HBRule(Rule):
    """Rule class for `hb` rule."""

    dist: PositiveFloat
    ang: Annotated[float, Field(ge=0.0, le=180.0)]
    type: RuleType = Field(RuleType.HB, frozen=True)


class MolClsConfig(BaseSettings, _RuleInput):
    """Settings class for input parameters."""

    model_config = SettingsConfigDict(
        env_prefix="MOLCLS_", use_enum_values=True, extra="ignore"
    )

    _rules: SymmetricDict[str, Rule] = PrivateAttr(default=SymmetricDict())

    solute: list[str] | None = None
    solvent: list[str] | None = None

    nucleus: list[str] | None = None

    follow: list[str | PositiveInt] | None = None

    ignore_composition: list[list[str]] | None = None

    _ignore_composition: dict[tuple[str], Literal[True]] | None = PrivateAttr(
        default=None
    )

    @model_validator(mode="after")
    def _build_rules_internal(self) -> Self:
        for mi, val in self.rules.items():
            for mj, rule in val.items():
                rl = rule.split()
                op = rl[0].lower()

                if mi == "solute" or mj == "solute":
                    raise ValueError(
                        "'solute' is not supported in 'rules' configuration."
                    )

                match op:
                    case RuleType.CM:
                        d = float(rl[1])
                        logger.trace(
                            f"Using {op} rule between {mi} and {mj} with dist={d}A"
                        )
                        self._rules[mi, mj] = CMRule(dist=d)
                    case RuleType.HB:
                        d = float(rl[rl.index("d") + 1]) if "d" in rl else 3.5
                        a = float(rl[rl.index("a") + 1]) if "a" in rl else 150.0
                        logger.trace(
                            f"Using {op} rule between {mi} and {mj} with dist={d}A and ang={a}º"
                        )
                        self._rules[mi, mj] = HBRule(dist=d, ang=a)

        return self

    @model_validator(mode="after")
    def _parse_solute(self) -> Self:
        if self.solute is not None:
            if "solute" in self.solute:
                raise ValueError("'solute' configuration cannot be solute.")

        errors = []

        attrs_solute_list = ("nucleus", "follow")

        for attr_nm in attrs_solute_list:
            attr = getattr(self, attr_nm)
            if attr is None:
                continue

            if "solute" in attr:
                if self.solute is None:
                    errors.append(
                        ValueError(
                            f"'solute' keywork used in {attr_nm}, but 'solute' not defined."
                        )
                    )

                id_s = attr.index("solute")
                attr.pop((id_s))
                attr[id_s:id_s] = self.solute
                setattr(self, attr_nm, list(set(attr)))

        attrs_solute_list_list = ("ignore_composition",)

        for attr_nm in attrs_solute_list_list:
            attr = getattr(self, attr_nm)
            if attr is None:
                continue

            for i, val in enumerate(attr):
                if "solute" in val:
                    if self.solute is None:
                        errors.append(
                            ValueError(
                                f"'solute' keywork used in {attr_nm}, but 'solute' not defined."
                            )
                        )

                    id_s = val.index("solute")
                    attr[i].pop((id_s))
                    attr[i][id_s:id_s] = self.solute
                    attr[i] = list(set(attr[i]))

        if errors:
            raise ExceptionGroup("Errors in input file", errors)

        return self

    @model_validator(mode="after")
    def _parse_ignore_composition(self) -> Self:
        if self.ignore_composition is not None:
            self._ignore_composition = {
                tuple(sorted(i)): True for i in self.ignore_composition
            }

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
