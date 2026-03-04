# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Logging tools and configuration."""

import functools
import sys
from enum import Enum
from pathlib import Path
from typing import Callable, Literal

from loguru import logger

from .version import version


class STATUS(Enum):
    """Status code, can be OFF (0) or ON (1)."""

    OFF = 0
    ON = 1


LOG_STATUS = STATUS.OFF


def start_logging(
    level: Literal["TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL"]
    | int = "INFO",
    *,
    filename: Path = Path("molclusters.log"),
) -> None:
    """Start log handlers for the molclusters library.

    Parameters
    ----------
    level : Literal[ "TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL" ] | int, optional
        Log level to start the logger, by default "INFO"
    """
    global LOG_STATUS

    if LOG_STATUS == STATUS.OFF:
        logger.remove(0)
        logger.add(sys.stderr, level=level)
        logger.add(filename.with_suffix(".json"), serialize=True, level=level)
        logger.add(filename, level=level)

        logger.enable("molclusters")
        logger.info(f"Package molclusters | Version {version} | Start logging.")

        LOG_STATUS = STATUS.ON


def _logger_wraps[**P, R](
    *,
    entry: bool = True,
    exit: bool = True,
    level: Literal["TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL"]
    | int = "DEBUG",
) -> Callable[P, R]:
    def wrapper(func: Callable[P, R]) -> Callable[P, R]:
        name = func.__name__

        @functools.wraps(func)
        def wrapped(*args: P.args, **kwargs: P.kwargs) -> R:
            global LOG_STATUS
            if LOG_STATUS == STATUS.OFF:
                return func(*args, **kwargs)

            logger_ = logger.opt(depth=1)

            if entry:
                logger_.log(
                    level, "Entering '{}' (args={}, kwargs={})", name, args, kwargs
                )

            try:
                result = func(*args, **kwargs)
            except Exception as e:
                logger_.log(
                    "ERROR", "Function '{}' raised an exception {}", name, str(e)
                )
                raise

            if exit:
                logger_.log(level, "Exiting '{}' (result={})", name, result)

            return result

        return wrapped

    return wrapper
