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
from tqdm import tqdm

from .version import version

# The terminal only needs what happened; the log files keep loguru's full format,
# with the module, function and line of each message.
_TERMINAL_FORMAT = "<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>"


class STATUS(Enum):
    """Status code, can be OFF (0) or ON (1)."""

    OFF = 0
    ON = 1


LOG_STATUS = STATUS.OFF


def start_logging(
    level: Literal["TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL"]
    | int = "INFO",
    *,
    filename: Path | str = Path("molclusters.log"),
) -> None:
    """Start log handlers for the molclusters library.

    Parameters
    ----------
    level : Literal[ "TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL" ] | int, optional
        Log level to start the logger, by default "INFO"
    """
    global LOG_STATUS

    if LOG_STATUS == STATUS.OFF:
        filename = Path(filename)

        logger.remove(0)
        # written through tqdm so a message doesn't break an active progress bar
        logger.add(
            lambda msg: tqdm.write(msg, end="", file=sys.stderr),
            level=level,
            format=_TERMINAL_FORMAT,
            colorize=sys.stderr.isatty(),
        )
        logger.add(filename.with_suffix(".json"), serialize=True, level=level)
        logger.add(filename, level=level)

        logger.enable("molclusters")
        logger.info(f"MolClusters {version} | Logging to {filename.resolve()}")

        LOG_STATUS = STATUS.ON


def _logger_wraps[**P, R](
    func: Callable[P, R] | None = None,
    *,
    entry: bool = True,
    exit: bool = True,
    level: Literal["TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL"]
    | int = "DEBUG",
) -> Callable[P, R]:
    def wrapper(inner_func: Callable[P, R]) -> Callable[P, R]:
        name = inner_func.__name__

        @functools.wraps(inner_func)
        def wrapped(*args: P.args, **kwargs: P.kwargs) -> R:
            global LOG_STATUS
            if LOG_STATUS == STATUS.OFF:
                return inner_func(*args, **kwargs)

            logger_ = logger.opt(depth=1)

            if entry:
                logger_.log(
                    level, "Entering '{}' (args={}, kwargs={})", name, args, kwargs
                )

            try:
                result = inner_func(*args, **kwargs)
            except Exception as e:
                logger_.log(
                    "ERROR", "Function '{}' raised an exception {}", name, str(e)
                )
                raise

            if exit:
                logger_.log(level, "Exiting '{}' (result={})", name, result)

            return result

        return wrapped

    if func is None:
        return wrapper
    else:
        return wrapper(func)
