# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Logging tools and configuration."""

import inspect
import logging
import platform
import sys
import warnings
from collections.abc import Callable
from enum import Enum
from pathlib import Path
from typing import Literal, TextIO

import MDAnalysis as mda
import networkx as nx
import numpy as np
import pandas as pd
from loguru import logger
from tqdm import tqdm

from .version import version

# The terminal only needs what happened; the log files keep loguru's full format,
# with the module, function and line of each message. `{analysis}` marks where the
# name of the analysis that logged a message goes (see `_formatter`).
_TERMINAL_FORMAT = "<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{analysis}{message}</level>"
_FILE_FORMAT = "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{analysis}{message}</level>"


def _formatter(template: str) -> Callable[[dict], str]:
    """Build a sink's format that names the analysis a message was logged by.

    A run tags what its analyses log with their name (``extra["analysis"]``, see
    `Run`), shown as a ``[Name]`` prefix; other messages have none.

    Parameters
    ----------
    template : str
        A loguru format with an ``{analysis}`` placeholder for the prefix.

    Returns
    -------
    Callable[[dict], str]
        The format, as a function of the record.
    """
    # a format function must end the line itself, which a format string doesn't
    tagged = template.replace("{analysis}", "[{extra[analysis]}] ") + "\n{exception}"
    untagged = template.replace("{analysis}", "") + "\n{exception}"
    return lambda record: tagged if record["extra"].get("analysis") else untagged


# Messages bound with `logger.bind(**FILE_ONLY)` only reach the log files: detail
# the terminal doesn't need, such as tracebacks or progress the bar already shows.
FILE_ONLY = {"file_only": True}


def _not_file_only(record: dict) -> bool:
    """Filter for the terminal sink.

    Returns
    -------
    bool
        False for messages bound with `FILE_ONLY`.
    """
    return not record["extra"].get("file_only", False)


_SHOWN_WARNINGS: set[tuple[type[Warning], str]] = set()


def _short_path(filename: str) -> str:
    """Shorten a source path to its package-relative part, to show warning origins.

    Returns
    -------
    str
        The path after the last ``site-packages``, or just the file name.
    """
    parts = Path(filename).parts
    if "site-packages" in parts:
        start = len(parts) - parts[::-1].index("site-packages")
        return "/".join(parts[start:])
    return Path(filename).name


def _showwarning(
    message: Warning | str,
    category: type[Warning],
    filename: str,
    lineno: int,
    file: TextIO | None = None,
    line: str | None = None,
) -> None:
    """Log a Python warning (e.g. MDAnalysis' "no dt information") through loguru.

    Replaces `warnings.showwarning`, so warnings reach the log files and don't break
    the progress bar; which warnings are shown is still up to the warning filters.
    Each message is only logged once: Python's filters repeat a warning for every
    line that triggers it, so MDAnalysis' "no dt information" (raised on every time
    lookup) would otherwise flood the log.
    """
    key = (category, str(message))
    if key in _SHOWN_WARNINGS:
        return
    _SHOWN_WARNINGS.add(key)

    logger.warning(
        "{} ({}:{}): {}", category.__name__, _short_path(filename), lineno, message
    )


class _InterceptHandler(logging.Handler):
    """Send standard-library log records (MDAnalysis' own logger) to loguru."""

    def emit(self, record: logging.LogRecord) -> None:
        try:
            level: str | int = logger.level(record.levelname).name
        except ValueError:
            level = record.levelno

        # find the caller that logged the message, so the log shows its origin
        frame, depth = inspect.currentframe(), 0
        while frame and (depth == 0 or frame.f_code.co_filename == logging.__file__):
            frame = frame.f_back
            depth += 1

        logger.opt(depth=depth, exception=record.exc_info).log(
            level, record.getMessage()
        )


def _route_third_party_logs() -> None:
    """Send Python warnings and stdlib log records at WARNING and up to loguru.

    MDAnalysis reports through both: `warnings.warn` (printed straight to stderr,
    bypassing the log files) and the `logging` module (silenced by the
    `NullHandler` MDAnalysis attaches to its logger).
    """
    warnings.showwarning = _showwarning
    logging.basicConfig(
        handlers=[_InterceptHandler()], level=logging.WARNING, force=True
    )


def format_duration(seconds: float) -> str:
    """Format a duration for humans, e.g. ``4.2s``, ``12m 03s`` or ``1h 02m 03s``.

    Returns
    -------
    str
        The formatted duration.
    """
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes, secs = divmod(round(seconds), 60)
    if minutes < 60:
        return f"{minutes}m {secs:02d}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h {minutes:02d}m {secs:02d}s"


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
    filename : Path | str, optional
        The log file, by default "molclusters.log"; a serialized copy is written
        next to it with a ``.json`` suffix. Warnings from other packages
        (MDAnalysis) are logged as well.
    """
    global LOG_STATUS

    if LOG_STATUS == STATUS.OFF:
        filename = Path(filename)

        logger.remove(0)
        # written through tqdm so a message doesn't break an active progress bar
        logger.add(
            lambda msg: tqdm.write(msg, end="", file=sys.stderr),
            level=level,
            format=_formatter(_TERMINAL_FORMAT),
            colorize=sys.stderr.isatty(),
            filter=_not_file_only,
        )
        logger.add(filename.with_suffix(".json"), serialize=True, level=level)
        logger.add(filename, level=level, format=_formatter(_FILE_FORMAT))
        _route_third_party_logs()

        logger.enable("molclusters")
        logger.info(f"MolClusters {version} | Logging to {filename.resolve()}")
        logger.info(
            f"Python {platform.python_version()}, MDAnalysis {mda.__version__}, "
            f"NumPy {np.__version__}, NetworkX {nx.__version__}, "
            f"pandas {pd.__version__}"
        )

        LOG_STATUS = STATUS.ON
