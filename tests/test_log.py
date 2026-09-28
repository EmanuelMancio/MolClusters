# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
import warnings

import pytest
from loguru import logger

from molclusters.log import (
    _FILE_FORMAT,
    _TERMINAL_FORMAT,
    FILE_ONLY,
    _formatter,
    _InterceptHandler,
    _not_file_only,
    _short_path,
    _showwarning,
    format_duration,
)


@pytest.mark.parametrize(
    ("seconds", "expected"),
    [(4.24, "4.2s"), (723, "12m 03s"), (3723, "1h 02m 03s"), (59.99, "60.0s")],
)
def test_format_duration(seconds: float, expected: str):
    assert format_duration(seconds) == expected


@pytest.mark.parametrize("template", [_TERMINAL_FORMAT, _FILE_FORMAT])
def test_messages_name_the_analysis_that_logged_them(template: str):
    logged: list[str] = []
    handler_id = logger.add(logged.append, format=_formatter(template))
    try:
        with logger.contextualize(analysis="SizeEvolution"):
            logger.info("summary")
        logger.info("run")
    finally:
        logger.remove(handler_id)

    tagged, untagged = logged
    assert tagged.endswith(" [SizeEvolution] summary\n")
    assert untagged.endswith(" run\n")
    assert "[" not in untagged


def test_file_only_messages_are_kept_off_the_terminal():
    terminal: list[str] = []
    handler_id = logger.add(terminal.append, format="{message}", filter=_not_file_only)
    try:
        logger.bind(**FILE_ONLY).info("progress")
        logger.info("summary")
    finally:
        logger.remove(handler_id)

    assert [m.strip() for m in terminal] == ["summary"]


@pytest.mark.parametrize(
    ("filename", "expected"),
    [
        ("/venv/lib/site-packages/MDAnalysis/coordinates/base.py", "MDAnalysis/coordinates/base.py"),
        ("C:\\scripts\\run.py", "run.py"),
    ],
)  # fmt: skip
def test_short_path(filename: str, expected: str):
    assert _short_path(filename) == expected


def test_python_warnings_are_logged(captured_logs: list[str]):
    with warnings.catch_warnings():
        warnings.simplefilter("always")
        warnings.showwarning = _showwarning
        warnings.warn("Reader has no dt information", UserWarning, stacklevel=1)

    (message,) = captured_logs
    assert message.startswith("UserWarning (test_log.py:")
    assert message.rstrip().endswith("): Reader has no dt information")


def test_repeated_warnings_are_logged_once(captured_logs: list[str]):
    with warnings.catch_warnings():
        warnings.simplefilter("always")
        warnings.showwarning = _showwarning
        warnings.warn("same message, first line", UserWarning, stacklevel=1)
        warnings.warn("same message, first line", UserWarning, stacklevel=1)
        warnings.warn("another message", UserWarning, stacklevel=1)

    assert len(captured_logs) == 2


def test_stdlib_log_records_are_logged(captured_logs: list[str]):
    std_logger = logging.getLogger("MDAnalysis.test")
    handler = _InterceptHandler()
    std_logger.addHandler(handler)
    try:
        std_logger.warning("Autocorrelation: step > 1")
    finally:
        std_logger.removeHandler(handler)

    assert [m.rstrip() for m in captured_logs] == ["Autocorrelation: step > 1"]
