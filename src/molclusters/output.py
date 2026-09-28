# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides the output files of a run: where they go and how they're declared.

Classes:
--------
- OutputFile: The declaration of a file (or family of files) an analysis writes.
- RunOutput: The service analyses write their files through during a run.
"""

import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path, PurePosixPath


@dataclass(frozen=True, slots=True)
class OutputFile:
    """A file, or a family of files, that an analysis writes.

    Attributes
    ----------
    name : str
        The file name, or a pattern where each ``<placeholder>`` stands for an
        integer, e.g. ``cls-id<id>.gro`` for ``cls-id1.gro``, ``cls-id2.gro``, ...
        It may start with folders, separated by ``/`` (placeholders only go in
        the file name), relative to the output directory.
    """

    name: str

    def matches(self, filename: str) -> bool:
        """Tell whether `filename` is this file, or one of this family of files.

        Parameters
        ----------
        filename : str
            A file's path relative to the output directory, with ``/`` separators.

        Returns
        -------
        bool
            True if `filename` is (one of) the declared file(s).
        """
        literal_parts = re.split(r"<[^<>]*>", self.name)
        pattern = r"\d+".join(re.escape(part) for part in literal_parts)
        return re.fullmatch(pattern, filename) is not None

    @property
    def is_pattern(self) -> bool:
        """Whether `name` is a pattern, standing for a family of files.

        Returns
        -------
        bool
            True if `name` has a ``<placeholder>``.
        """
        return re.search(r"<[^<>]*>", self.name) is not None

    def existing(self, directory: Path) -> list[Path]:
        """Find the files of this declaration that already exist in `directory`.

        Parameters
        ----------
        directory : Path
            The output directory.

        Returns
        -------
        list[Path]
            The existing files, sorted.
        """
        folder = directory / PurePosixPath(self.name).parent
        if not folder.is_dir():
            return []
        return sorted(
            file
            for file in folder.iterdir()
            if file.is_file() and self.matches(file.relative_to(directory).as_posix())
        )


class RunOutput:
    """Where a run's analyses write their files.

    Files written once, at the end of a run, are written to `path(name)`. Files
    appended to frame by frame go through `append`, which buffers them: opening a
    file costs milliseconds on Windows (antivirus scanning, on any file size),
    which dominated runs that appended every frame to their files. The buffer is
    shared by every analysis of the run, and used as a context manager it is
    flushed on exit, even when the run is interrupted.

    Every file starts over in a new run: the first flush of a file overwrites
    it, so an earlier run's results are never mixed with this one's.

    Attributes
    ----------
    directory : Path
        The directory the files are written to.
    max_chars : int
        The buffered size, in characters, that triggers a flush.
    """

    __slots__ = ["directory", "max_chars", "_pending", "_size", "_written"]

    def __init__(self, directory: Path, max_chars: int = 64 * 2**20) -> None:
        """Initialize the output of a run.

        Parameters
        ----------
        directory : Path
            The directory the files are written to.
        max_chars : int
            The buffered size, in characters, that triggers a flush.
        """
        self.directory = directory
        self.max_chars = max_chars
        self._pending: defaultdict[str, list[str]] = defaultdict(list)
        self._size = 0
        self._written: set[str] = set()

    @property
    def written(self) -> frozenset[str]:
        """The files written so far, named as given to `path` or `append`.

        Returns
        -------
        frozenset[str]
            The names of the files handed out by `path` or flushed by `flush`.
        """
        return frozenset(self._written)

    def path(self, name: str) -> Path:
        """Give the path to write the file `name` to.

        Parameters
        ----------
        name : str
            The file name, possibly in folders (see `OutputFile.name`).

        Returns
        -------
        Path
            The file's path, in `directory`; its folders are created if missing.
        """
        self._written.add(name)
        path = self.directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def append(self, name: str, text: str) -> None:
        """Queue `text` to be appended to the file `name`, flushing if the buffer is full.

        Parameters
        ----------
        name : str
            The file name.
        text : str
            The text to append.
        """
        self._pending[name].append(text)
        self._size += len(text)
        if self._size >= self.max_chars:
            self.flush()

    def flush(self) -> None:
        """Append the queued text to its files, opening each file once.

        A file flushed for the first time is overwritten instead: whatever it held
        came from an earlier run.
        """
        for name, chunks in self._pending.items():
            mode = "a" if name in self._written else "w"
            with self.path(name).open(mode) as out:
                out.writelines(chunks)
        self._pending.clear()
        self._size = 0

    def __enter__(self) -> "RunOutput":
        """Use the output for a run, flushed when the run ends.

        Returns
        -------
        RunOutput
            This output.
        """
        return self

    def __exit__(self, *exc_info: object) -> None:
        """Flush the buffered text, even if the run ended with an error."""
        self.flush()
