# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Provides the output files of a run: where they go and how they're declared.

Classes:
--------
- OutputFile: The declaration of a file (or family of files) an analysis writes.
- RunOutput: The service analyses write their files through during a run.

Constants:
----------
- COMPRESSED_SUFFIXES: The file suffixes `RunOutput.append` compresses, and how.
"""

import re
import zlib
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

import zstandard

# Measured on an 8.4 MB report: gzip 6 takes 5x as long as zstd 3 for a similar
# size, and higher levels of either cost far more time than they save space.
GZIP_LEVEL = 6
ZSTD_LEVEL = 3


class _GzipStream:
    """Compresses what's appended to a file into a single gzip member."""

    __slots__ = ["_zobj"]

    def __init__(self) -> None:
        # 16 + MAX_WBITS: the gzip header and trailer, not zlib's
        self._zobj = zlib.compressobj(GZIP_LEVEL, zlib.DEFLATED, 16 + zlib.MAX_WBITS)

    def compress(self, data: bytes) -> bytes:
        # a sync flush ends each chunk on a byte boundary, so what's been written
        # so far can be decompressed even if the stream is never finished
        return self._zobj.compress(data) + self._zobj.flush(zlib.Z_SYNC_FLUSH)

    def finish(self) -> bytes:
        return self._zobj.flush(zlib.Z_FINISH)


class _ZstdStream:
    """Compresses what's appended to a file into a single zstd frame."""

    __slots__ = ["_zobj"]

    def __init__(self) -> None:
        self._zobj = zstandard.ZstdCompressor(level=ZSTD_LEVEL).compressobj()

    def compress(self, data: bytes) -> bytes:
        # ending a block makes what's been written so far decompressible
        return self._zobj.compress(data) + self._zobj.flush(
            zstandard.COMPRESSOBJ_FLUSH_BLOCK
        )

    def finish(self) -> bytes:
        return self._zobj.flush(zstandard.COMPRESSOBJ_FLUSH_FINISH)


COMPRESSED_SUFFIXES: dict[str, type[_GzipStream | _ZstdStream]] = {
    ".gz": _GzipStream,
    ".zst": _ZstdStream,
}


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
    closed on exit (see `close`), even when the run is interrupted.

    Every file starts over in a new run: the first flush of a file overwrites
    it, so an earlier run's results are never mixed with this one's.

    Text is written in text mode, bytes as they are; a file takes one or the
    other. A file named with one of `COMPRESSED_SUFFIXES` (``.gz``, ``.zst``) is
    compressed as it's appended to, into one stream per run (text as UTF-8), which
    `close` completes. Every flush ends a compressed block, so a run killed before
    `close` still leaves everything flushed readable.

    Attributes
    ----------
    directory : Path
        The directory the files are written to.
    max_chars : int
        The buffered size, in characters (or bytes), that triggers a flush.
    """

    __slots__ = [
        "directory",
        "max_chars",
        "_pending",
        "_size",
        "_written",
        "_binary",
        "_streams",
    ]

    def __init__(self, directory: Path, max_chars: int = 64 * 2**20) -> None:
        """Initialize the output of a run.

        Parameters
        ----------
        directory : Path
            The directory the files are written to.
        max_chars : int
            The buffered size, in characters (or bytes), that triggers a flush.
        """
        self.directory = directory
        self.max_chars = max_chars
        self._pending: defaultdict[str, list[str | bytes]] = defaultdict(list)
        self._size = 0
        self._written: set[str] = set()
        # file name -> whether bytes (not text) are appended to it
        self._binary: dict[str, bool] = {}
        # compressed file name -> its compression stream, until `close`
        self._streams: dict[str, _GzipStream | _ZstdStream] = {}

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

    def append(self, name: str, data: str | bytes) -> None:
        """Queue `data` to be appended to the file `name`, flushing if the buffer is full.

        Parameters
        ----------
        name : str
            The file name; a compressed one ends with one of `COMPRESSED_SUFFIXES`.
        data : str | bytes
            The text, or bytes, to append.

        Raises
        ------
        TypeError
            If an uncompressed file is given text after bytes, or bytes after text.
        """
        binary = isinstance(data, bytes)
        if (
            self._binary.setdefault(name, binary) != binary
            and PurePosixPath(name).suffix not in COMPRESSED_SUFFIXES
        ):
            first, then = ("bytes", "text") if self._binary[name] else ("text", "bytes")
            raise TypeError(f"{name} was appended {first} to, and can't take {then}.")

        self._pending[name].append(data)
        self._size += len(data)
        if self._size >= self.max_chars:
            self.flush()

    def flush(self) -> None:
        """Append the queued data to its files, opening each file once.

        A file flushed for the first time is overwritten instead: whatever it held
        came from an earlier run.
        """
        for name, chunks in self._pending.items():
            mode = "a" if name in self._written else "w"
            stream = self.__stream(name)
            if stream is not None:
                data = b"".join(
                    chunk.encode() if isinstance(chunk, str) else chunk
                    for chunk in chunks
                )
                with self.path(name).open(mode + "b") as out:
                    out.write(stream.compress(data))
            else:
                with self.path(name).open(mode + "b" * self._binary[name]) as out:
                    out.writelines(chunks)
        self._pending.clear()
        self._size = 0

    def close(self) -> None:
        """Flush, and complete the compressed files' streams.

        Appending to a compressed file after this starts another stream after the
        first (a second gzip member or zstd frame), which readers of either format
        read on as one.
        """
        self.flush()
        for name, stream in self._streams.items():
            with self.path(name).open("ab") as out:
                out.write(stream.finish())
        self._streams.clear()

    def __stream(self, name: str) -> _GzipStream | _ZstdStream | None:
        """Give the compression stream of the file `name`, started on first use.

        Parameters
        ----------
        name : str
            The file name.

        Returns
        -------
        _GzipStream | _ZstdStream | None
            The stream, or None for a file that isn't compressed.
        """
        kind = COMPRESSED_SUFFIXES.get(PurePosixPath(name).suffix)
        if kind is None:
            return None
        if name not in self._streams:
            self._streams[name] = kind()
        return self._streams[name]

    def __enter__(self) -> "RunOutput":
        """Use the output for a run, closed when the run ends.

        Returns
        -------
        RunOutput
            This output.
        """
        return self

    def __exit__(self, *exc_info: object) -> None:
        """Close the output (see `close`), even if the run ended with an error."""
        self.close()
