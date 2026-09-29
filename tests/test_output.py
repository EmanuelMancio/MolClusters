# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

import gzip
import io
import zlib
from pathlib import Path

import pytest
import zstandard

from molclusters.output import FLUSH_THREADS, OutputFile, RunOutput


class TestOutputFile:
    def test_a_plain_name_matches_only_itself(self):
        out = OutputFile("evo.txt")

        assert out.matches("evo.txt")
        assert not out.matches("evo.txt.bak")
        assert not out.matches("evoXtxt")  # the dot is literal

    @pytest.mark.parametrize("name", ["cls-id1.gro", "cls-id1042.gro"])
    def test_a_placeholder_matches_an_integer(self, name: str):
        assert OutputFile("cls-id<id>.gro").matches(name)

    @pytest.mark.parametrize(
        "name", ["cls-id.gro", "cls-idX.gro", "cls-id2-mine.gro", "cls-id2.gro.bak"]
    )
    def test_a_placeholder_matches_nothing_else(self, name: str):
        assert not OutputFile("cls-id<id>.gro").matches(name)

    def test_a_name_in_a_folder_matches_only_in_that_folder(self):
        out = OutputFile("gro/cls-id<id>.gro")

        assert out.matches("gro/cls-id3.gro")
        assert not out.matches("cls-id3.gro")
        assert not out.matches("other/cls-id3.gro")

    def test_a_placeholder_makes_it_a_pattern(self):
        assert OutputFile("cls-id<id>.gro").is_pattern
        assert not OutputFile("evo.txt").is_pattern

    def test_existing_finds_the_files_in_their_folder(self, tmp_path: Path):
        (tmp_path / "gro").mkdir()
        for name in ["gro/cls-id3.gro", "gro/cls-id10.gro", "gro/x.gro", "cls-id4.gro"]:
            (tmp_path / name).write_text("")
        (tmp_path / "gro/cls-id5.gro").mkdir()  # named like one, but not a file

        found = OutputFile("gro/cls-id<id>.gro").existing(tmp_path)

        assert found == [tmp_path / "gro/cls-id10.gro", tmp_path / "gro/cls-id3.gro"]
        assert OutputFile("cls-id<id>.gro").existing(tmp_path) == [
            tmp_path / "cls-id4.gro"
        ]
        assert OutputFile("none/cls-id<id>.gro").existing(tmp_path) == []
        # the folder is a file
        assert OutputFile("cls-id4.gro/cls-id<id>.gro").existing(tmp_path) == []


class TestRunOutput:
    def test_paths_are_in_the_output_directory(self, tmp_path: Path):
        assert RunOutput(tmp_path).path("evo.txt") == tmp_path / "evo.txt"

    def test_folders_of_a_path_are_created(self, tmp_path: Path):
        output = RunOutput(tmp_path)

        output.append("a/b/c.gro", "x\n")
        output.flush()

        assert (tmp_path / "a" / "b" / "c.gro").read_text() == "x\n"

    def test_a_folder_removed_after_its_first_file_is_not_made_again(
        self, tmp_path: Path
    ):
        # folders are created once per run, not once per file
        output = RunOutput(tmp_path)
        output.path("gro/a.gro")
        (tmp_path / "gro").rmdir()

        assert output.path("gro/b.gro") == tmp_path / "gro/b.gro"
        assert not (tmp_path / "gro").exists()

    def test_nothing_is_appended_before_a_flush(self, tmp_path: Path):
        output = RunOutput(tmp_path)

        output.append("a.gro", "frame 1\n")

        assert not (tmp_path / "a.gro").exists()

    def test_flush_appends_each_file_in_order(self, tmp_path: Path):
        output = RunOutput(tmp_path)

        for text in ["1\n", "2\n"]:
            output.append("a.gro", f"a{text}")
            output.append("b.gro", f"b{text}")
            output.flush()
        output.flush()  # nothing left to write

        assert (tmp_path / "a.gro").read_text() == "a1\na2\n"
        assert (tmp_path / "b.gro").read_text() == "b1\nb2\n"

    def test_the_first_flush_overwrites_an_earlier_runs_file(self, tmp_path: Path):
        (tmp_path / "a.gro").write_text("earlier run\n")
        output = RunOutput(tmp_path)

        output.append("a.gro", "1\n")
        output.flush()
        output.append("a.gro", "2\n")
        output.flush()

        assert (tmp_path / "a.gro").read_text() == "1\n2\n"

    def test_written_names_the_files_of_the_run(self, tmp_path: Path):
        output = RunOutput(tmp_path)

        output.path("whole.txt")
        output.append("gro/a.gro", "x")
        output.flush()
        output.append("unflushed.gro", "x")

        assert output.written == {"whole.txt", "gro/a.gro"}

    def test_a_full_buffer_flushes_itself(self, tmp_path: Path):
        output = RunOutput(tmp_path, max_chars=6)
        target = tmp_path / "a.gro"

        output.append("a.gro", "abc")
        assert not target.exists()
        output.append("a.gro", "def")
        assert target.read_text() == "abcdef"

        output.append("a.gro", "ghi")
        output.flush()
        assert target.read_text() == "abcdefghi"

    def test_exiting_the_context_flushes_even_on_errors(self, tmp_path: Path):
        with pytest.raises(KeyboardInterrupt), RunOutput(tmp_path) as output:
            output.append("a.gro", "frame\n")
            raise KeyboardInterrupt

        assert (tmp_path / "a.gro").read_text() == "frame\n"

    def test_bytes_are_written_as_they_are(self, tmp_path: Path):
        output = RunOutput(tmp_path)

        output.append("a.bin", b"1\n")
        output.append("a.bin", b"2\r\n")
        output.flush()

        assert (tmp_path / "a.bin").read_bytes() == b"1\n2\r\n"

    def test_many_files_are_flushed_side_by_side_each_in_order(self, tmp_path: Path):
        # more files than FLUSH_THREADS, of every kind, each appended over flushes
        names = [f"gro/{i}.gro" for i in range(3 * FLUSH_THREADS)]
        names += ["a.bin", "a.jsonl.zst", "a.jsonl.gz"]
        output = RunOutput(tmp_path)
        for frame in range(3):
            for name in names:
                text = f"{name} {frame}\n"
                output.append(name, text.encode() if name == "a.bin" else text)
            output.flush()
        output.close()

        for name in names:
            file = tmp_path / name
            if name.endswith((".zst", ".gz")):
                text = decompress(file.read_bytes(), file.suffix).decode()
            else:
                text = file.read_text()
            assert text == "".join(f"{name} {i}\n" for i in range(3))

    def test_a_write_that_fails_raises(self, tmp_path: Path):
        output = RunOutput(tmp_path)
        output.append("a.gro", "x")
        output.append("b.gro", "x")
        (tmp_path / "b.gro").mkdir()  # a folder can't be opened as a file

        with pytest.raises(OSError):  # PermissionError on Windows
            output.flush()

        assert (tmp_path / "a.gro").read_text() == "x"

    @pytest.mark.parametrize(
        ("first", "then"), [("text", b"bytes"), (b"bytes", "text")]
    )
    def test_a_file_takes_either_text_or_bytes(
        self, tmp_path: Path, first: str | bytes, then: str | bytes
    ):
        output = RunOutput(tmp_path)
        output.append("a.txt", first)

        with pytest.raises(TypeError, match="a.txt was appended"):
            output.append("a.txt", then)


def decompress(data: bytes, suffix: str, *, across_frames: bool = True) -> bytes:
    """Decompress a whole .gz or .zst file, every gzip member or zstd frame of it.

    Returns
    -------
    bytes
        The decompressed data (of the first zstd frame only, unless `across_frames`).
    """
    if suffix == ".gz":
        return gzip.decompress(data)
    reader = zstandard.ZstdDecompressor().stream_reader(
        io.BytesIO(data), read_across_frames=across_frames
    )
    return reader.read()


def decompress_unfinished(data: bytes, suffix: str) -> bytes:
    """Decompress what a stream that was never finished holds so far.

    Returns
    -------
    bytes
        The data flushed into the stream.
    """
    if suffix == ".gz":
        return zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(data)
    return zstandard.ZstdDecompressor().decompressobj().decompress(data)


@pytest.mark.parametrize("suffix", [".gz", ".zst"])
class TestCompressedOutput:
    def test_appends_are_compressed_into_one_stream(self, tmp_path: Path, suffix: str):
        name = f"report.jsonl{suffix}"
        with RunOutput(tmp_path) as output:
            for i in range(3):
                output.append(name, f"{i}\n".encode())
                output.flush()

        data = (tmp_path / name).read_bytes()
        assert decompress(data, suffix, across_frames=False) == b"0\n1\n2\n"
        if suffix == ".gz":
            assert data[:2] == b"\x1f\x8b"
        else:
            assert data[:4] == b"\x28\xb5\x2f\xfd"

    def test_text_is_compressed_as_utf8(self, tmp_path: Path, suffix: str):
        name = f"a.txt{suffix}"
        with RunOutput(tmp_path) as output:
            output.append(name, "1 Å\n")
            output.append(name, b"2\n")

        assert decompress((tmp_path / name).read_bytes(), suffix) == "1 Å\n2\n".encode()

    def test_what_was_flushed_is_readable_without_closing(
        self, tmp_path: Path, suffix: str
    ):
        # a killed run never closes its output
        name = f"a{suffix}"
        output = RunOutput(tmp_path)
        output.append(name, b"frame 0\n")
        output.flush()
        output.append(name, b"frame 1\n")
        output.flush()

        data = (tmp_path / name).read_bytes()
        assert decompress_unfinished(data, suffix) == b"frame 0\nframe 1\n"

    def test_the_first_flush_overwrites_an_earlier_runs_file(
        self, tmp_path: Path, suffix: str
    ):
        name = f"a{suffix}"
        with RunOutput(tmp_path) as output:
            output.append(name, b"earlier\n")
        with RunOutput(tmp_path) as output:
            output.append(name, b"later\n")

        assert decompress((tmp_path / name).read_bytes(), suffix) == b"later\n"

    def test_appending_after_closing_starts_another_stream(
        self, tmp_path: Path, suffix: str
    ):
        name = f"a{suffix}"
        output = RunOutput(tmp_path)
        output.append(name, b"1\n")
        output.close()
        output.append(name, b"2\n")
        output.close()

        assert decompress((tmp_path / name).read_bytes(), suffix) == b"1\n2\n"
