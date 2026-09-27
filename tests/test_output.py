# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

from pathlib import Path

import pytest

from molclusters.output import OutputFile, RunOutput


class TestOutputFile:
    def test_a_plain_name_matches_only_itself(self):
        out = OutputFile("evo.txt")

        assert out.matches("evo.txt")
        assert not out.matches("evo.txt.bak")
        assert not out.matches("evoXtxt")  # the dot is literal

    @pytest.mark.parametrize("name", ["cls-id1.gro", "cls-id1042.gro"])
    def test_a_placeholder_matches_an_integer(self, name: str):
        assert OutputFile("cls-id<id>.gro", append=True).matches(name)

    @pytest.mark.parametrize(
        "name", ["cls-id.gro", "cls-idX.gro", "cls-id2-mine.gro", "cls-id2.gro.bak"]
    )
    def test_a_placeholder_matches_nothing_else(self, name: str):
        assert not OutputFile("cls-id<id>.gro", append=True).matches(name)


class TestRunOutput:
    def test_paths_are_in_the_output_directory(self, tmp_path: Path):
        assert RunOutput(tmp_path).path("evo.txt") == tmp_path / "evo.txt"

    def test_nothing_is_appended_before_a_flush(self, tmp_path: Path):
        output = RunOutput(tmp_path)

        output.append("a.gro", "frame 1\n")

        assert not (tmp_path / "a.gro").exists()

    def test_flush_appends_each_file_in_order(self, tmp_path: Path):
        (tmp_path / "a.gro").write_text("earlier\n")
        output = RunOutput(tmp_path)

        for text in ["1\n", "2\n"]:
            output.append("a.gro", f"a{text}")
            output.append("b.gro", f"b{text}")
        output.flush()
        output.flush()  # nothing left to write

        assert (tmp_path / "a.gro").read_text() == "earlier\na1\na2\n"
        assert (tmp_path / "b.gro").read_text() == "b1\nb2\n"

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
