# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

import textwrap
from pathlib import Path

import pytest
from loguru import logger
from pydantic import ValidationError

from molclusters.config import CMRule, HBRule, MolClsConfig, read_config


@pytest.fixture
def captured_logs():
    # The package disables its logger by default (see `molclusters/__init__.py`);
    # `start_logging()` re-enables it for a real run. A dedicated sink is used
    # instead of capsys/capfd, since loguru's default handler binds its own
    # stderr reference ahead of pytest's output capture.
    messages: list[str] = []
    logger.enable("molclusters")
    handler_id = logger.add(messages.append, format="{message}")
    yield messages
    logger.remove(handler_id)
    logger.disable("molclusters")


class TestFollowSolute:
    def test_follow_solute_flag_set_and_expanded(self):
        config = MolClsConfig(
            rules={"A": {"A": "cm 5.0"}},
            solute=["A"],
            follow=["solute"],
        )

        assert config._follow_solute is True
        assert config.follow == ["A"]

    def test_follow_solute_flag_unset_without_keyword(self):
        config = MolClsConfig(
            rules={"A": {"A": "cm 5.0"}},
            solute=["A"],
            follow=["A"],
        )

        assert config._follow_solute is False
        assert config.follow == ["A"]

    def test_follow_solute_flag_unset_without_follow(self):
        config = MolClsConfig(rules={"A": {"A": "cm 5.0"}})

        assert config._follow_solute is False


class TestLammpsResnames:
    def test_no_mapping_configured(self):
        config = MolClsConfig(rules={"A": {"A": "cm 5.0"}})

        assert config.lammps_resnames is None
        assert config.resname_for_resid(1) is None

    def test_single_ids_and_ranges_are_expanded(self):
        config = MolClsConfig(
            rules={"SOL": {"SOL": "cm 5.0"}},
            lammps_resnames={"SOL": "1-3", "NA": [4, "5"], "CL": 6},
        )

        assert [config.resname_for_resid(i) for i in range(1, 7)] == [
            "SOL",
            "SOL",
            "SOL",
            "NA",
            "NA",
            "CL",
        ]
        assert config.resname_for_resid(7) is None

    def test_same_id_repeated_for_same_name_is_allowed(self):
        config = MolClsConfig(
            rules={"SOL": {"SOL": "cm 5.0"}},
            lammps_resnames={"SOL": ["1-3", 2]},
        )

        assert config.resname_for_resid(2) == "SOL"

    def test_same_id_assigned_to_two_names_raises(self):
        with pytest.raises(ValueError, match="assigned to both"):
            MolClsConfig(
                rules={"SOL": {"SOL": "cm 5.0"}},
                lammps_resnames={"SOL": "1-3", "NA": "3-4"},
            )

    def test_reversed_range_raises(self):
        with pytest.raises(ValueError, match="first id must not be greater"):
            MolClsConfig(
                rules={"SOL": {"SOL": "cm 5.0"}},
                lammps_resnames={"SOL": "5-1"},
            )

    def test_malformed_range_raises(self):
        with pytest.raises(ValueError, match="expected 'first-last'"):
            MolClsConfig(
                rules={"SOL": {"SOL": "cm 5.0"}},
                lammps_resnames={"SOL": "a-b"},
            )

    def test_non_numeric_id_raises(self):
        with pytest.raises(ValueError, match="expected a whole number"):
            MolClsConfig(
                rules={"SOL": {"SOL": "cm 5.0"}},
                lammps_resnames={"SOL": "abc"},
            )

    def test_nested_same_name_ranges_cover_the_full_outer_range(self):
        # A narrower range listed after a wider one for the same name must not
        # shadow the wider one's coverage.
        config = MolClsConfig(
            rules={"SOL": {"SOL": "cm 5.0"}},
            lammps_resnames={"SOL": ["1-100", "50-60"]},
        )

        assert config.resname_for_resid(70) == "SOL"

    def test_overlapping_same_name_ranges_warn_how_they_are_merged(
        self, captured_logs: list[str]
    ):
        MolClsConfig(
            rules={"SOL": {"SOL": "cm 5.0"}},
            lammps_resnames={"SOL": ["1-100", "50-60"]},
        )

        warning = "\n".join(captured_logs)
        assert "overlapping entries" in warning
        assert "50-60" in warning
        assert "single range, 1-100" in warning

    def test_adjacent_same_name_ranges_do_not_warn(self, captured_logs: list[str]):
        # Touching-but-not-overlapping ranges (e.g. two chunks "1-100"/"101-200"
        # of the same species) are a normal way to write a config, not a mistake.
        MolClsConfig(
            rules={"SOL": {"SOL": "cm 5.0"}},
            lammps_resnames={"SOL": ["1-100", "101-200"]},
        )

        assert captured_logs == []

    def test_large_range_is_stored_as_a_range_not_expanded(self):
        # A config spanning millions of molecule ids must stay a handful of
        # stored ranges, not one dict/list entry per id.
        config = MolClsConfig(
            rules={"SOL": {"SOL": "cm 5.0"}},
            lammps_resnames={"SOL": "1-2000000"},
        )

        assert len(config._lammps_resid_ranges) == 1
        assert config.resname_for_resid(1) == "SOL"
        assert config.resname_for_resid(1_500_000) == "SOL"
        assert config.resname_for_resid(2_000_001) is None


class TestReadConfigYaml:
    """Exercises `read_config` end to end, through an actual YAML file on disk.

    `MolClsConfig(...)` built directly from a Python dict (as the other tests in
    this module do) skips PyYAML's own type resolution, e.g. whether an unquoted
    ``1-3`` becomes a string, or a list mixing a quoted and an unquoted id
    survives intact. These tests write the YAML a user would actually type.
    """

    def test_reads_rules_and_solute(self, tmp_path: Path):
        path = tmp_path / "input.yaml"
        path.write_text(
            textwrap.dedent("""\
                rules:
                  SOL:
                    SOL: cm 5.0
                solute: [SOL]
                """)
        )

        config = read_config(path)

        assert config.solute == ["SOL"]
        assert config.rules == {"SOL": {"SOL": "cm 5.0"}}

    def test_reads_lammps_resnames_written_the_way_a_user_would(self, tmp_path: Path):
        # Unquoted range, a plain int, and a list mixing a quoted and an
        # unquoted id -- exactly as someone editing the file by hand would
        # write it, with no reason to know PyYAML would parse them differently.
        path = tmp_path / "input.yaml"
        path.write_text(
            textwrap.dedent("""\
                rules:
                  SOL:
                    SOL: cm 5.0
                lammps_resnames:
                  SOL: 1-3
                  NA: 4
                  CL: ["5", 6]
                """)
        )

        config = read_config(path)

        assert [config.resname_for_resid(i) for i in range(1, 7)] == [
            "SOL",
            "SOL",
            "SOL",
            "NA",
            "CL",
            "CL",
        ]

    def test_invalid_lammps_range_from_yaml_raises(self, tmp_path: Path):
        path = tmp_path / "input.yaml"
        path.write_text(
            textwrap.dedent("""\
                rules:
                  SOL:
                    SOL: cm 5.0
                lammps_resnames:
                  SOL: 5-1
                """)
        )

        with pytest.raises(ValueError, match="first id must not be greater"):
            read_config(path)


class TestRules:
    def test_cm_rule(self):
        rule = MolClsConfig(rules={"A": {"B": "cm 4.5"}})._rules["B", "A"]

        assert rule == CMRule(dist=4.5)

    def test_hb_rule_defaults(self):
        rule = MolClsConfig(rules={"A": {"A": "hb"}})._rules["A", "A"]

        assert rule == HBRule(dist=3.5, ang=150.0)

    def test_hb_rule_flags_in_any_order(self):
        rule = MolClsConfig(rules={"A": {"A": "hb a 120 d 3.0"}})._rules["A", "A"]

        assert rule == HBRule(dist=3.0, ang=120.0)

    @pytest.mark.parametrize(
        ("spec", "message"),
        [
            ("", "Empty rule"),
            ("xx 1.0", "either 'cm' or 'hb'"),
            ("cm", "cm <number>"),
            ("cm 1 2", "cm <number>"),
            ("cm far", "Invalid number for 'cm'"),
            ("hb d", "pairs"),
            ("hb d 3 d 4", "only be used once"),
            ("hb x 3", "only supports flags"),
            ("hb a wide", "Invalid number for 'a'"),
        ],
    )
    def test_invalid_rules_raise(self, spec: str, message: str):
        with pytest.raises(ValidationError, match=message) as err:
            MolClsConfig(rules={"A": {"B": spec}})

        # the note is attached to the original error, which pydantic wraps
        original = err.value.errors()[0]["ctx"]["error"]
        assert original.__notes__ == ["From rule A:B"]

    @pytest.mark.parametrize("spec", ["hb a 200", "hb d -1", "cm 0"])
    def test_out_of_range_values_raise(self, spec: str):
        with pytest.raises(ValidationError):
            MolClsConfig(rules={"A": {"B": spec}})

    @pytest.mark.xfail(
        strict=True,
        reason="pydantic drops the 'From rule A:B' note from the message it prints",
    )
    def test_invalid_rule_message_names_the_rule(self):
        with pytest.raises(ValidationError, match="A:B"):
            MolClsConfig(rules={"A": {"B": "cm far"}})

    def test_solute_keyword_not_allowed_in_rules(self):
        with pytest.raises(ValueError, match="not supported in 'rules'"):
            MolClsConfig(rules={"solute": {"A": "cm 1.0"}})


class TestSoluteKeyword:
    def test_solute_cannot_be_named_solute(self):
        with pytest.raises(ValueError, match="cannot be solute"):
            MolClsConfig(rules={"A": {"A": "cm 1.0"}}, solute=["solute"])

    def test_nucleus_is_expanded(self):
        config = MolClsConfig(
            rules={"A": {"A": "cm 1.0"}}, solute=["A", "B"], nucleus=["solute", "C"]
        )

        assert sorted(config.nucleus) == ["A", "B", "C"]

    def test_keyword_without_solute_defined_raises(self):
        with pytest.raises(ExceptionGroup) as err:
            MolClsConfig(
                rules={"A": {"A": "cm 1.0"}},
                nucleus=["solute"],
                ignore_composition=[["solute"]],
            )

        assert len(err.value.exceptions) == 2

    def test_ignore_composition(self):
        config = MolClsConfig(
            rules={"A": {"A": "cm 1.0"}},
            solute=["A"],
            ignore_composition=[["B"], ["solute", "B"]],
        )

        assert config.is_ignored_composition(["B", "B"])
        assert config.is_ignored_composition(["B", "A", "B"])
        assert not config.is_ignored_composition(["A"])


class TestReadConfigFormats:
    def test_json(self, tmp_path: Path):
        path = tmp_path / "input.json"
        path.write_text('{"rules": {"SOL": {"SOL": "cm 5.0"}}, "nucleus": ["SOL"]}')

        config = read_config(path)

        assert config.nucleus == ["SOL"]

    def test_toml(self, tmp_path: Path):
        path = tmp_path / "input.toml"
        path.write_text('solute = ["SOL"]\n[rules.SOL]\nSOL = "hb d 3.0"\n')

        config = read_config(str(path))

        assert config.solute == ["SOL"]
        assert config._rules["SOL", "SOL"] == HBRule(dist=3.0)

    def test_unsupported_extension(self, tmp_path: Path):
        path = tmp_path / "input.ini"
        path.write_text("")

        with pytest.raises(ValueError, match="Unsupported config type"):
            read_config(path)
