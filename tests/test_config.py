# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

import pytest
from loguru import logger

from molclusters.config import MolClsConfig


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
