# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

import pytest

from molclusters.config import MolClsConfig


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
