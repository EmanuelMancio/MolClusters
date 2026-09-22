# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

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
