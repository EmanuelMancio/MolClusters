# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

import textwrap
from pathlib import Path

import pytest
from pydantic import ValidationError

from molclusters.config import CMRule, HBRule, MolClsConfig, read_config
from molclusters.output import FLUSH_THREADS


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

    def test_follow_solute_alone_does_not_warn(self, captured_logs: list[str]):
        MolClsConfig(rules={"A": {"A": "cm 5.0"}}, solute=["A"], follow=["solute"])

        assert captured_logs == []

    @pytest.mark.parametrize(
        ("follow", "unused"),
        [
            (["solute", 3], "[3]"),
            # without the keyword, listing the solute by name follows nothing
            (["A"], "['A']"),
        ],
    )
    def test_entries_without_effect_warn(
        self, captured_logs: list[str], follow: list, unused: str
    ):
        MolClsConfig(rules={"A": {"A": "cm 5.0"}}, solute=["A"], follow=follow)

        (warning,) = captured_logs
        assert f"'follow' entries {unused} have no effect" in warning


class TestSolvent:
    def test_defaults_to_every_non_solute_resname_in_the_rules(self):
        config = MolClsConfig(
            rules={"MOL": {"MOL": "cm 5", "SOL": "cm 5"}, "ION": {"SOL": "cm 5"}},
            solute=["MOL"],
        )

        assert config.solvent == ["ION", "SOL"]

    def test_given_solvent_is_kept(self):
        config = MolClsConfig(
            rules={"MOL": {"SOL": "cm 5", "ION": "cm 5"}},
            solute=["MOL"],
            solvent=["SOL"],
        )

        assert config.solvent == ["SOL"]

    def test_no_default_without_solute(self):
        assert MolClsConfig(rules={"A": {"B": "cm 5"}}).solvent is None


class TestDescribe:
    def test_shows_what_the_analysis_runs_with(self):
        config = MolClsConfig(
            rules={"B": {"A": "hb", "B": "cm 5"}},
            solute=["A"],
            solvent=["B"],
            nucleus=["solute"],
            follow=["solute"],
            ignore_composition=[["B", "A"], ["B"]],
            lammps_resnames={"A": "1-10", "B": [11, "12-20"]},
            lammps_timestep="2 fs",
            distance_backend="OpenMP",
            report_compression="gzip",
            flush_threads=16,
        )

        assert config.describe().splitlines() == [
            "rules (distances in angstrom, angles in degrees):",
            "  A - B: hb d 3.5 a 150.0 hmin 0.9 hmax 1.1 hq 0.3 aq -0.5",
            "  B - B: cm 5.0",
            "solute: A",
            "solvent: B",
            "nucleus: A",
            "follow: solute (one solute-<resid>.gro per solute)",
            "ignore_composition: A + B; B",
            "distance_backend: OpenMP (cm rules only)",
            "analyses: SizeEvolution, Lineage, SoluteSolvent, ClusterCoordinates, "
            "Nucleus, JsonReport",
            "report_compression: gzip (molclusters.jsonl.gz)",
            "flush_threads: 16 (files written at once)",
            "lammps_resnames: A = 1-10, B = 11-20",
            "lammps_timestep: 0.002 ps",
        ]

    def test_unset_options(self):
        config = MolClsConfig(rules={"A": {"A": "cm 5.0"}})

        assert config.describe().splitlines() == [
            "rules (distances in angstrom, angles in degrees):",
            "  A - A: cm 5.0",
            "solute: none",
            "solvent: none",
            "nucleus: none",
            "follow: none",
            "ignore_composition: none",
            "distance_backend: serial (cm rules only)",
            "analyses: SizeEvolution, Lineage, JsonReport",
            "report_compression: zstd (molclusters.jsonl.zst)",
            "flush_threads: 4 (files written at once)",
            "lammps_resnames: none",
            "lammps_timestep: none (LAMMPS dump times are step numbers)",
        ]

    def test_every_option_is_described(self):
        # a newly added config option must show up in the log too, even when unset
        described = MolClsConfig(rules={"A": {"A": "cm 5.0"}}).describe()

        for field in MolClsConfig.model_fields:
            assert any(line.startswith(f"{field}") for line in described.splitlines())


class TestAnalyses:
    def test_all_that_apply_run_by_default(self):
        config = MolClsConfig(rules={"A": {"A": "cm 5.0"}})

        assert config.analyses == {}

    def test_turned_off_ones_are_described(self):
        config = MolClsConfig(
            rules={"A": {"A": "cm 5.0"}},
            solute=["A"],
            analyses={"ClusterCoordinates": False, "JsonReport": False},
        )

        (line,) = [
            ln for ln in config.describe().splitlines() if ln.startswith("analyses")
        ]
        assert line == (
            "analyses: SizeEvolution, Lineage, SoluteSolvent "
            "(turned off: ClusterCoordinates, JsonReport)"
        )

    def test_unknown_names_raise(self):
        with pytest.raises(ValueError, match=r"Unknown analyses \['JSONReport'\]") as e:
            MolClsConfig(rules={"A": {"A": "cm 5.0"}}, analyses={"JSONReport": False})

        # the message lists the names to pick from
        assert "'JsonReport'" in str(e.value)

    def test_turned_on_without_what_they_need_raises(self):
        with pytest.raises(
            ValueError,
            match=r"\[\"SoluteSolvent needs 'solute'\", \"Nucleus needs 'nucleus'\"\]",
        ):
            MolClsConfig(
                rules={"A": {"A": "cm 5.0"}},
                analyses={"Nucleus": True, "SoluteSolvent": True, "JsonReport": True},
            )

    def test_read_from_yaml_the_way_a_user_would_write_it(self, tmp_path: Path):
        path = tmp_path / "config.yml"
        path.write_text(
            "rules:\n  A:\n    A: cm 5.0\nanalyses:\n  JsonReport: false\n"
            "  SizeEvolution: off\n"
        )

        config = read_config(path)

        assert config.analyses == {"JsonReport": False, "SizeEvolution": False}


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


class TestLammpsTimestep:
    def test_optional(self):
        config = MolClsConfig(rules={"A": {"A": "cm 5.0"}})

        assert config.lammps_timestep is None
        assert config.lammps_timestep_ps is None

    @pytest.mark.parametrize(
        ("spec", "ps"),
        [
            ("2 fs", 0.002),
            ("0.5fs", 0.0005),
            ("0.001 ps", 0.001),
            (" 1 NS ", 1000.0),
            ("1e-3 ps", 0.001),
        ],
    )
    def test_is_converted_to_ps(self, spec: str, ps: float):
        config = MolClsConfig(rules={"A": {"A": "cm 5.0"}}, lammps_timestep=spec)

        assert config.lammps_timestep_ps == pytest.approx(ps)

    @pytest.mark.parametrize("spec", [2, 0.002])
    def test_needs_a_unit(self, spec: float):
        with pytest.raises(ValueError, match="needs its unit"):
            MolClsConfig(rules={"A": {"A": "cm 5.0"}}, lammps_timestep=spec)

    @pytest.mark.parametrize("spec", ["2 min", "fs", "two fs", "2"])
    def test_invalid(self, spec: str):
        with pytest.raises(ValueError, match="expected a number and its unit"):
            MolClsConfig(rules={"A": {"A": "cm 5.0"}}, lammps_timestep=spec)

    @pytest.mark.parametrize("spec", ["0 fs", "-1 fs"])
    def test_must_be_positive(self, spec: str):
        with pytest.raises(ValueError, match="must be positive"):
            MolClsConfig(rules={"A": {"A": "cm 5.0"}}, lammps_timestep=spec)

    def test_read_from_yaml_the_way_a_user_would_write_it(self, tmp_path: Path):
        path = tmp_path / "input.yaml"
        path.write_text(
            textwrap.dedent("""\
                rules:
                  SOL:
                    SOL: cm 5.0
                lammps_timestep: 2 fs
                """)
        )

        assert read_config(path).lammps_timestep_ps == pytest.approx(0.002)


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

    def test_hb_rule_charge_and_mass_flags(self):
        spec = "hb aq -0.4 hq 0.2 hmax 2.1 hmin 0.5"
        rule = MolClsConfig(rules={"A": {"A": spec}})._rules["A", "A"]

        assert rule == HBRule(
            h_mass_min=0.5, h_mass_max=2.1, h_charge_min=0.2, a_charge_max=-0.4
        )
        assert str(rule) == "hb d 3.5 a 150.0 hmin 0.5 hmax 2.1 hq 0.2 aq -0.4"

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
            ("hb hq big", "Invalid number for 'hq'"),
            ("hb hmin 1.1", "hmin < hmax"),
            ("hb hmin 2 hmax 1", "hmin < hmax"),
        ],
    )
    def test_invalid_rules_raise(self, spec: str, message: str):
        with pytest.raises(ValidationError, match=message):
            MolClsConfig(rules={"A": {"B": spec}})

    @pytest.mark.parametrize(
        "spec", ["hb a 200", "hb d -1", "cm 0", "hb hmin -1", "hb aq nan", "hb hq inf"]
    )
    def test_out_of_range_values_raise(self, spec: str):
        with pytest.raises(ValidationError):
            MolClsConfig(rules={"A": {"B": spec}})

    def test_invalid_rule_message_names_the_rule(self):
        with pytest.raises(ValidationError, match="A:B"):
            MolClsConfig(rules={"A": {"B": "cm far"}})

    def test_conflicting_rules_for_a_pair_given_both_ways_raise(self):
        rules = {"A": {"B": "cm 5.0"}, "B": {"A": "hb"}}

        with pytest.raises(
            ValidationError,
            match=r"Conflicting rules A:B \('cm 5.0'\) and B:A \('hb'\)",
        ):
            MolClsConfig(rules=rules)

    def test_same_rule_for_a_pair_given_both_ways_is_accepted(self):
        rules = {"A": {"B": "hb"}, "B": {"A": "hb a 150 d 3.5"}}

        config = MolClsConfig(rules=rules)

        assert dict(config._rules) == {("A", "B"): HBRule(dist=3.5, ang=150.0)}

    def test_solute_keyword_not_allowed_in_rules(self):
        with pytest.raises(ValueError, match="not supported in 'rules'"):
            MolClsConfig(rules={"solute": {"A": "cm 1.0"}})


class TestDistanceBackend:
    def test_defaults_to_serial(self):
        config = MolClsConfig(rules={"A": {"A": "cm 5.0"}})

        assert config.distance_backend == "serial"

    def test_openmp_is_accepted(self):
        config = MolClsConfig(rules={"A": {"A": "cm 5.0"}}, distance_backend="OpenMP")

        assert config.distance_backend == "OpenMP"

    def test_unknown_backend_raises(self):
        with pytest.raises(ValidationError):
            MolClsConfig(rules={"A": {"A": "cm 5.0"}}, distance_backend="cuda")


class TestFlushThreads:
    def test_defaults_to_the_output_default(self):
        config = MolClsConfig(rules={"A": {"A": "cm 5.0"}})

        assert config.flush_threads == FLUSH_THREADS == 4

    def test_a_count_is_accepted(self):
        config = MolClsConfig(rules={"A": {"A": "cm 5.0"}}, flush_threads=16)

        assert config.flush_threads == 16

    @pytest.mark.parametrize("count", [0, -1, 1.5, "many"])
    def test_anything_but_a_positive_count_raises(self, count: object):
        with pytest.raises(ValidationError, match="flush_threads"):
            MolClsConfig(rules={"A": {"A": "cm 5.0"}}, flush_threads=count)


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

    def test_unknown_keys_warn(self, tmp_path: Path, captured_logs: list[str]):
        path = tmp_path / "input.json"
        path.write_text('{"rules": {"SOL": {"SOL": "cm 5.0"}}, "nucleous": ["SOL"]}')

        read_config(path)

        (warning,) = [m for m in captured_logs if m.startswith("Ignoring unknown")]
        assert warning.startswith("Ignoring unknown config key(s) ['nucleous']")

    def test_unsupported_extension(self, tmp_path: Path):
        path = tmp_path / "input.ini"
        path.write_text("")

        with pytest.raises(ValueError, match="Unsupported config type"):
            read_config(path)
