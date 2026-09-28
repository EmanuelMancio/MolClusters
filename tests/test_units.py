# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""The units of everything MolClusters reports, pinned against physical facts.

Each quantity is checked against a value known independently of how MolClusters
computes it (a model's literature value, a textbook case, a geometry built to
size), never against a copy of its formula, so a wrong conversion factor fails
here. The units are:

========================  ==============================================
Quantity                  Unit
========================  ==============================================
Mass                      amu (g/mol)
Charge                    elementary charges, e
Dipole moment             Debye, D
Radius, diameter, Rg      angstrom, Å
Volume                    Å³
Density                   g/cm³
Connection distance       Å (center-of-mass or donor-acceptor distance)
H-bond angle              degrees
Time, birth time, age     ps (for LAMMPS dumps, given `lammps_timestep`)
Death time, lifetime      ps
========================  ==============================================

The first classes check the properties; `TestOutputUnits` checks that a run
writes them to the report (molclusters.jsonl), the CSV tables and evo.txt in the
same units.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from MDAnalysis import Universe
from MDAnalysis.coordinates.memory import MemoryReader

import molclusters.main as cli
from molclusters.cluster import MolGroup
from molclusters.config import MolClsConfig
from molclusters.conntable import ConnectionTable
from molclusters.molclusters import MolClusters
from molclusters.report import read_report

# SPC/E water: O-H 1.0 A, H-O-H 109.47 deg, q(H) = +0.4238 e. Literature: 2.35 D.
SPCE_MASSES = [15.9994, 1.008, 1.008]
SPCE_CHARGES = [-0.8476, 0.4238, 0.4238]
SPCE_DIPOLE = 2.35
SPCE_MASS = sum(SPCE_MASSES)
_HALF_ANGLE = np.deg2rad(109.47 / 2)
# O at the origin, the H-O-H bisector along +y, so the dipole points along +y
SPCE_SITES = np.array(
    [
        [0.0, 0.0, 0.0],
        [np.sin(_HALF_ANGLE), np.cos(_HALF_ANGLE), 0.0],
        [-np.sin(_HALF_ANGLE), np.cos(_HALF_ANGLE), 0.0],
    ]
)
BOX = 50.0

SPHERE_SPACING = 1.0


def waters(
    frames: list[list[np.ndarray]], resnames: list[str], dt: float = 1.0
) -> Universe:
    """Build SPC/E waters, each given by its atoms' positions in every frame.

    Parameters
    ----------
    frames : list[list[np.ndarray]]
        Per frame, each water's (3, 3) atom positions (O, H, H), in angstroms.
    resnames : list[str]
        Each water's residue name.
    dt : float
        The time between frames, in ps.

    Returns
    -------
    Universe
        The in-memory Universe.
    """
    n_res = len(resnames)
    uni = Universe.empty(
        3 * n_res,
        n_residues=n_res,
        atom_resindex=np.repeat(np.arange(n_res), 3),
        trajectory=True,
    )
    uni.add_TopologyAttr("names", ["OW", "HW1", "HW2"] * n_res)
    uni.add_TopologyAttr("elements", ["O", "H", "H"] * n_res)
    uni.add_TopologyAttr("resnames", resnames)
    uni.add_TopologyAttr("resids", np.arange(1, n_res + 1))
    uni.add_TopologyAttr("masses", SPCE_MASSES * n_res)
    uni.add_TopologyAttr("charges", SPCE_CHARGES * n_res)
    uni.add_TopologyAttr(
        "bonds", [(3 * i, 3 * i + k) for i in range(n_res) for k in (1, 2)]
    )
    coords = np.array([np.concatenate(frame) for frame in frames], dtype=np.float32)
    uni.load_new(
        coords, format=MemoryReader, dimensions=[BOX, BOX, BOX, 90, 90, 90], dt=dt
    )
    # JsonReport records both paths
    uni.filename = "waters.top"
    uni.trajectory.filename = "waters.traj"
    return uni


def point_charges(positions: list, masses: list, charges: list, elements: list):
    """A molecule of atoms at the given positions, one residue.

    Returns
    -------
    MolGroup
        The group of the one molecule.
    """
    n_atoms = len(positions)
    uni = Universe.empty(
        n_atoms, n_residues=1, atom_resindex=np.zeros(n_atoms, int), trajectory=True
    )
    uni.add_TopologyAttr("resids", [1])
    uni.add_TopologyAttr("elements", elements)
    uni.add_TopologyAttr("masses", masses)
    uni.add_TopologyAttr("charges", charges)
    uni.add_TopologyAttr("bonds", [(i, i + 1) for i in range(n_atoms - 1)])
    uni.atoms.positions = np.asarray(positions) + BOX / 2
    uni.dimensions = [BOX, BOX, BOX, 90.0, 90.0, 90.0]
    return MolGroup(uni, [1])


def lone_atom(element: str) -> MolGroup:
    """A molecule of one atom.

    Returns
    -------
    MolGroup
        The group of the one molecule.
    """
    return point_charges([[0.0, 0.0, 0.0]], [1.0], [0.0], [element])


def grid_sphere(atom_mass: float = 1.0) -> MolGroup:
    """A molecule filling a sphere of radius 10 A with carbon atoms on a cubic grid.

    Each atom fills a cube of side SPHERE_SPACING, so the space the atoms' centers
    fill is the number of atoms times SPHERE_SPACING cubed (see `filled_radius`).

    Returns
    -------
    MolGroup
        The group of the one molecule.
    """
    axis = np.arange(-10.0, 10.0 + SPHERE_SPACING / 2, SPHERE_SPACING)
    grid = np.stack(np.meshgrid(axis, axis, axis), axis=-1).reshape(-1, 3)
    points = grid[np.linalg.norm(grid, axis=1) <= 10.0]
    n_atoms = len(points)
    return point_charges(
        points, [atom_mass] * n_atoms, [0.0] * n_atoms, ["C"] * n_atoms
    )


def filled_radius(sphere: MolGroup) -> float:
    """The radius of the sphere filling the space of a `grid_sphere`'s atom centers.

    Returns
    -------
    float
        The radius, in angstroms.
    """
    return (3 * len(sphere.atoms) * SPHERE_SPACING**3 / (4 * np.pi)) ** (1 / 3)


class TestPropertyUnits:
    def test_mass_is_in_amu(self):
        uni = waters([[SPCE_SITES + 10.0]], ["WAT"])

        assert MolGroup(uni, [1]).mass == pytest.approx(18.0154)

    def test_charge_is_in_elementary_charges(self):
        sodium = point_charges([[0.0, 0.0, 0.0]], [22.99], [1.0], ["Na"])

        assert sodium.charge == pytest.approx(1.0)

    def test_dipole_is_in_debye(self):
        # textbook: charges of +-1 e 1 A apart make 4.803 D
        pair = point_charges(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], [1.0, 1.0], [-1.0, 1.0], ["C", "C"]
        )

        assert pair.dipole_moment == pytest.approx(4.803, abs=1e-3)
        np.testing.assert_allclose(pair.dipole, [4.803, 0.0, 0.0], atol=1e-3)

    def test_dipole_of_spce_water_is_its_literature_value(self):
        uni = waters([[SPCE_SITES + 10.0]], ["WAT"])

        assert MolGroup(uni, [1]).dipole_moment == pytest.approx(SPCE_DIPOLE, abs=0.005)

    def test_radius_of_gyration_is_in_angstrom(self):
        # two equal masses 2 A apart: each is 1 A from the center of mass
        pair = point_charges(
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], [1.0, 1.0], [0.0, 0.0], ["C", "C"]
        )

        assert pair.radius_of_gyration == pytest.approx(1.0)

    def test_center_of_mass_is_in_angstrom(self):
        pair = point_charges(
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], [3.0, 1.0], [0.0, 0.0], ["C", "C"]
        )

        np.testing.assert_allclose(pair.center_of_mass, [BOX / 2 + 0.5] + [BOX / 2] * 2)

    def test_radius_diameter_and_volume_are_in_angstrom(self):
        sphere = grid_sphere()
        # half a carbon atom's van der Waals radius, 1.7 A
        radius = filled_radius(sphere) + 0.85

        assert sphere.radius == pytest.approx(radius, rel=1e-3)
        assert sphere.diameter == pytest.approx(2 * radius, rel=1e-3)
        assert sphere.volume == pytest.approx(4 / 3 * np.pi * radius**3, rel=3e-3)

    def test_lone_atom_has_half_its_van_der_waals_radius(self):
        assert lone_atom("O").radius == pytest.approx(1.52 / 2)
        assert lone_atom("Cl").radius == pytest.approx(1.75 / 2)

    def test_density_is_in_g_per_cm3(self):
        # liquid water: 18.015 amu per 29.915 A^3, i.e. 1 g/cm^3 within the space
        # the atoms fill, which the radius buffer then spreads over a larger sphere
        sphere = grid_sphere(atom_mass=18.015 / 29.915 * SPHERE_SPACING**3)

        diluted = (filled_radius(sphere) / sphere.radius) ** 3
        assert sphere.density == pytest.approx(diluted, rel=3e-3)

    def test_radius_needs_elements(self):
        uni = Universe.empty(1, n_residues=1, atom_resindex=[0], trajectory=True)
        uni.add_TopologyAttr("resids", [1])

        with pytest.raises(ValueError, match="guess_TopologyAttrs"):
            _ = MolGroup(uni, [1]).radius_buffer

    def test_radius_needs_known_elements(self):
        with pytest.raises(ValueError, match="Xx"):
            _ = lone_atom("Xx").radius_buffer


class TestConnectionUnits:
    def test_center_of_mass_distance_is_in_angstrom(self):
        uni = waters(
            [[SPCE_SITES + 10.0, SPCE_SITES + [13.0, 10.0, 10.0]]], ["WAT"] * 2
        )
        config = MolClsConfig(rules={"WAT": {"WAT": "cm 3.5"}})

        table = ConnectionTable(uni, config._rules, {"WAT": uni.atoms})

        assert table[1, 2] == pytest.approx(3.0, abs=1e-4)

    def test_hbond_distance_is_in_angstrom_and_angle_in_degrees(self):
        # water 1 points an H straight at water 2's O, 2.8 A from its own O; water
        # 2's Hs point away, so this is the only H-bond, and a linear one
        turn = _HALF_ANGLE - np.pi / 2  # brings the first H onto +x
        rotation = np.array([[np.cos(turn), -np.sin(turn), 0],
                             [np.sin(turn), np.cos(turn), 0],
                             [0, 0, 1]])  # fmt: skip
        donor = SPCE_SITES @ rotation.T
        assert donor[1] == pytest.approx([1.0, 0.0, 0.0], abs=1e-12)
        # (x, y) -> (y, -x): the bisector, and so both Hs, now point along +x
        acceptor = SPCE_SITES @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        assert acceptor[1:, 0].min() > 0
        uni = waters([[donor + 10.0, acceptor + [12.8, 10.0, 10.0]]], ["WAT"] * 2)
        config = MolClsConfig(rules={"WAT": {"WAT": "hb d 3.5 a 150"}})

        table = ConnectionTable(uni, config._rules, {"WAT": uni.atoms})

        (edge,) = table.conntab.edges(data=True)
        assert edge[2]["distance"] == pytest.approx(2.8, abs=1e-4)
        assert edge[2]["angle"] == pytest.approx(180.0, abs=1e-2)


class TestTimeUnits:
    APART = [SPCE_SITES + 10.0, SPCE_SITES + [30.0, 10.0, 10.0]]
    TOGETHER = [SPCE_SITES + 10.0, SPCE_SITES + [13.0, 10.0, 10.0]]

    def test_birth_time_and_age_are_in_ps(self, tmp_path: Path):
        # frame 0 the waters are apart; from frame 1 (t = 2 ps) they form a cluster
        apart, together = self.APART, self.TOGETHER
        uni = waters([apart, together, together, together], ["WAT"] * 2, dt=2.0)
        molcls = MolClusters(uni, MolClsConfig(rules={"WAT": {"WAT": "cm 3.5"}}))
        molcls.run(output_dir=tmp_path)

        (cluster,) = molcls.tracker.clusters.values()
        uni.trajectory[3]  # clusters describe the Universe's current frame
        assert uni.trajectory.time == pytest.approx(6.0)
        assert cluster.birth_time == pytest.approx(2.0)
        assert cluster.age == pytest.approx(4.0)

    def test_event_times_and_lifetimes_are_in_ps(self, tmp_path: Path):
        # the waters pair up at t = 2 ps and part at t = 6 ps: they lived 4 ps
        apart, together = self.APART, self.TOGETHER
        uni = waters([apart, together, together, apart], ["WAT"] * 2, dt=2.0)
        config = MolClsConfig(rules={"WAT": {"WAT": "cm 3.5"}})
        MolClusters(uni, config).run(output_dir=tmp_path)

        events = pd.read_csv(tmp_path / "cluster_events.csv")
        assert events["Event"].tolist() == ["formation", "dissolution"]
        assert events["Time"].tolist() == pytest.approx([2.0, 6.0])
        (life,) = pd.read_csv(tmp_path / "cluster_lifetimes.csv").itertuples()
        assert (life.BirthTime, life.DeathTime) == pytest.approx((2.0, 6.0))
        assert life.Lifetime == pytest.approx(4.0)
        frames = list(read_report(tmp_path / "molclusters.jsonl.zst").frames())
        (cluster,) = frames[2]["clusters"]
        assert cluster["birth_time"] == pytest.approx(2.0)


class TestLammpsTimeUnits:
    """LAMMPS dumps record step numbers; `lammps_timestep` makes their times ps."""

    DATA = Path(__file__).parent / "data" / "lammps_mini.data"

    def run_cli(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        *options: str,
        steps: tuple[int, ...] = (0, 500, 1000),
    ) -> Path:
        # a real-units run (2 fs steps), dumped at `steps`: by default every 500
        # steps, so at 0, 1 and 2 ps
        snapshot = Path(__file__).parent / "data" / "lammps_mini.lammpsdump"
        body = snapshot.read_text().split("\n", 2)[2]  # without the TIMESTEP item
        dump = tmp_path / "traj.lammpsdump"
        dump.write_text("".join(f"ITEM: TIMESTEP\n{step}\n{body}" for step in steps))
        config = tmp_path / "input.yaml"
        config.write_text(
            "rules:\n  SOL:\n    SOL: cm 5.0\nlammps_resnames:\n  SOL: 1-2\n"
            "lammps_timestep: 2 fs\n"
        )
        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr(cli, "start_logging", lambda **_: None)
        monkeypatch.setattr(
            sys, "argv",
            ["molclusters", str(dump), str(self.DATA), str(config), *options],
        )  # fmt: skip
        cli.main()
        return tmp_path

    @pytest.mark.parametrize(
        ("options", "times"),
        [
            ((), [0.0, 1.0, 2.0]),
            (("--traj-memory",), [0.0, 1.0, 2.0]),
            (("--traj-memory", "--in-memory-step", "2"), [0.0, 2.0]),
        ],
    )
    def test_times_are_in_ps(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        options: tuple,
        times: list[float],
    ):
        out = self.run_cli(tmp_path, monkeypatch, *options)

        frames = read_report(out / "molclusters.jsonl.zst").frames()
        assert [f["time"] for f in frames] == pytest.approx(times)
        np.testing.assert_allclose(np.loadtxt(out / "evo.txt", ndmin=2)[:, 0], times)

    @pytest.mark.parametrize(
        ("options", "times"),
        [
            ((), [2.0, 3.0, 5.0]),
            (("--traj-memory",), [2.0, 3.0, 5.0]),
            (("--traj-memory", "--in-memory-step", "2"), [2.0, 5.0]),
        ],
    )
    def test_a_continued_unevenly_dumped_run_keeps_its_times(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        options: tuple,
        times: list[float],
    ):
        # continued from step 1000, then dumped 500 and 1000 steps apart
        out = self.run_cli(tmp_path, monkeypatch, *options, steps=(1000, 1500, 2500))

        frames = read_report(out / "molclusters.jsonl.zst").frames()
        assert [f["time"] for f in frames] == pytest.approx(times)
        np.testing.assert_allclose(np.loadtxt(out / "evo.txt", ndmin=2)[:, 0], times)


class TestOutputUnits:
    """A run writes the properties in the units above.

    Two SPC/E waters with parallel dipoles and centers of mass 3 A apart form one
    cluster: WA is the solute and the nucleus, WB the solvent. Two frames, 2 ps
    apart.
    """

    @pytest.fixture
    def run(self, tmp_path: Path) -> tuple[Universe, Path]:
        pair = [SPCE_SITES + 10.0, SPCE_SITES + [13.0, 10.0, 10.0]]
        uni = waters([pair, pair], ["WA", "WB"], dt=2.0)
        config = MolClsConfig(
            rules={"WA": {"WB": "cm 3.5"}},
            solute=["WA"],
            solvent=["WB"],
            nucleus=["WA"],
        )
        MolClusters(uni, config).run(output_dir=tmp_path)
        return uni, tmp_path

    def test_report(self, run: tuple[Universe, Path]):
        uni, out = run
        report = read_report(out / "molclusters.jsonl.zst")
        frames = list(report.frames())
        cluster_group = MolGroup(uni, [1, 2])
        water = MolGroup(uni, [1])

        # the units the header declares are those the values are in
        assert report.header["units"] == {
            "time": "ps",
            "birth_time": "ps",
            "mass": "amu",
            "volume": "angstrom^3",
            "radius": "angstrom",
            "diameter": "angstrom",
            "density": "g/cm^3",
            "charge": "e",
            "dipole_moment": "D",
            "distance": "angstrom",
            "angle": "degrees",
        }
        assert [f["time"] for f in frames] == pytest.approx([0.0, 2.0])
        for frame in frames:
            (cluster,) = frame["clusters"]
            assert cluster["mass"] == pytest.approx(2 * SPCE_MASS)
            assert cluster["charge"] == pytest.approx(0.0, abs=1e-6)
            # parallel dipoles add up
            assert cluster["dipole_moment"] == pytest.approx(2 * SPCE_DIPOLE, abs=0.01)
            ((_, _, distance, _, _),) = cluster["connections"]
            assert distance == pytest.approx(3.0, abs=1e-4)
            # the geometric ones as the properties pinned above give them
            for key in ["radius", "diameter", "volume", "density"]:
                assert cluster[key] == pytest.approx(getattr(cluster_group, key))

            nuclei = cluster["Nucleus"]
            (nucleus,) = nuclei["nuclei"]
            assert nucleus["mass"] == pytest.approx(SPCE_MASS)
            assert nucleus["dipole_moment"] == pytest.approx(SPCE_DIPOLE, abs=0.005)
            combined = nuclei["combined_dipole_moment"]
            assert combined == pytest.approx(SPCE_DIPOLE, abs=0.005)
            assert nucleus["radius"] == pytest.approx(water.radius)
            assert nucleus["density"] == pytest.approx(water.density)

    def test_solute_solvent_table(self, run: tuple[Universe, Path]):
        uni, out = run
        table = pd.read_csv(out / "solute_solvent.csv")
        cluster_group = MolGroup(uni, [1, 2])

        assert list(table["Time"]) == pytest.approx([0.0, 2.0])
        assert list(table["Dipole"]) == pytest.approx([2 * SPCE_DIPOLE] * 2, abs=0.01)
        assert list(table["Charge"]) == pytest.approx([0.0] * 2, abs=1e-6)
        assert list(table["Radius"]) == pytest.approx([cluster_group.radius] * 2)
        assert list(table["Density"]) == pytest.approx([cluster_group.density] * 2)

    def test_nucleus_table(self, run: tuple[Universe, Path]):
        uni, out = run
        table = pd.read_csv(out / "nucleus_data.csv")
        water = MolGroup(uni, [1])

        assert list(table["Time"]) == pytest.approx([0.0, 2.0])
        assert list(table["Dipole"]) == pytest.approx([SPCE_DIPOLE] * 2, abs=0.005)
        assert list(table["Radius"]) == pytest.approx([water.radius] * 2)
        assert list(table["Density"]) == pytest.approx([water.density] * 2)

    def test_size_evolution_time(self, run: tuple[Universe, Path]):
        _, out = run
        evo = np.loadtxt(out / "evo.txt")

        np.testing.assert_allclose(evo[:, 0], [0.0, 2.0])
