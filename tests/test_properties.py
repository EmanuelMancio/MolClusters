# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Every reported property, checked against an independent numpy reference.

The system is a compact cluster of charged and neutral three-site residues whose
whole (unwrapped) coordinates the test builds itself. The Universe only ever sees
those coordinates wrapped into a small periodic box, so a property that is not
computed on a correctly made-whole group, or that is disturbed by another group
being analyzed over the same atoms (a nucleus inside its cluster), drifts away
from the reference computed here straight from the whole coordinates.
"""

import json
from itertools import permutations
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd
import pytest
from MDAnalysis.coordinates.memory import MemoryReader

from molclusters.cluster import MolGroup
from molclusters.config import MolClsConfig
from molclusters.molclusters import MolClusters

BOX = 20.0
CUTOFF = 3.0
# 1 e·A in Debye, from SI: e * 1e-10 m / (1e-21 / c C·m)
E_ANGSTROM_IN_DEBYE = 1.602176634e-19 * 1e-10 / (1e-21 / 299792458)
# residue centers of mass sit on a grid SPACING apart: axis neighbours connect,
# diagonal ones (SPACING * sqrt(2) = 3.5 A) don't
SPACING = 2.5
GRID = (3, 2, 2)

# (masses, charges) of the three sites of each residue name
SPECIES = {
    "SOL": ([15.999, 1.008, 1.008], [-0.834, 0.417, 0.417]),
    "CAT": ([15.999, 1.008, 1.008], [-0.4, 0.7, 0.7]),  # +1
    "ANI": ([32.06, 15.999, 15.999], [0.2, -0.6, -0.6]),  # -1
}
# the first residue is solvent and the ions are spread out, so no nucleus shares
# the cluster's first residue and the ions only connect through the solvent
RESNAMES = ["SOL", "SOL", "CAT", "SOL", "SOL", "SOL",
            "SOL", "ANI", "SOL", "SOL", "SOL", "CAT"]  # fmt: skip
IONS = ["CAT", "ANI"]
RULES = {"SOL": {name: f"cm {CUTOFF}" for name in ["SOL", *IONS]}}

ALL = tuple(range(1, len(RESNAMES) + 1))
SOLVENT = tuple(r for r in ALL if RESNAMES[r - 1] == "SOL")
NUCLEI = tuple((r,) for r in ALL if RESNAMES[r - 1] in IONS)
ALL_IONS = tuple(r for (r,) in NUCLEI)

# cluster centers: in the middle of the box, and on a box corner, so that the
# cluster is split across all three periodic boundaries
CENTERED, ON_CORNER = 0, 1
CENTERS = [(BOX / 2, BOX / 2, BOX / 2), (0.3, 0.2, 0.1)]

# (JSON key, analyzer attribute) of every scalar property the outputs report
PROPERTIES = [
    ("Mass", "mass"),
    ("Charge", "charge"),
    ("Radius", "radius"),
    ("Diameter", "diameter"),
    ("Volume", "volume"),
    ("Density", "density"),
    ("Dipole Moment", "dipole_moment"),
    ("Sphericity", "sphericity"),
    ("Shape", "shape_parameter"),
]
# column of solute_solvent.csv / nucleus_data.csv -> JSON key
CSV_COLUMNS = {
    "Radius": "Radius",
    "Density": "Density",
    "Charge": "Charge",
    "Dipole": "Dipole Moment",
    "Spher": "Sphericity",
    "Shape": "Shape",
}


type Group = tuple[int, ...]
# the Universe, and each frame's whole coordinates
type System = tuple[mda.Universe, list[np.ndarray]]


def reference(uni: mda.Universe, whole: np.ndarray, resids: Group) -> dict:
    """Compute every property from first principles on the whole coordinates.

    Returns
    -------
    dict
        Property values keyed as in molclusters.json, plus ``"dipole"`` (vector)
        and ``"bsphere"`` (bounding-sphere radius).
    """
    atoms = uni.residues[np.array(resids) - 1].atoms
    m, q = atoms.masses, atoms.charges
    pos = whole[atoms.ix].astype(np.float64)

    mass = m.sum()
    rel = pos - m @ pos / mass
    rg = np.sqrt(m @ np.sum(rel**2, axis=1) / mass)
    moments = np.linalg.eigvalsh((rel * m[:, None]).T @ rel / mass)
    dev = moments - moments.mean()
    # the uniform sphere with that radius of gyration: Rg = sqrt(3/5) R
    volume = 4 / 3 * np.pi * (np.sqrt(5 / 3) * rg) ** 3
    dipole = q @ rel * E_ANGSTROM_IN_DEBYE

    return {
        "Mass": mass,
        "Charge": q.sum(),
        "Radius": rg,
        "Diameter": 2 * rg,
        "Volume": volume,
        # amu/A^3 -> g/cm^3: 1 amu = 1.66053906660e-24 g, 1 A^3 = 1e-24 cm^3
        "Density": mass / volume * 1.66053906660,
        "Dipole Moment": np.linalg.norm(dipole),
        "Sphericity": 1 - 1.5 * np.sum(dev**2) / moments.sum() ** 2,
        "Shape": 27 * np.prod(dev) / moments.sum() ** 3,
        "dipole": dipole,
        "bsphere": np.max(np.linalg.norm(pos - pos.mean(axis=0), axis=1)),
    }


def assert_matches(values: dict, expected: dict, label: object = "") -> None:
    """Compare JSON-keyed property values with the reference."""
    for key, _ in PROPERTIES:
        assert values[key] == pytest.approx(expected[key], rel=1e-4, abs=1e-4), (
            f"{key} of {label}"
        )


def assert_row_matches(row: pd.Series, groups: list[dict], frame: object) -> None:
    """Compare a CSV row with the reference averaged over `groups`."""
    for column, key in CSV_COLUMNS.items():
        expected = np.mean([group[key] for group in groups])
        assert row[column] == pytest.approx(expected, rel=1e-4, abs=1e-4), (
            f"{column}, frame {frame}"
        )


def analyzer_values(analyzer: MolGroup) -> dict:
    """Read every scalar property of `analyzer`, keyed as in molclusters.json.

    Returns
    -------
    dict
        The analyzer's property values.
    """
    return {key: getattr(analyzer, attr) for key, attr in PROPERTIES}


def build_whole_frame(rng: np.random.Generator, center: tuple) -> np.ndarray:
    """Place the residues on the grid around `center`, each randomly oriented.

    Returns
    -------
    np.ndarray
        Whole coordinates, one row per atom.
    """
    grid = np.array(np.meshgrid(*map(np.arange, GRID), indexing="ij"))
    nodes = grid.reshape(3, -1).T * SPACING
    nodes = nodes - nodes.mean(axis=0) + center

    frame = []
    for node, name in zip(nodes, RESNAMES, strict=True):
        masses = np.array(SPECIES[name][0])
        # first site in the middle, the other two ~1 A away in random directions
        offsets = np.zeros((3, 3))
        for k in (1, 2):
            direction = rng.normal(size=3)
            offsets[k] = rng.uniform(0.9, 1.1) * direction / np.linalg.norm(direction)
        offsets -= masses @ offsets / masses.sum()  # center of mass on the node
        frame.append(node + offsets)

    return np.concatenate(frame)


def wrap_by_residue(whole: np.ndarray) -> np.ndarray:
    """Wrap each residue into the box by its first atom, keeping residues whole.

    Returns
    -------
    np.ndarray
        Wrapped coordinates, like a ``-pbc mol`` trajectory.
    """
    first_atoms = np.repeat(whole[0::3], 3, axis=0)
    return whole - BOX * np.floor(first_atoms / BOX)


def make_universe(frames: list[np.ndarray]) -> mda.Universe:
    """Build a Universe of the three-site residues over the given frames.

    Returns
    -------
    mda.Universe
        The in-memory Universe.
    """
    n_res = len(RESNAMES)
    uni = mda.Universe.empty(
        3 * n_res,
        n_residues=n_res,
        atom_resindex=np.repeat(np.arange(n_res), 3),
        trajectory=True,
    )
    uni.add_TopologyAttr("names", ["A1", "A2", "A3"] * n_res)
    uni.add_TopologyAttr("elements", ["X"] * 3 * n_res)
    uni.add_TopologyAttr("resnames", RESNAMES)
    uni.add_TopologyAttr("resids", np.arange(1, n_res + 1))
    uni.add_TopologyAttr("masses", [m for n in RESNAMES for m in SPECIES[n][0]])
    uni.add_TopologyAttr("charges", [q for n in RESNAMES for q in SPECIES[n][1]])
    uni.add_TopologyAttr(
        "bonds", [(3 * i, 3 * i + k) for i in range(n_res) for k in (1, 2)]
    )
    uni.load_new(
        np.array(frames, dtype=np.float32),
        format=MemoryReader,
        dimensions=[BOX, BOX, BOX, 90.0, 90.0, 90.0],
        dt=1.0,
    )
    # JsonReport records both paths, and chokes on None
    uni.filename = "synthetic.top"
    uni.trajectory.filename = "synthetic.traj"
    return uni


@pytest.fixture
def system() -> System:
    """The cluster centered (frame 0) and split across a box corner (frame 1).

    Returns
    -------
    System
        The Universe, and each frame's whole coordinates (in float32, as stored).
    """
    rng = np.random.default_rng(2026)
    whole = [build_whole_frame(rng, center) for center in CENTERS]
    uni = make_universe([wrap_by_residue(w) for w in whole])
    return uni, [w.astype(np.float32) for w in whole]


GROUPS = {
    "cluster": ALL,
    "solvent": SOLVENT,
    "one-ion": NUCLEI[0],
    "all-ions": ALL_IONS,  # charged and not contiguous
}


class TestSystem:
    """The synthetic system is what the other tests assume it is."""

    def test_the_corner_frame_is_actually_split_across_the_box(self, system: System):
        uni, _ = system
        uni.trajectory[ON_CORNER]

        first_atoms = uni.atoms.positions[0::3]

        assert np.all(np.ptp(first_atoms, axis=0) > BOX / 2)
        assert np.all((first_atoms >= 0) & (first_atoms < BOX))

    def test_the_nuclei_are_charged_and_far_apart(self, system: System):
        uni, whole = system

        charges = [reference(uni, whole[0], nuc)["Charge"] for nuc in NUCLEI]
        dipole = reference(uni, whole[0], ALL_IONS)["Dipole Moment"]

        assert charges == pytest.approx([1.0, -1.0, 1.0])
        # the ions sit several angstrom apart: a large dipole is expected
        assert dipole > 5


class TestAnalyzerProperties:
    @pytest.mark.parametrize("frame", [CENTERED, ON_CORNER], ids=["centered", "corner"])
    @pytest.mark.parametrize("group", GROUPS.values(), ids=GROUPS.keys())
    def test_properties_match_the_reference(
        self, system: System, frame: int, group: Group
    ):
        uni, whole = system
        uni.trajectory[frame]
        expected = reference(uni, whole[frame], group)

        analyzer = MolGroup(uni, list(group))

        assert_matches(analyzer_values(analyzer), expected, group)
        np.testing.assert_allclose(analyzer.dipole, expected["dipole"], atol=1e-3)
        assert analyzer.bsphere[0] == pytest.approx(expected["bsphere"], rel=1e-4)

    def test_molecules_broken_across_the_box_are_made_whole(self, system: System):
        uni, whole = system
        uni.trajectory[ON_CORNER]
        # wrap atom by atom, splitting the residues themselves across the box
        uni.atoms.positions = whole[ON_CORNER] % BOX
        expected = reference(uni, whole[ON_CORNER], ALL)

        analyzer = MolGroup(uni, list(ALL))

        assert_matches(analyzer_values(analyzer), expected)

    @pytest.mark.parametrize(
        "order",
        list(permutations(range(len(NUCLEI) + 2))),
        ids=lambda order: "-".join(map(str, order)),
    )
    def test_overlapping_groups_read_in_any_order_keep_their_values(
        self, system: System, order: tuple[int, ...]
    ):
        # the cluster, each nucleus inside it and all nuclei combined, as the run
        # analyzes them, over the same atoms: reading one must not disturb another
        uni, whole = system
        uni.trajectory[ON_CORNER]
        groups = [ALL, *NUCLEI, ALL_IONS]
        analyzers = [MolGroup(uni, list(g)) for g in groups]

        values = {i: analyzer_values(analyzers[i]) for i in order}
        # and once more after everything was read
        for i, group in enumerate(groups):
            expected = reference(uni, whole[ON_CORNER], group)
            assert_matches(values[i], expected, group)
            assert_matches(analyzer_values(analyzers[i]), expected, group)

    def test_values_follow_the_frame(self, system: System):
        uni, whole = system
        analyzer = MolGroup(uni, list(ALL))

        for frame in [CENTERED, ON_CORNER, CENTERED]:
            uni.trajectory[frame]
            expected = reference(uni, whole[frame], ALL)
            assert_matches(analyzer_values(analyzer), expected, frame)


class TestRunOutputs:
    """The values written by a full run are the reference values."""

    @pytest.fixture
    def run_dir(
        self, system: System, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> Path:
        uni, _ = system
        config = MolClsConfig(rules=RULES, solute=IONS, solvent=["SOL"], nucleus=IONS)
        monkeypatch.chdir(tmp_path)
        MolClusters(uni, config).run()
        return tmp_path

    @pytest.fixture
    def json_frames(self, run_dir: Path) -> list[dict]:
        data = json.loads((run_dir / "molclusters.json").read_text())
        return data["MolClusters"]

    def test_the_run_sees_one_cluster_with_one_nucleus_per_ion(
        self, json_frames: list[dict]
    ):
        for frame in json_frames:
            (cluster,) = frame["Clusters"]
            assert tuple(cluster["ResIDs"]) == ALL
            nuclei = sorted(tuple(n["ResIDs"]) for n in cluster["Nucleus"])
            assert nuclei == list(NUCLEI)

    def test_cluster_properties(self, system: System, json_frames: list[dict]):
        uni, whole = system

        for f, frame in enumerate(json_frames):
            (cluster,) = frame["Clusters"]
            assert_matches(cluster, reference(uni, whole[f], ALL), f"frame {f}")

    def test_nucleus_properties(self, system: System, json_frames: list[dict]):
        uni, whole = system

        for f, frame in enumerate(json_frames):
            (cluster,) = frame["Clusters"]
            for nucleus in cluster["Nucleus"]:
                expected = reference(uni, whole[f], tuple(nucleus["ResIDs"]))
                assert_matches(nucleus, expected, f"frame {f}, {nucleus['ResIDs']}")

    def test_nuclei_dipole_is_the_dipole_of_all_nuclei_together(
        self, system: System, json_frames: list[dict]
    ):
        uni, whole = system

        for f, frame in enumerate(json_frames):
            (cluster,) = frame["Clusters"]
            expected = reference(uni, whole[f], ALL_IONS)["Dipole Moment"]
            assert cluster["NucleiDipole"] == pytest.approx(expected, rel=1e-4)

    def test_solute_solvent_table(self, system: System, run_dir: Path):
        uni, whole = system
        table = pd.read_csv(run_dir / "solute_solvent.csv")

        for f, row in table.iterrows():
            # the only cluster holds both solute and solvent
            assert_row_matches(row, [reference(uni, whole[f], ALL)], f)

    def test_nucleus_table_averages_the_nuclei(self, system: System, run_dir: Path):
        uni, whole = system
        table = pd.read_csv(run_dir / "nucleus_data.csv")

        for f, row in table.iterrows():
            nuclei = [reference(uni, whole[f], nuc) for nuc in NUCLEI]
            assert_row_matches(row, nuclei, f)
