# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

from collections import Counter
from collections.abc import Callable, Sequence
from typing import Any

import pytest

from molclusters.config import MolClsConfig
from molclusters.tracker import ClusterTracker

from .conftest import ALL_PAIRS_RULES, CUTOFF, MOL_RULES, Groups, UniverseFactory

type Track = Callable[..., ClusterTracker]


def step(tracker: ClusterTracker, frame: int) -> None:
    """Advance to `frame` and reconcile clusters, as one iteration of a run does."""
    tracker.uni.trajectory[frame]
    tracker.update()


def snapshot(tracker: ClusterTracker) -> dict[int, set[int]]:
    """Cluster id -> set of resids, as plain ints for readable assertions.

    Returns
    -------
    dict[int, set[int]]
        The current clusters.
    """
    return {cid: {int(m) for m in cls} for cid, cls in tracker.clusters.items()}


def id_of(tracker: ClusterTracker, members: set[int]) -> int:
    """Return the id of the cluster made of exactly `members`.

    Returns
    -------
    int
        The matching cluster id.
    """
    matches = [cid for cid, mols in snapshot(tracker).items() if mols == members]
    assert len(matches) == 1, f"no single cluster {members} in {snapshot(tracker)}"
    return matches[0]


def id_of_in(clusters: dict[int, set[int]], members: set[int]) -> int:
    """Return the id of the cluster made of exactly `members` in a snapshot.

    Returns
    -------
    int
        The matching cluster id.
    """
    (cid,) = [cid for cid, mols in clusters.items() if mols == members]
    return cid


def assert_membership_is_consistent(tracker: ClusterTracker) -> None:
    """Check that `mol_clt` maps exactly the clustered molecules to their cluster."""
    clusters = snapshot(tracker)
    assert set(tracker.mol_clt) == set().union(*clusters.values())
    for cid, mols in clusters.items():
        assert {tracker.find(m) for m in mols} == {cid}


@pytest.fixture
def track(make_universe: UniverseFactory) -> Track:
    """Build a ClusterTracker over synthetic frames (frame 0 is loaded on creation).

    Returns
    -------
    Track
        ``track(frames, n_res, resnames=None, **config)``.
    """

    def factory(
        frames: Sequence[Groups],
        n_res: int,
        resnames: Sequence[str] | None = None,
        **config_kwargs: Any,  # noqa: ANN401
    ) -> ClusterTracker:
        uni = make_universe(frames, n_res, resnames)
        config_kwargs.setdefault("rules", MOL_RULES)
        return ClusterTracker(uni, MolClsConfig(**config_kwargs))

    return factory


class TestClusterIdentity:
    """The dominance algorithm: which cluster keeps its id from one frame to the next."""

    def test_initial_clusters_are_the_connected_components(self, track: Track):
        tracker = track([[[1, 2, 3], [4, 5]]], 6)

        assert sorted(snapshot(tracker).values(), key=len) == [{4, 5}, {1, 2, 3}]
        assert tracker.find(1) == id_of(tracker, {1, 2, 3})
        assert tracker.find(6) is None

    def test_unchanged_cluster_keeps_its_id(self, track: Track):
        tracker = track([[[1, 2, 3]], [[1, 2, 3]]], 4)
        cid = id_of(tracker, {1, 2, 3})

        step(tracker, 1)

        assert snapshot(tracker) == {cid: {1, 2, 3}}

    def test_growing_cluster_keeps_its_id(self, track: Track):
        tracker = track([[[1, 2, 3]], [[1, 2, 3, 4]]], 5)
        cid = id_of(tracker, {1, 2, 3})

        step(tracker, 1)

        assert snapshot(tracker) == {cid: {1, 2, 3, 4}}
        assert tracker.find(4) == cid

    def test_shrinking_cluster_keeps_its_id(self, track: Track):
        tracker = track([[[1, 2, 3, 4]], [[1, 2, 3]]], 4)
        cid = id_of(tracker, {1, 2, 3, 4})

        step(tracker, 1)

        assert snapshot(tracker) == {cid: {1, 2, 3}}
        assert tracker.find(4) is None

    def test_formation_creates_a_new_cluster(self, track: Track):
        tracker = track([[], [[1, 2]]], 3)
        assert tracker.clusters == {}

        step(tracker, 1)

        assert list(snapshot(tracker).values()) == [{1, 2}]

    def test_dissolution_removes_the_cluster(self, track: Track):
        tracker = track([[[1, 2]], []], 3)

        step(tracker, 1)

        assert tracker.clusters == {}
        assert tracker.mol_clt == {}

    def test_split_larger_fragment_keeps_the_id(self, track: Track):
        tracker = track([[[1, 2, 3, 4, 5]], [[1, 2, 3], [4, 5]]], 5)
        cid = id_of(tracker, {1, 2, 3, 4, 5})

        step(tracker, 1)

        assert id_of(tracker, {1, 2, 3}) == cid
        assert id_of(tracker, {4, 5}) > cid

    def test_merge_larger_cluster_keeps_its_id_even_if_younger(self, track: Track):
        tracker = track([[[1, 2]], [[1, 2], [3, 4, 5]], [[1, 2, 3, 4, 5]]], 5)
        step(tracker, 1)
        small, big = id_of(tracker, {1, 2}), id_of(tracker, {3, 4, 5})
        assert small < big

        step(tracker, 2)

        assert snapshot(tracker) == {big: {1, 2, 3, 4, 5}}

    def test_merge_of_equal_sizes_keeps_the_older_id(self, track: Track):
        tracker = track([[[4, 5, 6]], [[4, 5, 6], [1, 2, 3]], [[1, 2, 3, 4, 5, 6]]], 6)
        step(tracker, 1)
        older, younger = id_of(tracker, {4, 5, 6}), id_of(tracker, {1, 2, 3})
        assert older < younger

        step(tracker, 2)

        assert snapshot(tracker) == {older: {1, 2, 3, 4, 5, 6}}

    def test_three_way_merge_of_equal_sizes_keeps_the_oldest_id(self, track: Track):
        # the oldest cluster holds the highest resids, so it is not the first
        # candidate seen when the merged component is scanned
        tracker = track(
            [
                [[7, 8]],
                [[7, 8], [1, 2]],
                [[7, 8], [1, 2], [4, 5]],
                [[1, 2, 4, 5, 7, 8]],
            ],
            8,
        )
        step(tracker, 1)
        step(tracker, 2)
        oldest = id_of(tracker, {7, 8})
        middle, youngest = id_of(tracker, {1, 2}), id_of(tracker, {4, 5})
        assert oldest < middle < youngest

        step(tracker, 3)

        assert snapshot(tracker) == {oldest: {1, 2, 4, 5, 7, 8}}

    @pytest.mark.parametrize("order", [[6, 7, 5], [7, 5, 6], [5, 6, 7]])
    def test_older_cluster_is_the_lowest_id_among_all_top_ties(self, order: list[int]):
        origin = Counter({cid: 3 for cid in order})

        older = ClusterTracker._get_older_cluster(order, origin)

        assert older == 5

    def test_older_cluster_ignores_lower_ids_with_smaller_counts(self):
        origin = Counter({6: 4, 7: 4, 5: 2})

        older = ClusterTracker._get_older_cluster([6, 7, 5], origin)

        assert older == 6

    def test_dimer_bridging_two_clusters_is_always_new(self, track: Track):
        tracker = track([[[1, 2], [3, 4]], [[2, 3]]], 4)
        old_ids = set(tracker.clusters)

        step(tracker, 1)

        (new_id,) = tracker.clusters
        assert new_id not in old_ids
        assert snapshot(tracker)[new_id] == {2, 3}

    def test_one_molecule_from_each_cluster_is_new(self, track: Track):
        tracker = track([[[1, 2], [3, 4], [5, 6]], [[1, 3, 5]]], 6)
        old_ids = set(tracker.clusters)

        step(tracker, 1)

        (new_id,) = tracker.clusters
        assert new_id not in old_ids

    def test_fragment_loses_to_a_later_fragment_holding_more_of_the_cluster(
        self, track: Track
    ):
        # {1,2,10..14} is processed first (larger), but {3,4,5,6} carries more of
        # the original cluster, so the id must go to the latter.
        tracker = track(
            [[[1, 2, 3, 4, 5, 6]], [[1, 2, 10, 11, 12, 13, 14], [3, 4, 5, 6]]], 14
        )
        cid = id_of(tracker, {1, 2, 3, 4, 5, 6})

        step(tracker, 1)

        assert id_of(tracker, {3, 4, 5, 6}) == cid
        assert id_of(tracker, {1, 2, 10, 11, 12, 13, 14}) > cid

    def test_fragment_loses_to_a_later_dimer_fragment(self, track: Track):
        tracker = track([[[1, 2, 3, 4, 5]], [[1, 2, 10, 11, 12, 13], [3, 4]]], 13)
        cid = id_of(tracker, {1, 2, 3, 4, 5})

        step(tracker, 1)

        assert id_of(tracker, {3, 4}) == cid
        assert id_of(tracker, {1, 2, 10, 11, 12, 13}) > cid

    def test_unrelated_dimer_does_not_take_the_id(self, track: Track):
        # {20, 21} pair up elsewhere in the box, sharing nothing with the cluster
        tracker = track([[[1, 2, 3, 4, 5, 6]], [[1, 2, 3, 4, 5, 10], [20, 21]]], 21)
        cid = id_of(tracker, {1, 2, 3, 4, 5, 6})

        step(tracker, 1)

        assert id_of(tracker, {1, 2, 3, 4, 5, 10}) == cid
        assert id_of(tracker, {20, 21}) > cid

    def test_larger_fragment_keeps_the_id_over_a_dimer_fragment(self, track: Track):
        # the cluster splits into seven members (plus a newcomer) and a dimer; the
        # dimer holds less of the cluster, so it must not inherit the id
        tracker = track(
            [[[1, 2, 3, 4, 5, 6, 7, 8, 9]], [[1, 2, 3, 4, 5, 6, 7, 20], [8, 9]]], 20
        )
        cid = id_of(tracker, set(range(1, 10)))

        step(tracker, 1)

        assert id_of(tracker, {1, 2, 3, 4, 5, 6, 7, 20}) == cid
        assert id_of(tracker, {8, 9}) > cid

    def test_surviving_fragment_with_new_members_keeps_the_id(self, track: Track):
        # half the cluster dissolves while the other half picks up a newcomer
        tracker = track([[[1, 2, 3, 4, 5, 6]], [[1, 2, 3, 10]]], 10)
        cid = id_of(tracker, {1, 2, 3, 4, 5, 6})

        step(tracker, 1)

        assert snapshot(tracker) == {cid: {1, 2, 3, 10}}

    def test_cluster_already_continued_cannot_claim_a_second_fragment(
        self, track: Track
    ):
        tracker = track(
            [[[1, 2, 3, 4, 5, 6], [20, 21]], [[1, 2, 3, 4], [5, 6, 20]]], 21
        )
        big = id_of(tracker, {1, 2, 3, 4, 5, 6})
        small = id_of(tracker, {20, 21})

        step(tracker, 1)

        assert id_of(tracker, {1, 2, 3, 4}) == big
        assert id_of(tracker, {5, 6, 20}) not in {big, small}
        assert small not in tracker.clusters

    def test_ignored_composition_never_becomes_a_cluster(self, track: Track):
        resnames = ["MOL"] * 3 + ["SOL"] * 2
        tracker = track(
            [[[1, 2, 3], [4, 5]], [[1, 2, 3], [4, 5]]],
            5,
            resnames,
            rules={"MOL": {"MOL": f"cm {CUTOFF}"}, "SOL": {"SOL": f"cm {CUTOFF}"}},
            ignore_composition=[["SOL"]],
        )
        assert list(snapshot(tracker).values()) == [{1, 2, 3}]

        step(tracker, 1)

        assert list(snapshot(tracker).values()) == [{1, 2, 3}]
        assert tracker.find(4) is None

    def test_merge_loser_dies_whatever_the_group_order(self, track: Track):
        # A loses a merge to B in both cases, leaving a remnant elsewhere; only the
        # remnant's size relative to the merge group differs (free newcomers)
        a, b = list(range(1, 11)), list(range(11, 21))
        remnant_smaller = track(
            [[a, b], [list(range(1, 9)) + list(range(11, 20)), [9, 10, 30]]], 30
        )
        a_smaller = id_of(remnant_smaller, set(a))

        a, b = list(range(1, 11)), list(range(11, 18))
        remnant_larger = track(
            [[a, b], [[1, 2, 3, 4, *range(30, 42)], list(range(5, 11)) + b]], 41
        )
        a_larger = id_of(remnant_larger, set(a))

        step(remnant_smaller, 1)
        step(remnant_larger, 1)

        survives = {
            "remnant smaller": a_smaller in remnant_smaller.clusters,
            "remnant larger": a_larger in remnant_larger.clusters,
        }
        assert not any(survives.values()), f"A survives: {survives}"

    @pytest.mark.parametrize(
        ("frames", "n_res", "pure"),
        [
            pytest.param(
                [[[1, 2, 3, 4]], [[1, 2, 10], [3, 4]]], 10, {3, 4}, id="dimer"
            ),
            pytest.param(
                [[[1, 2, 3, 4, 5, 6]], [[1, 2, 3, 10, 11], [4, 5, 6]]],
                11,
                {4, 5, 6},
                id="trimer",
            ),
        ],
    )
    def test_pure_fragment_wins_a_tie_whatever_its_size(
        self, track: Track, frames: list[Groups], n_res: int, pure: set[int]
    ):
        # the cluster splits into two halves; one half picked up newcomers, the
        # other is made only of the cluster's own molecules and must keep the id
        tracker = track(frames, n_res)
        cid = id_of(tracker, set(frames[0][0]))

        step(tracker, 1)

        assert id_of(tracker, pure) == cid

    @pytest.mark.parametrize(
        "halves",
        [
            [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]],
            [[6, 7, 8, 9, 10], [1, 2, 3, 4, 5]],
            [[2, 4, 6, 8, 10], [1, 3, 5, 7, 9]],
        ],
        ids=["low-first", "high-first", "interleaved"],
    )
    def test_even_split_keeps_the_id_in_the_half_with_the_lowest_resid(
        self, track: Track, halves: list[list[int]]
    ):
        # pins the current tie-break: equal-size groups keep networkx's component
        # order (resid order), so the half holding resid 1 wins. Arbitrary rather
        # than physical; update if a spatial tie-break is adopted.
        tracker = track([[list(range(1, 11))], halves], 10)
        cid = id_of(tracker, set(range(1, 11)))

        step(tracker, 1)

        (low,) = [set(h) for h in halves if 1 in h]
        (high,) = [set(h) for h in halves if 1 not in h]
        assert id_of(tracker, low) == cid
        assert id_of(tracker, high) > cid

    @pytest.mark.parametrize(
        ("cluster", "pieces", "n_res"),
        [
            pytest.param(
                list(range(1, 10)),
                [[7, 8, 9], [1, 2, 3], [4, 5, 6]],
                9,
                id="three-way",
            ),
            pytest.param(
                list(range(1, 10)),
                [[2, 5, 8], [3, 6, 9], [1, 4, 7]],
                9,
                id="three-way-interleaved",
            ),
            pytest.param(
                list(range(1, 7)),
                [[4, 5, 6, 10], [1, 2, 3, 11]],
                11,
                id="halves-with-newcomers",
            ),
        ],
    )
    def test_even_split_into_equal_pieces_keeps_the_id_with_the_lowest_resid(
        self, track: Track, cluster: list[int], pieces: list[list[int]], n_res: int
    ):
        # same tie-break as above: every piece holds as many of the cluster's
        # molecules and is the same size (so equally pure), and the piece holding
        # resid 1 wins whatever its newcomers or listing order
        tracker = track([[cluster], pieces], n_res)
        cid = id_of(tracker, set(cluster))

        step(tracker, 1)

        ids = {min(p): id_of(tracker, set(p)) for p in pieces}
        assert ids.pop(1) == cid
        assert len(set(ids.values())) == len(ids)
        assert all(i > cid for i in ids.values())
        assert_membership_is_consistent(tracker)

    def test_one_frame_break_mints_a_transient_id(self, track: Track):
        # pins the current behaviour: no persistence window, so a contact broken for
        # a single frame creates an id that dies when the halves rejoin
        tracker = track(
            [[[1, 2, 3, 4, 5, 6]], [[1, 2, 3], [4, 5, 6]], [[1, 2, 3, 4, 5, 6]]], 6
        )
        cid = id_of(tracker, {1, 2, 3, 4, 5, 6})

        step(tracker, 1)
        transient = id_of(tracker, {4, 5, 6})
        assert transient > cid

        step(tracker, 2)

        assert snapshot(tracker) == {cid: {1, 2, 3, 4, 5, 6}}

    def test_three_way_split_keeps_the_id_in_the_largest_fragment(self, track: Track):
        tracker = track([[list(range(1, 10))], [[1, 2, 3, 4], [5, 6, 7], [8, 9]]], 9)
        cid = id_of(tracker, set(range(1, 10)))

        step(tracker, 1)

        assert id_of(tracker, {1, 2, 3, 4}) == cid
        trimer, dimer = id_of(tracker, {5, 6, 7}), id_of(tracker, {8, 9})
        assert cid < trimer < dimer
        assert_membership_is_consistent(tracker)

    def test_clusters_born_in_the_same_frame_are_numbered_largest_first(
        self, track: Track
    ):
        tracker = track([[], [[1, 2], [3, 4, 5, 6], [7, 8, 9]]], 9)

        step(tracker, 1)

        assert (
            id_of(tracker, {3, 4, 5, 6})
            < id_of(tracker, {7, 8, 9})
            < id_of(tracker, {1, 2})
        )

    def test_reformed_cluster_gets_a_new_id(self, track: Track):
        # no memory across a dissolution: the same molecules regrouping are a new
        # cluster, younger than the one that dissolved
        tracker = track([[[1, 2, 3]], [], [[1, 2, 3]]], 3)
        cid = id_of(tracker, {1, 2, 3})

        step(tracker, 1)
        step(tracker, 2)

        assert id_of(tracker, {1, 2, 3}) > cid

    def test_split_and_merge_in_one_frame_keeps_both_ids(self, track: Track):
        # a minority of A joins B while A's majority stays together
        tracker = track(
            [[[1, 2, 3, 4, 5, 6], [7, 8, 9, 10]], [[1, 2, 3, 4], [5, 6, 7, 8, 9, 10]]],
            10,
        )
        a, b = id_of(tracker, {1, 2, 3, 4, 5, 6}), id_of(tracker, {7, 8, 9, 10})

        step(tracker, 1)

        assert snapshot(tracker) == {a: {1, 2, 3, 4}, b: {5, 6, 7, 8, 9, 10}}
        assert_membership_is_consistent(tracker)

    def test_exchanging_molecules_keeps_both_ids(self, track: Track):
        tracker = track([[[1, 2, 3, 4], [5, 6, 7, 8]], [[1, 2, 3, 8], [4, 5, 6, 7]]], 8)
        a, b = id_of(tracker, {1, 2, 3, 4}), id_of(tracker, {5, 6, 7, 8})

        step(tracker, 1)

        assert snapshot(tracker) == {a: {1, 2, 3, 8}, b: {4, 5, 6, 7}}
        assert_membership_is_consistent(tracker)

    def test_group_goes_to_a_cluster_that_chose_it_over_a_bigger_contributor(
        self, track: Track
    ):
        # B gives 4 molecules to the first group but 5 to the second, so B continues
        # in the second and doesn't compete for the first, which A (3) keeps
        tracker = track(
            [
                [[1, 2, 3], list(range(4, 13))],
                [[1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12]],
            ],
            12,
        )
        a, b = id_of(tracker, {1, 2, 3}), id_of(tracker, set(range(4, 13)))

        step(tracker, 1)

        assert snapshot(tracker) == {a: {1, 2, 3, 4, 5, 6, 7}, b: {8, 9, 10, 11, 12}}
        assert_membership_is_consistent(tracker)

    def test_cluster_absorbed_one_molecule_at_a_time_dies(self, track: Track):
        # A's molecules join B and C one each (or go free), so A has no best group
        tracker = track(
            [[[1, 2, 3], [4, 5, 6], [7, 8, 9]], [[1, 4, 5, 6], [2, 7, 8, 9]]], 9
        )
        a = id_of(tracker, {1, 2, 3})
        b, c = id_of(tracker, {4, 5, 6}), id_of(tracker, {7, 8, 9})

        step(tracker, 1)

        assert snapshot(tracker) == {b: {1, 4, 5, 6}, c: {2, 7, 8, 9}}
        assert a not in tracker.clusters
        assert tracker.find(3) is None
        assert_membership_is_consistent(tracker)

    def test_cluster_scattered_into_other_clusters_dies(self, track: Track):
        # each of A's molecules joins a different cluster, none goes free
        tracker = track(
            [[[1, 2, 3], [4, 5], [6, 7], [8, 9]], [[1, 4, 5], [2, 6, 7], [3, 8, 9]]], 9
        )
        a = id_of(tracker, {1, 2, 3})
        b, c, d = (id_of(tracker, g) for g in ({4, 5}, {6, 7}, {8, 9}))

        step(tracker, 1)

        assert snapshot(tracker) == {b: {1, 4, 5}, c: {2, 6, 7}, d: {3, 8, 9}}
        assert a not in tracker.clusters
        assert_membership_is_consistent(tracker)

    def test_cluster_scattered_into_new_clusters_dies(self, track: Track):
        # each of A's molecules pairs with a free molecule, so no pair carries A
        tracker = track([[[1, 2, 3]], [[1, 10], [2, 11], [3, 12]]], 12)
        a = id_of(tracker, {1, 2, 3})

        step(tracker, 1)

        ids = [id_of(tracker, g) for g in ({1, 10}, {2, 11}, {3, 12})]
        assert a not in tracker.clusters
        assert len(set(ids)) == 3
        assert all(i > a for i in ids)
        assert_membership_is_consistent(tracker)

    def test_merge_loser_remnant_is_a_new_cluster(self, track: Track):
        tracker = track(
            [
                [list(range(1, 7)), list(range(7, 15))],
                [[1, 2, 3, 4, *range(7, 15)], [5, 6, 20]],
            ],
            20,
        )
        a, b = id_of(tracker, set(range(1, 7))), id_of(tracker, set(range(7, 15)))

        step(tracker, 1)

        assert id_of(tracker, {1, 2, 3, 4, *range(7, 15)}) == b
        assert id_of(tracker, {5, 6, 20}) > max(a, b)
        assert a not in tracker.clusters
        assert_membership_is_consistent(tracker)

    def test_ignored_composition_group_does_not_claim_the_id(self, track: Track):
        # the cluster's solvent shell {3, 4, 5} drifts off as a pure-solvent group,
        # which is ignored; it holds more of the cluster than {1, 2, 6}, but ignored
        # groups are invisible to the dominance check, so the id stays with the solute
        resnames = ["MOL"] + ["SOL"] * 4 + ["MOL"]
        tracker = track(
            [[[1, 2, 3, 4, 5]], [[1, 2, 6], [3, 4, 5]]],
            6,
            resnames,
            rules=ALL_PAIRS_RULES,
            ignore_composition=[["SOL"]],
        )
        cid = id_of(tracker, {1, 2, 3, 4, 5})

        step(tracker, 1)

        assert snapshot(tracker) == {cid: {1, 2, 6}}


def assert_flows_add_up(tracker: ClusterTracker, before: dict[int, set[int]]) -> None:
    """Check that the flows account for every molecule of both frames' clusters."""
    t = tracker.transition
    for cid, mols in before.items():
        assert sum(t.destinations(cid).values()) == len(mols)
    for cid, mols in snapshot(tracker).items():
        assert sum(t.sources(cid).values()) == len(mols)


class TestTransition:
    """What `ClusterTracker.transition` says happened between two frames."""

    def run(
        self,
        track: Track,
        frames: list[Groups],
        n_res: int,
        **config: Any,  # noqa: ANN401
    ) -> tuple[ClusterTracker, dict[int, set[int]], dict[int, set[int]]]:
        """Track two frames; give the tracker and each frame's clusters.

        Returns
        -------
        tuple[ClusterTracker, dict[int, set[int]], dict[int, set[int]]]
            The tracker, at the second frame, and both frames' clusters.
        """
        tracker = track(frames, n_res, **config)
        before = snapshot(tracker)
        step(tracker, 1)
        assert_flows_add_up(tracker, before)
        return tracker, before, snapshot(tracker)

    def test_first_frame_clusters_are_born_of_free_molecules(self, track: Track):
        tracker = track([[[1, 2, 3], [4, 5]]], 6)
        a, b = id_of(tracker, {1, 2, 3}), id_of(tracker, {4, 5})

        t = tracker.transition

        assert t.born == {a, b}
        assert dict(t.flows) == {(0, a): 3, (0, b): 2}
        assert (dict(t.merged), t.dissolved) == ({}, frozenset())

    def test_unchanged_cluster_flows_to_itself(self, track: Track):
        tracker, _, _ = self.run(track, [[[1, 2, 3]], [[1, 2, 3]]], 3)
        (cid,) = tracker.clusters

        t = tracker.transition

        assert dict(t.flows) == {(cid, cid): 3}
        assert (t.born, t.ended) == (frozenset(), frozenset())

    def test_growth_and_shrinking_flow_from_and_to_no_cluster(self, track: Track):
        tracker, _, _ = self.run(
            track, [[[1, 2, 3], [4, 5, 6]], [[1, 2, 3, 7], [4, 5]]], 7
        )
        a, b = id_of(tracker, {1, 2, 3, 7}), id_of(tracker, {4, 5})

        t = tracker.transition

        assert dict(t.flows) == {(a, a): 3, (0, a): 1, (b, b): 2, (b, 0): 1}
        assert t.sources(a) == {a: 3, 0: 1}
        assert t.destinations(b) == {b: 2, 0: 1}

    def test_formation_is_born_of_free_molecules(self, track: Track):
        tracker, _, _ = self.run(track, [[], [[1, 2]]], 3)
        (cid,) = tracker.clusters

        t = tracker.transition

        assert t.born == {cid}
        assert t.sources(cid) == {0: 2}

    def test_dissolution_ends_the_cluster(self, track: Track):
        tracker, before, _ = self.run(track, [[[1, 2]], []], 3)
        (cid,) = before

        t = tracker.transition

        assert t.dissolved == {cid}
        assert t.ended == {cid}
        assert t.destinations(cid) == {0: 2}

    def test_split_fragment_is_born_of_the_cluster(self, track: Track):
        tracker, before, _ = self.run(
            track, [[[1, 2, 3, 4, 5]], [[1, 2, 3], [4, 5]]], 5
        )
        (cid,) = before
        new = id_of(tracker, {4, 5})

        t = tracker.transition

        assert t.born == {new}
        assert t.sources(new) == {cid: 2}
        assert t.destinations(cid) == {cid: 3, new: 2}
        assert t.ended == frozenset()

    def test_split_into_many_gives_each_fragment_the_cluster_as_source(
        self, track: Track
    ):
        tracker, before, _ = self.run(
            track, [[list(range(1, 10))], [[1, 2, 3, 4], [5, 6, 7], [8, 9]]], 9
        )
        (cid,) = before
        trimer, dimer = id_of(tracker, {5, 6, 7}), id_of(tracker, {8, 9})

        t = tracker.transition

        assert t.born == {trimer, dimer}
        assert t.destinations(cid) == {cid: 4, trimer: 3, dimer: 2}

    def test_merge_maps_the_absorbed_cluster_to_the_survivor(self, track: Track):
        tracker, before, _ = self.run(
            track, [[[1, 2], [3, 4, 5]], [[1, 2, 3, 4, 5]]], 5
        )
        small, big = id_of_in(before, {1, 2}), id_of_in(before, {3, 4, 5})

        t = tracker.transition

        assert dict(t.merged) == {small: big}
        assert t.dissolved == frozenset()
        assert t.sources(big) == {big: 3, small: 2}
        assert t.absorbed(big) == [small]

    def test_merge_of_many_maps_every_absorbed_cluster_to_the_survivor(
        self, track: Track
    ):
        tracker, before, _ = self.run(
            track, [[[1, 2], [4, 5, 6], [7, 8]], [[1, 2, 4, 5, 6, 7, 8]]], 8
        )
        a, b, c = (id_of_in(before, g) for g in ({1, 2}, {4, 5, 6}, {7, 8}))

        t = tracker.transition

        assert dict(t.merged) == {a: b, c: b}
        assert t.absorbed(b) == sorted([a, c])
        assert t.sources(b) == {a: 2, b: 3, c: 2}
        assert t.ended == {a, c}

    def test_new_cluster_of_pieces_of_many_has_each_as_source(self, track: Track):
        # a trimer of one molecule from each cluster; the rest of each dissolves
        tracker, before, _ = self.run(track, [[[1, 2], [3, 4], [5, 6]], [[1, 3, 5]]], 6)
        (new,) = tracker.clusters

        t = tracker.transition

        assert t.sources(new) == {cid: 1 for cid in before}
        assert t.born == {new}
        assert t.dissolved == set(before)
        assert dict(t.merged) == {}

    def test_merge_and_split_of_several_clusters_in_one_frame(self, track: Track):
        # A (6) and B (4) merge, keeping B's id, with 4 of A's molecules; A's other
        # two leave with a newcomer, and C (3) splits into two new clusters
        tracker, before, _ = self.run(
            track,
            [
                [list(range(1, 7)), list(range(7, 15)), [20, 21, 22, 23]],
                [[1, 2, 3, 4, *range(7, 15)], [5, 6, 30], [20, 21], [22, 23]],
            ],
            30,
        )
        a = id_of_in(before, set(range(1, 7)))
        b = id_of_in(before, set(range(7, 15)))
        c = id_of_in(before, {20, 21, 22, 23})
        remnant = id_of(tracker, {5, 6, 30})

        t = tracker.transition

        assert dict(t.merged) == {a: b}
        assert t.sources(remnant) == {a: 2, 0: 1}
        assert id_of(tracker, {20, 21}) == c
        (piece,) = t.born - {remnant}
        assert t.destinations(c) == {c: 2, piece: 2}

    def test_cluster_absorbed_one_molecule_at_a_time_dissolves(self, track: Track):
        tracker, before, _ = self.run(
            track, [[[1, 2, 3], [4, 5, 6], [7, 8, 9]], [[1, 4, 5, 6], [2, 7, 8, 9]]], 9
        )
        a = id_of_in(before, {1, 2, 3})
        b, c = id_of(tracker, {1, 4, 5, 6}), id_of(tracker, {2, 7, 8, 9})

        t = tracker.transition

        assert t.dissolved == {a}
        assert t.destinations(a) == {b: 1, c: 1, 0: 1}

    def test_cluster_scattered_into_other_clusters_dissolves(self, track: Track):
        # its molecules all join other clusters, but none of them takes A in
        tracker, before, _ = self.run(
            track,
            [[[1, 2, 3], [4, 5], [6, 7], [8, 9]], [[1, 4, 5], [2, 6, 7], [3, 8, 9]]],
            9,
        )
        a = id_of_in(before, {1, 2, 3})
        b, c, d = (id_of_in(before, g) for g in ({4, 5}, {6, 7}, {8, 9}))

        t = tracker.transition

        assert t.dissolved == {a}
        assert dict(t.merged) == {}
        assert t.born == frozenset()
        assert t.destinations(a) == {b: 1, c: 1, d: 1}

    def test_cluster_scattered_into_new_clusters_dissolves(self, track: Track):
        tracker, before, _ = self.run(
            track, [[[1, 2, 3]], [[1, 10], [2, 11], [3, 12]]], 12
        )
        (a,) = before

        t = tracker.transition

        assert t.dissolved == {a}
        assert t.born == set(tracker.clusters)
        assert t.destinations(a) == {new: 1 for new in t.born}
        assert all(t.sources(new) == {a: 1, 0: 1} for new in t.born)

    def test_ignored_composition_counts_as_no_cluster(self, track: Track):
        resnames = ["MOL"] + ["SOL"] * 4 + ["MOL"]
        tracker, before, _ = self.run(
            track,
            [[[1, 2, 3, 4, 5]], [[1, 2, 6], [3, 4, 5]]],
            6,
            resnames=resnames,
            rules=ALL_PAIRS_RULES,
            ignore_composition=[["SOL"]],
        )
        (cid,) = before

        assert tracker.transition.destinations(cid) == {cid: 2, 0: 3}

    def test_is_read_only(self, track: Track):
        tracker = track([[[1, 2]]], 2)

        with pytest.raises(TypeError):
            tracker.transition.flows[0, 99] = 1  # type: ignore[index]
        with pytest.raises(AttributeError):
            tracker.transition.born = frozenset()  # type: ignore[misc]


class TestIds:
    def test_each_tracker_numbers_its_clusters_from_one(self, track: Track):
        # both clusters carry on in frame 1, and {6, 7} forms
        frames = [[[1, 2, 3], [4, 5]], [[1, 2], [3, 4, 5], [6, 7]]]

        first, second = track(frames, 7), track(frames, 7)
        for tracker in (first, second):
            step(tracker, 1)

        assert snapshot(first) == snapshot(second)
        assert snapshot(first) == {1: {1, 2}, 2: {3, 4, 5}, 3: {6, 7}}
