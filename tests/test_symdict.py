# SPDX-FileCopyrightText: © 2026 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: LGPL-3.0-or-later

import pytest

from molclusters.symdict import SymmetricDict


@pytest.fixture
def sd() -> SymmetricDict[str, int]:
    d = SymmetricDict()
    d["b", "a"] = 1
    d["a", "c"] = 2
    d["c", "c"] = 3
    return d


def test_pair_access_is_order_independent(sd: SymmetricDict):
    assert sd["a", "b"] == sd["b", "a"] == 1
    assert ("b", "a") in sd
    assert ("a", "z") not in sd
    assert sd.get(("a", "z"), -1) == -1


def test_pairs_are_stored_once_in_normalized_order(sd: SymmetricDict):
    assert len(sd) == 3
    assert list(sd) == [("a", "b"), ("a", "c"), ("c", "c")]
    assert repr(sd) == "{('a', 'b'): 1, ('a', 'c'): 2, ('c', 'c'): 3}"


def test_overwriting_a_reversed_pair_keeps_one_entry(sd: SymmetricDict):
    sd["a", "b"] = 10

    assert len(sd) == 3
    assert sd["b", "a"] == 10


def test_single_key_lists_its_partners(sd: SymmetricDict):
    assert sd["a"] == {"b": 1, "c": 2}
    assert sd["c"] == {"a": 2, "c": 3}
    assert sd["z"] == {}


def test_single_key_assignment_sets_every_pair():
    sd = SymmetricDict()

    sd["x"] = {"y": 1, "x": 2}

    assert sd["y", "x"] == 1
    assert sd["x", "x"] == 2


def test_single_key_assignment_requires_a_dict():
    sd = SymmetricDict()

    with pytest.raises(ValueError, match="dict"):
        sd["x"] = 1


def test_delete_pair(sd: SymmetricDict):
    del sd["c", "a"]

    assert ("a", "c") not in sd
    assert len(sd) == 2


def test_delete_single_key_drops_all_its_pairs(sd: SymmetricDict):
    del sd["c"]

    assert list(sd) == [("a", "b")]


def test_all_keys(sd: SymmetricDict):
    assert sd.all_keys() == {"a", "b", "c"}
