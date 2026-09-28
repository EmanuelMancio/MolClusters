# The tracking algorithm

This page describes `ClusterTracker.update` (`tracker.py`) in enough detail to
change it. The user-facing summary is [How clusters are tracked](../user-guide/concepts.md).

## Terms

- A **connected group** (`subconn` in the code, a `ConnectionTable._SubConnTable`)
  is a connected component of the current frame's connection graph: molecules linked
  to each other, directly or indirectly, and to nothing else. It has no id yet.
- A **cluster** is a connected group once `update` gives it an id, new or inherited.
- `mol_clt` maps each molecule in a cluster to its cluster id. During `update` it
  still holds the *previous* frame's assignment until it is overwritten.

## Steps

1. `conntab.update()` rebuilds the connection graph for the current frame.
2. For each connected group whose composition isn't ignored,
   `_count_origins` counts where its molecules came from: previous cluster id →
   number of molecules, id `0` for molecules that were free. For example, 4
   molecules of cluster 7, 2 of cluster 9 and 1 free molecule give
   `Counter({7: 4, 9: 2, 0: 1})`.
3. **Each previous cluster picks its best group** (`_best_groups`): the group holding
   the most of its molecules; on a tie the smallest group (the purest, largest share
   of its molecules from the cluster); on a further tie the first group in the
   frame's order. Contributions of a single molecule are ignored, so a cluster
   with no contribution of two or more molecules picks nothing and dies. Every group
   is looked at before any choice is made, so the result doesn't depend on the order
   groups are processed in.
4. **Each group continues one of the clusters that picked it**
   (`_get_older_cluster`): the one that contributed the most molecules; on a tie the
   oldest (lowest id, since ids are handed out in increasing order). A group nobody
   picked becomes a new cluster (`_create_new_cluster`) with the next id.
5. The continuing cluster gets the group's frozen graph (`Cluster._update`);
   `mol_clt` is rewritten for the group's molecules. Molecules in no group are
   dropped from `mol_clt` (a flow to id 0), and previous clusters no group continues
   are removed.
6. A `Transition` records the frame: `flows[(prev, cur)]` molecule counts (0 for no
   cluster), `born`, `merged` (every other cluster that picked a group → the one that
   continued in it) and `dissolved` (ended without merging).

## Resulting events

| Event | How it arises |
| --- | --- |
| formation | a group of free molecules is picked by nobody → new |
| growth, shrinking | a group picked by one cluster continues it |
| split | the cluster continues in its best fragment; the other fragments, picked by nobody, are new |
| merge | several clusters pick the same group; one continues, the others end, even if they left a remnant elsewhere |
| dissolution | a cluster that contributed at most one molecule to every group picks nothing and ends |

So a mixed dimer (one molecule from each of two origins) is always new, and each id
continues in at most one group.

## Known behaviours

These are deterministic but arbitrary choices, worth knowing before changing
anything:

- **Even splits** keep the id in the fragment that comes first in
  `nx.connected_components` order (resid order). Tests:
  `test_even_split_keeps_the_id_in_the_half_with_the_lowest_resid`,
  `test_even_split_into_equal_pieces_keeps_the_id_with_the_lowest_resid`.
- **No hysteresis**: a cluster that splits and re-merges one frame later mints an id
  that lives for one frame (the original keeps its id on the re-merge).
- Ids depend on the frames analysed: skipping frames (`--in-memory-step`) changes
  which events are seen and so the ids.
- Results from different engines (e.g. the same system in GROMACS and LAMMPS) can
  differ frame by frame where a distance sits right at a cutoff; through id
  continuity one such difference can relabel later clusters. Aggregate statistics
  should match.

## Changing it

Any change to which ids are assigned changes computed output. Whether that is a
breaking change depends on whether the previous result was valid (see
[Commits and releases](releases.md#breaking-changes)). The tracker tests
(`tests/test_tracker.py`) build frames from explicit groups of resids (the `track`
fixture, on top of `conftest.py`'s `make_universe`), so a new rule or tie-break gets a test that states the
layout and the expected ids directly.
