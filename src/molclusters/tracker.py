# SPDX-FileCopyrightText: © 2024 Emanuel Mancio <emanuelmancio@usp.br>
#
# SPDX-License-Identifier: GPL-3.0-only

"""Provides the `ClusterTracker` class, which follows molecular clusters over a trajectory.

The tracker only keeps the clusters and their ids up to date, frame by frame; the
analyses built on them (and all file output) live in `molclusters.molclusters`.

Classes:
--------
- ClusterTracker: Finds the clusters of the current frame and keeps their ids stable
  across formation, growth, split, merge and dissolution events.
"""

from collections import Counter, defaultdict
from itertools import count

import MDAnalysis as mda
from MDAnalysis import core

from .cluster import Cluster
from .config import MolClsConfig
from .conntable import ConnectionTable


class ClusterTracker:
    """Follows the molecular clusters of a trajectory, keeping their ids stable.

    The clusters of the frame the Universe is on are found on creation; after
    moving the trajectory to another frame, `update` reconciles them with it.

    Attributes
    ----------
    uni : mda.Universe
        The MDAnalysis Universe object associated with the simulation.
    config : MolClsConfig
        The analysis configuration (only its rules, distance backend and ignored
        compositions are used).
    sels : dict[str, core.groups.AtomGroup]
        Atom groups for each residue type.
    conntab : ConnectionTable
        The connectivity table for molecular clusters.
    clusters : dict[int, Cluster]
        A dictionary of detected clusters, keyed by cluster ID.
    mol_clt : dict[int, int]
        A mapping of molecule IDs to their respective cluster IDs.
    """

    __slots__ = ["uni", "config", "sels", "conntab", "clusters", "mol_clt", "_ids"]

    def __init__(self, universe: mda.Universe, config: MolClsConfig) -> None:
        """Initialize the tracker with the clusters of the current frame.

        Parameters
        ----------
        universe : mda.Universe
            The MDAnalysis Universe object associated with the simulation.
        config : MolClsConfig
            The analysis configuration.
        """
        self.uni = universe
        self.config = config

        self.sels: dict[str, core.groups.AtomGroup] = {
            res: self.uni.select_atoms(f"resname {res}")
            for res in config._rules.all_keys()
        }

        self.conntab = ConnectionTable(
            self.uni,
            self.config._rules,
            self.sels,
            backend=self.config.distance_backend,
        )
        self.clusters: dict[int, Cluster] = {}
        self.mol_clt: dict[int, int] = {}
        # each tracker numbers its own clusters, from 1
        self._ids = count(1)

        self._start_clusters()

    def _start_clusters(self) -> None:
        """Initialize clusters at the beginning of the analysis."""
        for subconn in self.conntab.subconntables():
            if self.config.is_ignored_composition(subconn.resnames):
                continue

            cls_id = self._create_new_cluster(subconn)

            for mol in subconn:
                self.mol_clt[mol] = cls_id

    # TODO: make a better name for this function
    def _gen_origin_cluster_counter(
        self, subconn: ConnectionTable._SubConnTable
    ) -> Counter:
        """Count where the molecules of a new connected group came from.

        Each molecule of `subconn` (a connected group of the current frame) is
        looked up in `mol_clt`, which still holds the *previous* frame's assignment,
        so the result says how many of the connected group's molecules each
        previous cluster contributed. Molecules that were free in the previous frame
        are counted under id ``0`` (safe as a sentinel because cluster ids start
        at 1).

        For example, a connected group made of 4 molecules of cluster 7, 2 of
        cluster 9 and 1 free molecule gives ``Counter({7: 4, 9: 2, 0: 1})``.

        Parameters
        ----------
        subconn : ConnectionTable._SubConnTable
            A connected group of the current frame.

        Returns
        -------
        Counter
            Previous-frame cluster id (``0`` for free molecules) -> number of the
            connected group's molecules that came from it.
        """
        mols_origin_clusters = {mol: self.mol_clt.get(mol, 0) for mol in subconn}
        return Counter(mols_origin_clusters.values())

    @staticmethod
    def _best_groups(
        conn_info: list[tuple[ConnectionTable._SubConnTable, Counter]],
    ) -> dict[int, int]:
        """Find the connected group that best continues each previous cluster.

        A previous cluster's best group is the one holding the most of its
        molecules. On a tie, the smallest group wins, i.e. the purest one (the
        largest share of its molecules came from the cluster), so a fragment made
        only of the cluster's own molecules beats one mixed with newcomers. If the
        groups are also the same size, the first one in `conn_info` wins.

        A contribution of a single molecule is ignored: one molecule joining other
        molecules doesn't carry the cluster's identity. A cluster with no
        contribution of two or more molecules has no best group and dies.

        Every group is looked at before any choice is made, so the result doesn't
        depend on the order the groups are processed in.

        Parameters
        ----------
        conn_info : list[tuple[ConnectionTable._SubConnTable, Counter]]
            Every connected group of the current frame, largest first, each with its
            origin counter (see `_gen_origin_cluster_counter`).

        Returns
        -------
        dict[int, int]
            Previous-frame cluster id -> index in `conn_info` of its best group.
        """
        best: dict[int, tuple[tuple[int, int], int]] = {}

        for i, (subconn, origin_clusters) in enumerate(conn_info):
            for cls_id, n in origin_clusters.items():
                if not cls_id or n == 1:  # free molecules, or a single molecule
                    continue

                key = (n, -len(subconn))
                if cls_id not in best or key > best[cls_id][0]:
                    best[cls_id] = (key, i)

        return {cls_id: i for cls_id, (_, i) in best.items()}

    def _create_new_cluster(self, subconn: ConnectionTable._SubConnTable) -> int:
        """Create a new cluster from a subconnection table.

        Parameters
        ----------
        subconn : ConnectionTable._SubConnTable
            The subconnection table from which to create the new cluster.

        Returns
        -------
        int
            The ID of the newly created cluster.
        """
        cluster = Cluster(self.uni, subconn, cluster_id=next(self._ids))
        id = cluster.id
        self.clusters[id] = cluster

        return id

    @staticmethod
    def _get_older_cluster(candidates: list[int], origin_clusters: Counter) -> int:
        """Pick which candidate cluster keeps its id in a merge.

        The candidate that contributed the most molecules wins, so a large cluster
        absorbing a small one keeps its id even if it is younger. If several
        candidates tie for the most molecules, the oldest wins. Cluster ids are
        handed out in increasing order, so the oldest is the lowest id. The order of
        `candidates` doesn't matter.

        Parameters
        ----------
        candidates : list[int]
            Non-empty list of the previous clusters whose best group (see
            `_best_groups`) is this connected group.
        origin_clusters : Counter
            Origin counter of the connected group.

        Returns
        -------
        int
            Id of the cluster that continues as the connected group.
        """
        return max(candidates, key=lambda c: (origin_clusters[c], -c))

    def update(self) -> None:
        """Reconcile the current frame's connected groups with the previous clusters.

        This is the dominance algorithm, which keeps cluster ids stable while
        clusters grow, shrink, split and merge.

        A *connected group* (``subconn`` in the code, a
        `ConnectionTable._SubConnTable`) is a set of molecules linked to each
        other, directly or through other molecules, by the config's rules in the
        current frame, and to nothing outside the set. Free (unlinked) molecules
        form no connected group. A connected group has no id yet; a *cluster* is
        what it becomes once this method gives it one, new or inherited from the
        previous frame.

        The connection table is rebuilt for the current frame, and the previous
        clusters that the molecules of each connected group (ignored compositions
        excluded) came from are counted (see `_gen_origin_cluster_counter`). Ids
        are then assigned in two steps, each of which sees the whole frame:

        1. Each previous cluster picks its best group: the one holding the most of
           its molecules, the purest on a tie (see `_best_groups`).
        2. Each connected group chosen by one or more previous clusters continues
           the one that contributed the most molecules, the oldest on a tie (see
           `_get_older_cluster`). A group chosen by nobody is a new cluster.

        The usual events follow from these two rules:

        - formation: a group of free molecules is chosen by nobody -> new;
        - growth or shrinking: a group chosen by one cluster continues it;
        - split: the cluster continues in its best fragment and the other fragments,
          chosen by nobody, are new;
        - merge: several clusters choose the same group, one continues and the
          others die, even if they also left a remnant elsewhere;
        - dissolution: a cluster that contributed at most one molecule to every
          group has no best group and dies.

        So a mixed dimer (one molecule from each origin) is always new, and each id
        continues in at most one group. Previous clusters that no connected group
        continues are dropped, along with the `mol_clt` entries of molecules that
        are now free.
        """
        self.conntab.update()

        modified_mols = set()
        modified_clusters = set()

        conn_info = [
            (sub, self._gen_origin_cluster_counter(sub))
            for sub in self.conntab.subconntables()
            if not self.config.is_ignored_composition(sub.resnames)
        ]

        candidates = defaultdict(list)
        for cls_id, i in self._best_groups(conn_info).items():
            candidates[i].append(cls_id)

        for i, (subconn, origin_clusters) in enumerate(conn_info):
            if candidates[i]:
                id = self._get_older_cluster(candidates[i], origin_clusters)
                self.clusters[id]._update(subconn)
            else:
                id = self._create_new_cluster(subconn)

            for mol in subconn:
                self.mol_clt[mol] = id

            modified_clusters.add(id)
            modified_mols.update(subconn)

        for mol in set(self.mol_clt.keys()).difference(modified_mols):
            self.mol_clt.pop(mol)

        # TODO: deal with clusters that weren't modified. Needs to consider that some clusters merged (for log filing)  # noqa: E501
        for cls in set(self.clusters.keys()).difference(modified_clusters):
            self.clusters.pop(cls)

    def find(self, mol: int) -> int | bool:
        """Find the cluster ID for a given molecule.

        Parameters
        ----------
        mol : int
            The molecule ID to search for.

        Returns
        -------
        int | bool
            The cluster ID if the molecule is found, or False if not found.
        """
        return self.mol_clt.get(mol, False)
