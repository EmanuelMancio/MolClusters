import MDAnalysis as mda
import networkx as nx
import numpy as np
from MDAnalysis import core
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis


import warnings
from typing import Dict, List, Tuple, Type, Union


class ConnTable:
    __slots__ = ["uni", "clst_args", "sels", "conntab", "cms", "hbs"]

    class SubConnTable:
        __slots__ = ["conntab", "_graph", "_cm"]

        def __init__(
            self, subgraph: Type[nx.Graph], conntable: Type["ConnTable"]
        ) -> None:
            self.conntab = conntable
            self._graph = subgraph

            self._cm = self.__selection().center_of_mass()

        def __selection(self) -> Type[core.groups.ResidueGroup]:
            return core.groups.ResidueGroup(
                [m - 1 for m in self.graph], self.conntab.uni
            )

        def __getitem__(self, key):
            return self.conntab[key]

        def __len__(self):
            return len(self._graph)

        @property
        def cm(self) -> Type[np.ndarray]:
            return self._cm

        @property
        def graph(self):
            return self._graph

        def __iter__(self):
            return iter(self._graph)

    def __init__(
        self,
        universe: Type[mda.Universe],
        cluster_args: Dict[str, Dict[str, Tuple[str, float]]],
        selections: Dict[str, Type[core.groups.AtomGroup]],
    ) -> None:
        self.uni = universe
        self.clst_args = cluster_args
        self.sels = selections

        self.cms = {}
        self.hbs: Dict[str, Dict[str, Type[HydrogenBondAnalysis]]] = {}
        self.__get_hbonds()
        self.update()

    def __get_mass_centers(self):
        for res in self.clst_args:
            self.cms[res] = self.sels[res].center_of_mass(compound="residues")

    def __get_hbonds(self):
        for resi in self.clst_args:
            for resj in self.clst_args[resi]:
                if self.clst_args[resi][resj] == "cm":
                    continue

                if (
                    (resi not in self.hbs)
                    or (resj not in self.hbs[resi])
                    or (resj not in self.hbs)
                    or (resi not in self.hbs[resj])
                ):
                    hb = HydrogenBondAnalysis(
                        self.uni,
                        between=[f"resname {resi}", f"resname {resj}"],
                        d_a_cutoff=self.clst_args[resi][resj][1]["d"],
                        d_h_a_angle_cutoff=self.clst_args[resi][resj][1]["a"],
                        update_selections=False,
                    )

                    hb._prepare()

                    if resi not in self.hbs:
                        self.hbs[resi] = {resj: [hb, 0]}
                    else:
                        self.hbs[resi][resj] = [hb, 0]

                    if resj not in self.hbs:
                        self.hbs[resj] = {resi: [hb, 0]}
                    else:
                        self.hbs[resj][resi] = [hb, 0]

    def __get_pair_and_distances(
        self,
        resi: str,
        resj: str,
        box: Type[np.ndarray],
    ) -> Tuple[Type[np.ndarray], Type[np.ndarray]]:
        pairs: Type[np.ndarray]
        distances: Type[np.ndarray]

        if self.clst_args[resi][resj][0] == "cm":
            cm1: Type[np.ndarray] = self.cms[resi]
            cutoff = self.clst_args[resi][resj][1]
            if resi != resj:
                cm2: Type[np.ndarray] = self.cms[resj]
                pairs, distances = mda.lib.distances.capped_distance(
                    cm1, cm2, cutoff, box=box
                )
            else:
                pairs, distances = mda.lib.distances.self_capped_distance(
                    cm1, cutoff, box=box
                )

            for k, (moli, molj) in enumerate(pairs):
                pairs[k, 0] = self.sels[resi].residues[moli].resid
                pairs[k, 1] = self.sels[resj].residues[molj].resid
        else:
            hb = self.hbs[resi][resj][0]
            hb._ts = self.uni.trajectory.ts

            # suppress warnings when there are no HBonds
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                hb._single_frame()

            frame_id = self.hbs[resi][resj][1]
            res = (np.asarray(hb.results.hbonds).T)[frame_id:, -4:-1]
            distances = res[:, -2]
            pairs = np.empty((0, 2), int)
            for h_ati, a_ati, _ in res:
                h_ati, a_ati = int(h_ati), int(a_ati)
                moli = self.uni.atoms[h_ati].resid
                molj = self.uni.atoms[a_ati].resid
                pairs = np.append(pairs, [[moli, molj]], axis=0)

            self.hbs[resi][resj][1] += len(res)

        return pairs, distances

    def __construct_table(self) -> Type[nx.Graph]:
        self.conntab = nx.Graph()

        analyzed = dict.fromkeys(self.clst_args.keys())
        for k in analyzed:
            analyzed[k] = []

        for resi in self.clst_args:
            for resj in self.clst_args[resi]:
                if resj not in analyzed[resi]:
                    pairs: Type[np.ndarray]
                    distances: Type[np.ndarray]

                    pairs, distances = self.__get_pair_and_distances(
                        resi, resj, self.uni.dimensions
                    )

                    for k, (ri, rj) in enumerate(pairs):
                        self.conntab.add_edge(ri, rj, d=distances[k])

                    analyzed[resi].append(resj)
                    analyzed[resj].append(resi)

    def update(self) -> None:
        self.__get_mass_centers()
        self.__construct_table()

    def __getitem__(self, key: Union[Tuple[int, int], int]) -> Union[float, List[int]]:
        if isinstance(key, tuple):
            if len(key) > 2:
                raise KeyError(
                    "ConnTable accepts only one or two parameters to access data!"
                )

            if key[0] not in self:
                raise IndexError(f"{key[0]} not in ConnTable")

            if key[1] not in self:
                raise IndexError(f"{key[1]} not in ConnTable")

            return self.conntab[key[0]][key[1]]["d"]

        if key not in self:
            raise IndexError(f"{key} not in ConnTable")

        return self.mols_connected_to(key)

    def __contains__(self, item: int) -> bool:
        return item in self.conntab

    def __graph_from_mol(self, mol):
        desc = nx.node_connected_component(self.conntab, mol)
        return self.conntab.subgraph(desc)

    def __iter__(self):
        return iter(self.conntab)

    def connections_from(self, mol):
        return list(self.conntab.edges(mol))

    def connection_tree_from(self, mol):
        g = self.__graph_from_mol(mol)
        return list(g.edges)

    def mols_connected_to(self, mol):
        return list(self.conntab[mol])

    def mols_connected_tree_to(self, mol):
        g = self.__graph_from_mol(mol)
        return list(g.nodes)

    def __subgraphs(self):
        for c in nx.connected_components(self.conntab):
            yield self.conntab.subgraph(c)

    def subconntables(self):
        for s in sorted(
            self.__subgraphs(), key=lambda x: len(x), reverse=True
        ):  # sorts to guarantee that in case of separation the biggest cluster keeps
            # the id
            yield self.SubConnTable(s, self)
