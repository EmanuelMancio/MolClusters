import networkx as nx
import MDAnalysis as mda
import numpy as np

from typing import Optional, Type, Iterable, Iterator, List
from MDAnalysis import core

from .conntable import ConnTable

EA2D = 1/0.3934303

class Cluster:
    __cls_id = 1

    __slots__ = [
        "uni",
        "initial_time",
        "cluster",
        "_cm",
        "_id",
        "_ag",
        "__centered_time",
    ]

    def __init__(
        self,
        universe: Type[mda.Universe],
        subconntab: Optional[Type["ConnTable.SubConnTable"]] = None,
    ) -> None:
        self.uni = universe
        self.initial_time: float = universe.trajectory.time
        self.cluster = None

        if subconntab:
            self.cluster: Type[nx.Graph] = nx.Graph(subconntab.graph)
            self._cm = (
                subconntab._cm
            )  # TODO: Change to keep the sum of center of masses of molecules
        else:
            self.cluster = nx.Graph()
            self._cm = np.empty(3)

        self._ag = core.groups.ResidueGroup(np.array(self.cluster) - 1, self.uni)

        self._id = Cluster.__cls_id
        Cluster.__cls_id += 1
        self.__centered_time = -np.inf

    @classmethod
    def _from_graph(cls, uni: Type[mda.Universe], graph: Type[nx.Graph]) -> "Cluster":
        tmp_cls = cls(uni, None)
        tmp_cls.cluster = graph
        tmp_cls.__recalculate_cm()
        return tmp_cls

    @property
    def ag(self):
        return self._ag

    @property
    def id(self) -> int:
        return self._id

    @property
    def cm(self) -> Type[np.ndarray]:
        return self._cm

    @cm.setter
    def cm(self, value: Type[np.ndarray]):
        self._cm = value

    def __recalculate_cm(self):
        self._cm = self._ag.center_of_mass()

    def __update_ag(self) -> None:
        self._ag = core.groups.ResidueGroup(np.array(self.cluster) - 1, self.uni)

    def add_mol(
        self,
        ref_mol: int,
        mol: int,
        resname: str,
        dist: float,
    ) -> None:
        if ref_mol not in self:
            raise ValueError(f"ref_mol {ref_mol} not in the cluster")

        if mol in self:
            raise ValueError(f"mol {mol} already in the cluster, use add_con instead")

        self.cluster.add_node(mol, name=resname)
        self.add_con(ref_mol, mol, dist)
        self._ag += core.groups.ResidueGroup([mol - 1], self.uni)
        self.__recalculate_cm()

    def add_con(
        self,
        moli: int,
        molj: int,
        dist: float,
    ) -> None:
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster, use add_mol instead")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster, use add_mol instead")

        if dist == 0.0:
            raise ValueError("dist is zero. Check your trajectory")

        self.cluster.add_edge(moli, molj, distance=dist, weight=np.exp(1 / dist))

    def remove_mol(self, mol: int) -> None:
        if mol not in self:
            raise ValueError(f"mol {mol} not in the cluster")

        self.cluster.remove_node(mol)
        self._ag -= core.groups.ResidueGroup([mol - 1], self.uni)
        self.__recalculate_cm()

    def remove_con(self, moli: int, molj: int) -> None:
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster")

        self.cluster.remove_edge(moli, molj)

    def remove_cons(self, ref_mol: int, cons: Iterable) -> None:
        for con in cons:
            self.remove_con(ref_mol, con)

    def get_age(self, time: float) -> float:
        return time - self.initial_time

    def merge(self, other: "Cluster") -> None:
        self.cluster = nx.compose(self.cluster, other._graph)
        self.__recalculate_cm()

    def separate(self) -> List[Type[nx.Graph]]:
        sub_clusters: List[Type[nx.Graph]] = [
            self.cluster.subgraph(c).copy()
            for c in sorted(
                nx.connected_components(self.cluster), key=len, reverse=True
            )
        ]

        self.cluster = sub_clusters[0]
        self.__update_ag()
        self.__recalculate_cm()

        return [Cluster._from_graph(self.uni, nx.Graph(g)) for g in sub_clusters[1:]]

    @property
    def _graph(self) -> Type[nx.Graph]:
        return self.cluster

    @property
    def size(self) -> int:
        return len(self)

    def update_from_conntable(self, conn: Type["ConnTable.SubConnTable"]):
        # removed_mols = set(self).difference(conn)
        # removed_cons = self.cluster.edges - conn.graph.edges
        # new_cons = conn.graph.edges - self.cluster.edges

        self.cluster = nx.Graph(conn.graph)
        self._cm = conn.cm
        self.__update_ag()

    def get_dist(self, moli: int, molj: int) -> float:
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster")

        return self.cluster[moli][molj]["distance"]

    def set_dist(self, moli, molj, dist):
        if dist == 0.0:
            raise ValueError("dist is zero")

        self.cluster[moli][molj]["weight"] = np.exp(1 / dist)
        self.cluster[moli][molj]["distance"] = dist

    def get_weight(self, moli: int, molj: int) -> float:
        if moli not in self:
            raise ValueError(f"mol {moli} not in the cluster")
        if molj not in self:
            raise ValueError(f"mol {molj} not in the cluster")

        return self.cluster[moli][molj]["weight"]

    def __make_whole(self):
        #! This should NOT be used before printing
        if self.__centered_time != self.uni.trajectory.time:
            boxcenter = np.sum(self.uni.trajectory.ts.triclinic_dimensions, axis=0) / 2
            # self._ag.atoms.unwrap(compound="residues",reference="cog",inplace=True)

            ref_mol_cm = self._ag[:1].center_of_mass(unwrap=True)
            vector = boxcenter - ref_mol_cm
            self._ag.atoms.positions += vector
            self._ag.atoms.unwrap(compound="residues", reference="cog", inplace=True)
            # center_in_box(self._ag,point=ref_mol_cm)(self.uni.trajectory.ts)
            # center_in_box(self._ag)(self.uni.trajectory.ts)
            self.__centered_time = self.uni.trajectory.time

    @property
    def resnames(self):
        return self._ag.resnames

    @property
    def resids(self):
        return self._ag.resids

    @property
    def mass(self):
        return self._ag.total_mass()

    @property
    def sphericity(self):
        self.__make_whole()
        return 1 - self._ag.asphericity()

    @property
    def dipole_moment(self):
        self.__make_whole()
        return self._ag.atoms.dipole_moment() * EA2D

    @property
    def dipole(self):
        self.__make_whole()
        return self._ag.atoms.dipole_vector() * EA2D

    @property
    def shape_parameter(self):
        self.__make_whole()
        return self._ag.shape_parameter()

    @property
    def bsphere(self):
        self.__make_whole()
        return self._ag.bsphere()

    @property
    def radius_of_gyration(self):
        self.__make_whole()
        return self._ag.radius_of_gyration()

    @property
    def volume(self):
        r = self.radius_of_gyration
        return 4 * np.pi * r**3 / 3  # angstrom^3

    @property
    def density(self):
        return (self.mass / self.volume) * 0.602214076  # g/cm^3

    @property
    def charge(self):
        return self._ag.total_charge()

    def __contains__(self, item: int) -> bool:
        return item in self.cluster

    def __iter__(self) -> Iterator:
        return iter(self.cluster.copy())

    def __getitem__(self, key: int):
        return self.cluster[key]

    def __eq__(self, other) -> bool:
        if not isinstance(other, Cluster):
            return False

        return nx.utils.graphs_equal(self.cluster, other.cluster)

    def __len__(self) -> int:
        return len(self.cluster)

    def __str__(self):
        return self.cluster.edges.data().__str__()
