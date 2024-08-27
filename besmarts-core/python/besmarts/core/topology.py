"""
besmarts.core.topology

Topology in BESMARTS describe the structure of a graph, and are a necessary
component of structure graphs.
"""

from typing import Sequence, Tuple
from itertools import permutations


class structure_topology:
    """
    The structure topology describes several properties of a graph that allow
    it to be a structure. The data in a structure is agnostic to the graph and
    uses relative indexing of the nodes. All graphs that have the same
    structure, such as angles, use the same topology.
    """

    __slots__ = ("primary", "connect", "permutations")

    def __init__(
        self,
        primary: Sequence[int],
        connect: Sequence[Tuple[int, int]],
        permutations: Sequence[Sequence[int]],
    ):
        self.primary = primary
        self.connect = connect
        self.permutations = permutations

    def __eq__(self, o):
        return (
            self.primary == o.primary
            and self.connect == o.connect
            and self.permutations == o.permutations
        )

    def __neq__(self, o):
        return not self.__eq__(o)

    def __hash__(self):
        return hash((self.primary, self.connect, self.permutations))

class transcode_topology:
    """
    A pair of topologies that defines a transformation from A to B
    """

    __slots__ = ("A", "B", "transcode")

    def __init__(
        self,
        A: structure_topology,
        B: structure_topology,
        transcode,
    ):
        self.A: structure_topology= A
        self.B: structure_topology= B
        self.transcode: Dict[Sequence[Sequence[int]], Sequence[Sequence[int]]] = transcode

    def __eq__(self, o):
        return (
            self.A == o.A
            and self.B == o.B
            and self.transcode == o.transcode
        )

    def __neq__(self, o):
        return not self.__eq__(o)

    def __hash__(self):
        return hash((self.A, self.B, self.transcode))

def null_topology() -> structure_topology:
    """
    Return a topology that describes a null/undefined topology. This is used
    when the data/structures describe don't have a well defined graph topology,
    such as entire molecules where the permutations aren't known.

    Parameters
    ----------

    Returns
    -------
    structure_topology
    """
    return structure_topology((), (), ((),))


def atom_topology() -> structure_topology:
    """
    Return a topology that describes an atom

    Parameters
    ----------

    Returns
    -------
    structure_topology
    """
    return structure_topology((0,), (), ((0,),))


def bond_topology() -> structure_topology:
    """
    Return a topology that describes a bond

    Parameters
    ----------

    Returns
    -------
    structure_topology
    """
    return structure_topology((0, 1), ((0, 1),), ((0, 1), (1, 0)))


def angle_topology() -> structure_topology:
    """
    Return a topology that describes an angle

    Parameters
    ----------

    Returns
    -------
    structure_topology
    """
    return structure_topology(
        (0, 1, 2), ((0, 1), (1, 2)), ((0, 1, 2), (2, 1, 0))
    )


def torsion_topology() -> structure_topology:
    """
    Return a topology that describes a torsion dihedral

    Parameters
    ----------

    Returns
    -------
    structure_topology
    """
    return structure_topology(
        (0, 1, 2, 3), ((0, 1), (1, 2), (2, 3)), ((0, 1, 2, 3), (3, 2, 1, 0))
    )


def outofplane_topology() -> structure_topology:
    """
    Return a topology that describes an out-of-plane dihedral

    Parameters
    ----------

    Returns
    -------
    structure_topology
    """
    return structure_topology(
        (0, 1, 2, 3),
        ((0, 1), (1, 2), (1, 3)),
        (
            (0, 1, 2, 3),
            (0, 1, 3, 2),
            (2, 1, 0, 3),
            (2, 1, 3, 0),
            (3, 1, 0, 2),
            (3, 1, 2, 0),
        ),
    )


def chain_topology(n: int) -> structure_topology:
    """
    Return a topology that describes an arbitrary linear sequence of atoms.

    Parameters
    ----------
    n : int
        The length of the chain

    Returns
    -------
    structure_topology
    """
    return structure_topology(
        tuple(range(n)),
        tuple(((i, i + 1) for i in range(n))),
        (tuple(range(n)), tuple(reversed(range(n)))),
    )


def ring_topology(n: int) -> structure_topology:
    """
    Return a topology that describes a ring of atoms of arbitrary length.

    Parameters
    ----------
    n : int
        The length of the ring

    Returns
    -------
    structure_topology
    """
    cycles = tuple(tuple(i % n for i in range(j, n + j)) for j in range(n - 1))
    return structure_topology(
        tuple(range(n)),
        tuple(((i, (i + 1) % n) for i in range(n))),
        (*cycles, *(tuple(reversed(x)) for x in cycles)),
    )


def n_body_topology(n: int) -> structure_topology:
    """
    Return a topology that describes n nonbonded atoms

    Parameters
    ----------
    n : int
        The number of bodies

    Returns
    -------
    structure_topology
    """

    return structure_topology(
        tuple(range(n)),
        (),
        tuple((tuple(x) for x in permutations(range(n), n))),
    )


def pair_topology() -> structure_topology:
    """
    Return a topology that describes a pair, or a set of atoms that may not be
    bonded

    Parameters
    ----------

    Returns
    -------
    structure_topology
    """
    return n_body_topology(2)


# Singletons
null = null_topology()
atom = atom_topology()
bond = bond_topology()
angle = angle_topology()
torsion = torsion_topology()
outofplane = outofplane_topology()
pair = n_body_topology(2)
triplet = n_body_topology(3)

atom_to_atom = transcode_topology(
    atom,
    atom,
    {((0,),) : ((0,),)}
)

# 2x1
bond_to_atom = transcode_topology(
    bond,
    atom,
    [[1],[1]]
)

# 3x1 mapping
# each have an equal weight
angle_to_atom = transcode_topology(
    angle,
    atom,
    [[1],[1],[1]]
)

torsion_to_atom = transcode_topology(
    torsion,
    atom,
    {((0,), (1,), (2,), (3,)) : ((0,1,2,3),)}
)

torsion_to_bond = transcode_topology(
    torsion,
    bond,
    {0:0, 1:0, 2:1, 3:1}
)


topology_index = [
    null,
    atom,
    bond,
    angle,
    torsion,
    outofplane,
    pair,
    triplet,
]

def index_of(topo: structure_topology):
    return topology_index.index(topo)

