"""
besmarts.core.graph_bijections
"""

import functools

from besmarts.core import topology
from besmarts.core import graphs
from besmarts.core import mapper
from besmarts.core import graph_bitwise


# TODO convert to workspaces?
class map_vertices_ctx:
    cg = None
    o = None
    a = None
    b = None
    H = None
    strict = False
    equality = False


class bijection:
    """
    A pair of graphs with an bijective node mapping between them
    """

    __slots__ = "topology", "g", "h", "t"

    def __init__(
        self,
        topo: topology.structure_topology,
        g: graphs.subgraph,
        h: graphs.subgraph,
        t: dict,
    ):
        self.topology = topo
        self.g = g
        self.h = h
        self.t = t


def mapped_to_bijection(M: mapped_type) -> bijection:
    g = graphs.structure_to_subgraph(M.G)
    g = graphs.structure_to_subgraph(M.H)
    t = M.map
    F = bijection(G.topology, g, h, t)
    return F

def bijection_mcs(G: graphs.structure, H: graphs.structure) -> bijection:
    M = mapper.map_to(G, H, add_nodes=0)
    F = mapped_to_bijection(M)
    return F

def bijection_constant_left(
    G: graphs.structure, H: graphs.structure
) -> bijection:
    M = mapper.map_to(G, H, add_nodes=2, fill=True)
    F = mapped_to_bijection(M)
    return F

def bijection_constant_right(
    G: graphs.structure, H: graphs.structure
) -> bijection:
    M = mapper.map_to(G, H, add_nodes=2, fill=True)
    F = mapped_to_bijection(M)
    return F


def bijection_branch(G: graphs.structure, H: graphs.structure) -> bijection:
    M = mapper.map_to(G, H, add_nodes=3, fill=True)
    F = mapped_to_bijection(M)
    return F


def bijection_invert(F: bijection) -> bijection:
    """
    Transform the map from G -> H to G <- H

    Parameters
    ----------
    T : bijection
        The input map

    Returns
    -------
        A new mapped type
    """
    m = {v: k for k, v in F.t.items()}

    return bijection(F.topology, F.h, F.g, m)


def bijection_compose_pair(F1: bijection, F2: bijection) -> bijection:
    """
    Transform the map from G -> H to G <- H

    Parameters
    ----------
    T : bijection
        The input map

    Returns
    -------
        A new mapped type
    """
    # assert T1.H == T2.G
    assert F1.topology == F2.topology
    t = {i: F1.t[j] for i, j in F2.t.items()}
    return bijection(F1.topology, F1.g, F2.h, t)


def bijection_compose(F: list[bijection]) -> bijection:
    return functools.reduce(bijection_compose_pair, F)


def bijection_union(F: bijection) -> graphs.structure:
    return graph_bitwise.union(F.g, F.h, F.t)


def bijection_intersection(F: bijection) -> graphs.structure:
    return graph_bitwise.intersection(F.g, F.h, F.t)


def bijection_difference(F: bijection) -> graphs.structure:
    return graph_bitwise.difference(F.g, F.h, F.t)


def bijection_symmetric_difference(F: bijection) -> graphs.structure:
    return graph_bitwise.symmetric_difference(F.g, F.h, F.t)


def bijection_subset(F: bijection) -> graphs.structure:
    return graph_bitwise.symmetric_difference(F.g, F.h, F.t)


def bijection_equal(F: bijection) -> graphs.structure:
    return graph_bitwise.symmetric_difference(F.g, F.h, F.t)


def bijection_superset(F: bijection) -> graphs.structure:
    return graph_bitwise.symmetric_difference(F.g, F.h, F.t)


def bijection_subset_equal(F: bijection) -> graphs.structure:
    return graph_bitwise.symmetric_difference(F.g, F.h, F.t)


def bijection_superset_equal(F: bijection) -> graphs.structure:
    return graph_bitwise.symmetric_difference(F.g, F.h, F.t)


def structure_union_bijection_mcs(
    G: graphs.structure, H: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(G, H, add_nodes=0)
    F = mapped_to_bijection(M)
    g = graph_bitwise.union(F)
    return g


def structure_union_bijection_constant_left(
    G: graphs.structure, H: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(G, H, add_nodes=2)
    F = mapped_to_bijection(M)
    g = graph_bitwise.union(F)
    return g


def structure_union_bijection_constant_right(
    G: graphs.structure, H: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(G, H, add_nodes=1)
    F = mapped_to_bijection(M)
    g = graph_bitwise.union(F)
    return g


def structure_union_bijection_branch(
    G: graphs.structure, H: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(G, H, add_nodes=3)
    F = mapped_to_bijection(M)
    g = bijection_union(F)
    return g


def structure_intersection_bijection_mcs(
    G: graphs.structure, H: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(G, H, add_nodes=0, fill=True)
    F = mapped_to_bijection(M)
    g = bijection_intersection(F)
    return g


def structure_intersection_bijection_constant_left(
    G: graphs.structure, H: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(G, H, add_nodes=2, fill=True)
    F = mapped_to_bijection(M)
    g = bijection_intersection(F)
    return g


def structure_intersection_bijection_constant_right(
    G: graphs.structure, H: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(G, H, add_nodes=1, fill=True)
    F = mapped_to_bijection(M)
    g = bijection_intersection(F)
    return g


def structure_intersection_bijection_branch(
    G: graphs.structure, H: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(G, H, add_nodes=2, fill=True)
    F = mapped_to_bijection(M)
    g = bijection_intersection(F)
    return g


def structure_symmetric_difference_bijection_mcs(
    G: graphs.structure, H: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(G, H, add_nodes=0, fill=False)
    F = mapped_to_bijection(M)
    g = bijection_symmetric_difference(F)
    return g


def structure_symmetric_difference_bijection_constant_left(
    G: graphs.structure, H: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(G, H, add_nodes=2, fill=False)
    F = mapped_to_bijection(M)
    g = bijection_symmetric_difference(F)
    return g


def structure_symmetric_difference_bijection_constant_right(
    G: graphs.structure, H: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(G, H, add_nodes=1, fill=False)
    F = mapped_to_bijection(M)
    g = bijection_symmetric_difference(F)
    return g


def structure_symmetric_difference_bijection_branch(
    G: graphs.structure, H: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(G, H, add_nodes=2, fill=False)
    F = mapped_to_bijection(M)
    g = bijection_symmetric_difference(F)
    return g


def structure_difference_bijection_mcs(
    G: graphs.structure, H: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(G, H, add_nodes=0, fill=False)
    F = mapped_to_bijection(M)
    g = bijection_symmetric_difference(F)
    return g


def structure_difference_bijection_constant_left(
    G: graphs.structure, H: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(G, H, add_nodes=2, fill=False)
    F = mapped_to_bijection(M)
    g = bijection_symmetric_difference(F)
    return g


def structure_difference_bijection_constant_right(
    G: graphs.structure, H: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(G, H, add_nodes=1, fill=True)
    F = mapped_to_bijection(M)
    g = bijection_symmetric_difference(F)
    return g


def structure_difference_bijection_branch(
    G: graphs.structure, H: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(G, H, add_nodes=3, fill=True)
    F = mapped_to_bijection(M)
    g = bijection_symmetric_difference(F)
    return g
