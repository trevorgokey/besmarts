"""
besmarts.bijections

Definitions and functions to define a mapping between two input graphs that 
are potentially modified such that a bijection is formed.
"""

import functools

from besmarts import graphs, topology
# from besmarts.core import mapper
from besmarts.core import graph

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

def structure_bijection(
    g: graphs.structure,
    h: graphs.structure,
    add_nodes=0,
    fill=0
):
    T = mappings.map_to(g, h, add_nodes=add_nodes, fill=fill)
    F = mapping_to_bijection(T)
    return F

def mapping_to_bijection(M: mappings.mapping) -> bijection:
    g = graphs.structure_to_subgraph(M.g)
    h = graphs.structure_to_subgraph(M.h)
    t = M.map
    assert not (
        set(t.keys()).symmetric_difference(g.nodes) or 
        set(t.values()).symmetric_difference(h.nodes)
    )
    F = bijection(M.g.topology, g, h, t)
    return F

def structure_bijection_mcs(g: graphs.structure, h: graphs.structure) -> bijection:
    F = structure_bijection(g, h, add_nodes=0)
    return F

def structure_bijection_constant_right(
    g: graphs.structure, h: graphs.structure
) -> bijection:
    F = structure_bijection(g, h, add_nodes=1, fill=True)
    return F

def structure_bijection_constant_left(
    g: graphs.structure, h: graphs.structure
) -> bijection:
    F = structure_bijection(g, h, add_nodes=2, fill=True)
    return F

def structure_bijection_branch(g: graphs.structure, h: graphs.structure) -> bijection:
    F = structure_bijection(g, h, add_nodes=3, fill=True)
    return F

def bijection_invert(F: bijection) -> bijection:
    """
    Transform the map from g -> h to g <- h

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
    Transform the map from g -> h to g <- h

    Parameters
    ----------
    T : bijection
        The input map

    Returns
    -------
        A new mapped type
    """
    # assert T1.h == T2.g
    assert F1.topology == F2.topology
    t = {i: F1.t[j] for i, j in F2.t.items()}
    return bijection(F1.topology, F1.g, F2.h, t)


def bijection_compose(F: list[bijection]) -> bijection:
    return functools.reduce(bijection_compose_pair, F)

def bijection_union(F: bijection) -> graphs.subgraph:
    return graph.bitwise.subgraph_union(F.g, F.h, F.t)

def bijection_intersection(F: bijection) -> graphs.structure:
    return graph.bitwise.subgraph_intersection(F.g, F.h, F.t)

def bijection_xor(F: bijection) -> graphs.structure:
    return graph.bitwise.subgraph_xor(F.g, F.h, F.t)

def bijection_subtract(F: bijection) -> graphs.structure:
    return graph.bitwise.subgraph_subtract(F.g, F.h, F.t)

def bijection_symmetric_difference(F: bijection) -> graphs.structure:
    return graph.bitwise.symmetric_difference(F.g, F.h, F.t)


def bijection_subset(F: bijection) -> bool:
    diff = graph.bitwise.subgraph_subtract(F.g, F.h, F.t)
    bits = graphs.graph_bits(diff)
    return bits == 0


def bijection_equal(F: bijection) -> bool:
    diff = graph.bitwise.subgraph_xor(F.g, F.h, F.t)
    bits = graphs.graph_bits(diff)
    return bits == 0


def bijection_superset(F: bijection) -> graphs.structure:
    F = bijection_invert(F)
    return bijection_subset(F)


def bijection_subset_equal(F: bijection) -> graphs.structure:
    return bijection_subset(F) or bijection_equal(F)


def bijection_superset_equal(F: bijection) -> graphs.structure:
    F = bijection_invert(F)
    return bijection_subset_equal(F)


def structure_union_bijection_mcs(
    g: graphs.structure, h: graphs.structure
) -> graphs.structure:
    M = mappings.map_to(g, h, add_nodes=0)
    F = mapping_to_bijection(M)
    g = graph_bitwise.union(F)
    return g


def structure_union_bijection_branch(
    g: graphs.structure, h: graphs.structure
) -> graphs.structure:
    M = mappings.map_to(g, h, add_nodes=1, fill=False)
    F = mapping_to_bijection(M)
    g = bijection_union(F)
    return g


def structure_union_bijection_constant_right(
    g: graphs.structure, h: graphs.structure
) -> graphs.structure:
    M = mappings.map_to(g, h, add_nodes=2, fill=False)
    F = mapping_to_bijection(M)
    g = graph_bitwise.union(F)
    return g


def structure_union_bijection_constant_left(
    g: graphs.structure, h: graphs.structure
) -> graphs.structure:
    M = mappings.map_to(g, h, add_nodes=3, fill=False)
    F = mapping_to_bijection(M)
    g = graph_bitwise.union(F)
    return g



###
# Intersections
###


def structure_intersection_bijection_mcs(
    g: graphs.structure, h: graphs.structure
) -> graphs.structure:
    M = mappings.map_to(g, h, add_nodes=0, fill=True)
    F = mapping_to_bijection(M)
    g = bijection_intersection(F)
    return g


def structure_intersection_bijection_constant_left(
    g: graphs.structure, h: graphs.structure
) -> graphs.structure:
    M = mappings.map_to(g, h, add_nodes=3, fill=True)
    F = mapping_to_bijection(M)
    g = bijection_intersection(F)
    return g


def structure_intersection_bijection_constant_right(
    g: graphs.structure, h: graphs.structure
) -> graphs.structure:
    M = mappings.map_to(g, h, add_nodes=2, fill=True)
    F = mapping_to_bijection(M)
    g = bijection_intersection(F)
    return g


def structure_intersection_bijection_branch(
    g: graphs.structure, h: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(g, h, add_nodes=1, fill=True)
    F = mapping_to_bijection(M)
    g = bijection_intersection(F)
    return g

###
# Symmetric Difference: fill is True as it assumes missings nodes are anything
###

def structure_symmetric_difference_bijection_mcs(
    g: graphs.structure, h: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(g, h, add_nodes=0, fill=True)
    F = mapping_to_bijection(M)
    g = bijection_symmetric_difference(F)
    return g


def structure_symmetric_difference_bijection_branch(
    g: graphs.structure, h: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(g, h, add_nodes=1, fill=True)
    F = mapping_to_bijection(M)
    g = bijection_symmetric_difference(F)
    return g


def structure_symmetric_difference_bijection_constant_right(
    g: graphs.structure, h: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(g, h, add_nodes=2, fill=True)
    F = mapping_to_bijection(M)
    g = bijection_symmetric_difference(F)
    return g


def structure_symmetric_difference_bijection_constant_left(
    g: graphs.structure, h: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(g, h, add_nodes=3, fill=True)
    F = mapping_to_bijection(M)
    g = bijection_symmetric_difference(F)
    return g


###
# Difference: fill is True as it assumes missings nodes are anything
###

def structure_difference_bijection_mcs(
    g: graphs.structure, h: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(g, h, add_nodes=0, fill=True)
    F = mapping_to_bijection(M)
    g = bijection_symmetric_difference(F)
    return g

def structure_difference_bijection_constant_left(
    g: graphs.structure, h: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(g, h, add_nodes=2, fill=True)
    F = mapping_to_bijection(M)
    g = bijection_symmetric_difference(F)
    return g

def structure_difference_bijection_constant_right(
    g: graphs.structure, h: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(g, h, add_nodes=3, fill=True)
    F = mapping_to_bijection(M)
    g = bijection_symmetric_difference(F)
    return g

def structure_difference_bijection_branch(
    g: graphs.structure, h: graphs.structure
) -> graphs.structure:
    M = mapper.map_to(g, h, add_nodes=3, fill=True)
    F = mapping_to_bijection(M)
    g = bijection_symmetric_difference(F)
    return g
