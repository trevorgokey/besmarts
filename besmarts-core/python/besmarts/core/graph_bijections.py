"""
besmarts.core.graph_bijections
"""

from besmarts.core import topology
import functools

# TODO convert to workspaces?
class map_vertices_ctx:
    cg = None
    o = None
    a = None
    b = None
    H = None
    strict = False
    equality = False

"""
graph_bijection_mcs(G, H) 
graph_bijection_constant_left(G, H) 
graph_bijection_constant_right(G, H) 
graph_bijection_branch(G, H) 
"""

class bijection:
    """
    A pair of graphs with an bijective node mapping between them
    """

    __slots__ = "G", "H", "T"

    def __init__(
            self,
            topo: topology.structure_topology,
            g: graphs.subgraph,
            h: graphs.subgraph,
            t: dict
        ):
        self.topology = topo
        self.g = g
        self.h = h
        self.t = t


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

def bijection_union(F: bijection):
    return graph_bitwise.union(F.g, F.h, F.t)

def bijection_union(F: bijection):
    return graph_bitwise.union(F.g, F.h, F.t)

def map_vertices_parallel(permA, permB, a, b):
    cg = map_vertices_ctx.cg
    o = map_vertices_ctx.o

    H = map_vertices_ctx.H

    strict = map_vertices_ctx.strict
    equality = map_vertices_ctx.equality

    mapping = {}
    S = 0

    valid = True
    for i, j in itertools.zip_longest(permA, permB):
        # this makes sure that if A has no node, we need to ensure
        # that B is ~[*] since that is the only way a None node
        # will be "in" B
        if i is None and j is not None and strict and not equality:
            if not o.nodes[j].all():
                valid = False
                S = -1
                break
            for nbr_j in graphs.subgraph_connection(o, j):
                edge = o.edges[tuple(sorted((j, nbr_j)))]
                if not edge.all():
                    valid = False
                    S = -1
                    break
            if not valid:
                valid = False
                S = -1
                break

        if i is None or j is None:
            continue

        edge_a = cg.edges.get(tuple(sorted((a, i))), False)
        edge_b = o.edges.get(tuple(sorted((b, j))), False)
        if strict:
            if edge_a is False or edge_b is False:
                valid = False
                S = -1
                break
            if equality:
                if not (cg.nodes[i] == o.nodes[j] and edge_a == edge_b):
                    # score_cache[(permA, permB)] = (-1, None)
                    valid = False
                    S = -1
                    break
            else:
                if not (cg.nodes[i] in o.nodes[j] and edge_a in edge_b):
                    # score_cache[(permA, permB)] = (-1, None)
                    valid = False
                    S = -1
                    break

        if edge_a is False or edge_b is False:
            edge_score = 0
        else:
            edge_score = (edge_a + edge_b).bits(maxbits=True)

        mapping[i] = j
        # add 1 so that prefer when we have two node mapping with 0 overlap
        # over the case where there was only 1 node
        # if (i,j) not in H:
        #     scores.update(pairwise_overlap(cg, sucA, o, sucB))
        if graphs.structure_node_depth(cg, i) == graphs.structure_node_depth(
            o, j
        ):
            S += H[(i, j)] + edge_score + 1
    dprint("mapped vertices:", S, mapping)
    return permA, permB, S, mapping
