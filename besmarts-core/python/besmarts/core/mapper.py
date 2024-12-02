"""
besmarts.core.mapper

Functions for mapping between two structures of arbitrary size
"""

import datetime
import math
import os
import itertools
import multiprocessing

from typing import Dict, Sequence, List, Tuple

from besmarts.core import configs, chem, graphs, db, codecs
from besmarts.core.graphs import structure
from besmarts.core import arrays
from besmarts.core.logs import dprint, timestamp
from besmarts.core import compute

# TODO: remove this shim
from besmarts.core.graphs import structure_extend as mapper_smarts_extend

mapping = Dict[int, int]


class mapped_type:
    """
    A pair of graphs with an associated node mapping between them
    """

    def __init__(self, G: graphs.subgraph, H: graphs.subgraph, map: mapping):
        self.G: graphs.subgraph = G
        self.H: graphs.subgraph = H
        self.map = map


# TODO rework to use workspaces
class union_ctx:
    A = None
    topology = None
    reference = None
    icd = None
    config = None


class align_score_ctx:
    ref = None
    to_check = None

# TODO convert to workspaces?
class map_vertices_ctx:
    cg = None
    o = None
    a = None
    b = None
    H = None
    strict = False
    equality = False

verbose=False

def mapper_invert(T: mapped_type) -> mapped_type:
    """
    Transform the map from G -> H to G <- H

    Parameters
    ----------
    T : mapped_type
        The input map

    Returns
    -------
        A new mapped type
    """
    assert T.map
    m = {v: k for k, v in T.map.items() if v is not None}

    return mapped_type(T.H, T.G, m)

def mapper_compose(T1: mapped_type, T2: mapped_type) -> mapped_type:
    """
    Transform the map from G -> H to G <- H

    Parameters
    ----------
    T : mapped_type
        The input map

    Returns
    -------
        A new mapped type
    """
    assert T1.H == T2.G
    assert T1.map
    assert T2.map
    M = {i: T2.map.get(j) for i, j in T1.map.items() if j is not None}
    return mapped_type(T1.G, T2.H, M)

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
    dprint("mapped vertices:", S, mapping, on=verbose)
    return permA, permB, S, mapping

def mapper(
    G: graphs.structure,
    H: graphs.structure,
    strict=False,
    equality=False,
    skip=None,
    add_nodes=0,
    fill=0,
    mode="high",
    return_all=False,
    pool=None,
) -> mapped_type:
    """
    Determine the mapping between two graphs G and H

    Parameters
    ----------
    G : graphs.structure
        The first input structure that defines the domain of the map
    H : graphs.structure
        The second input structure the defines the range of the map
    strict : bool
        Whether the map must satsify that G is a subset of H
    equality : bool
        Whether the map must satsify that G is equal to H
    skip : Dict[node_id, node_id]
        A set of mapped nodes that should be held constant. The mapper will try to map
        any remaining nodes.
    mode: "high" or "low"
        Whether the mapper should prefer to map nodes with "high" overlap or "low"
        overlap.

    Returns
    -------
    mapped_type
        A new mapped type
    """

    return map_to(
        G,
        H,
        strict=strict,
        equality=equality,
        add_nodes=add_nodes,
        fill=fill,
        skip=skip,
        mode=mode,
        pool=pool,
        return_all=return_all,
    )


def mapper_force_equality(G, H, pool=None) -> mapped_type:
    """
    Determine the mapping between two graphs G and H where G must be equal to H

    Parameters
    ----------
    G : graphs.structure
        The first input structure that defines the domain of the map
    H : graphs.structure
        The second input structure the defines the range of the map

    Returns
    -------
    mapped_type
        A new mapped type
    """

    return mapper(G, H, strict=True, equality=True, pool=pool)


def mapper_force_subset(G, H) -> mapped_type:
    """
    Determine the mapping between two graphs G and H where G must be a subset of H

    Parameters
    ----------
    G : graphs.structure
        The first input structure that defines the domain of the map
    H : graphs.structure
        The second input structure the defines the range of the map

    Returns
    -------
    mapped_type
        A new mapped type
    """
    return mapper(G, H, strict=True, equality=False)


def mapper_match(G, H, pool=None) -> bool:
    """
    Determine whether G is a subset of H, where missing nodes in H are assumed
    to exist but unknown, and as such are treated as full (~[*]). This allows
    a SMARTS pattern such as [#6]~[*] to be a subset of [#6]. Coincidentally,
    the opposite is true as well for this case, but wouldn't be the case for
    e.g. [#6]-[#1].

    Parameters
    ----------
    G : graphs.structure
        The first input structure that defines the domain of the map
    H : graphs.structure
        The second input structure the defines the range of the map

    Returns
    -------
    bool
        Whether G matches H
    """
    # if len(G.nodes) < len(H.nodes):
    #     return False

    Tl = mapper(
        G,
        H,
        strict=True,
        equality=False,
        add_nodes=2,
        fill=True,
        pool=pool,
        return_all=False,
    )
    Tl = [Tl]
    t = []
    for T in Tl:
        if all((x in T.map.values() for x in H.nodes)):
            t.append(T)
    Tl = t
    if not Tl:
        return False
    # print(T.map)
    # print("**")
    # graphs.subgraph_print(G)
    # print("--")
    # graphs.subgraph_print(H)
    # print("||")
    # graphs.subgraph_print(T.G)
    # print("--")
    # graphs.subgraph_print(T.H)
    cfg = None  # configs.mapper_config(0, False, "high")
    for T in Tl:
        diff = difference(T.G, T.H, cfg, map=T.map, pool=pool)
        # graphs.subgraph_print(diff)
        if not graphs.subgraph_any(diff):
            return True
    return False


def isomorphic(G: graphs.structure, H: graphs.structure) -> mapped_type:
    """
    Determine a mapping between G and H where G is isomorphic, or equal to, H

    Parameters
    ----------
    G : graphs.structure
        The first input structure that defines the domain of the map
    H : graphs.structure
        The second input structure the defines the range of the map

    Returns
    -------
    mapped_type
        The mapping
    """
    T = mapper_force_equality(G, H)
    return T


def is_isomorphic(G: graphs.structure, H: graphs.structure) -> bool:
    """
    Determine a mapping between G and H where G is isomorphic, or equal to, H

    Parameters
    ----------
    G : graphs.structure
        The first input structure that defines the domain of the map
    H : graphs.structure
        The second input structure the defines the range of the map

    Returns
    -------
    bool
    """
    T = isomorphic(G, H)
    if not T.map:
        return False
    else:
        return True


def map_to(
    cg: graphs.structure,
    o: graphs.structure,
    strict=False,
    equality=False,
    skip=None,
    return_all=False,
    add_nodes=0,
    fill=0,
    mode="high",
    pool=None,
) -> mapped_type:
    """
    Return the node mapping of one structure to another

    Parameters
    ----------
    cg : graphs.structure
        The structure to map from
    o : graphs.structure
        The structure to map to
    strict : bool
        Whether the map must satisfy that G is a subset of H
    equality : bool
        Whether the map must satisfy that G is equal to H
    skip : Dict[node_id, node_id]
        A set of mapped nodes that should be held constant. The mapper will try to map
        any remaining nodes.
    add_nodes : 0|1|2|3
        The mode for adding nodes. 
            0 means MCS
            1 means add/remove G
            2 means add/remove H
            3 means add to both
    fill=0,
    mode: str, "high" or "low"
        Whether the mapper should prefer to map nodes with "high" overlap or "low"
        overlap.

    Returns
    -------
    Dict[node_id, node_id]
        A node mapping
    """

    if skip is None:
        skip = {}

    cg_orig = cg
    o_orig = o
    cg = graphs.structure_copy(cg)
    o = graphs.structure_copy(o)

    score_cache = {}

    cg_depth_cache = {}
    o_depth_cache = {}

    # prefer adding new nodes that are already in the graph rather than
    # adding empty or full nodes
    if add_nodes > 0:
        d_cg = graphs.structure_max_depth(cg)
        d_o = graphs.structure_max_depth(o)
        # print("BOTH DEPTHS", d_cg, d_o)
        if d_cg > d_o:
            mapper_smarts_extend(
                configs.smarts_extender_config(d_cg, d_cg, True), [o]
            )
        elif d_cg < d_o:
            mapper_smarts_extend(
                configs.smarts_extender_config(d_o, d_o, True), [cg]
            )

    nbr_cg = {
        i: [x for x in graphs.subgraph_connection(cg, i) if x in cg.select]
        for i in cg.select
    }
    nbr_o = {
        i: [x for x in graphs.subgraph_connection(o, i) if x in o.select]
        for i in o.select
    }

    scores = overlap_scores(cg, o, skip=skip)

    map_scores = {}

    cg_primary = tuple([cg.select[i] for i in cg.topology.primary])
    o_primary = tuple([o.select[i] for i in o.topology.primary])
    if not all(i in skip for i in cg_primary):
        for permA in cg.topology.permutations:
            valid = True

            A = tuple((cg_primary[i] for i in permA))

            B = o_primary

            # need this to remap the connections to the permutation
            # perm_map = {i: v + 1 for i, v in enumerate(permA, 1)}

            # eh, assume same IC type for now
            # this checks for edge mapping
            # tricky edge case for pairs which have no connects
            for edge_a, edge_b in zip(cg.topology.connect, o.topology.connect):
                # edge_a = (cg.select[permA[edge_a[0]]], cg.select[edge_a[1]])
                edge_a = (A[edge_a[0]], A[edge_a[1]])
                if edge_a[0] > edge_a[1]:
                    edge_a = edge_a[::-1]

                edge_b = (o.select[edge_b[0]], o.select[edge_b[1]])
                if edge_b[0] > edge_b[1]:
                    edge_b = edge_b[::-1]

                if edge_a[0] in skip or edge_a[1] in skip:
                    continue

                if strict:
                    if equality:
                        if cg.edges[edge_a] != o.edges[edge_b]:
                            valid = False
                            break
                    else:
                        if cg.edges[edge_a] not in o.edges[edge_b]:
                            valid = False
                            break

            if not valid:
                continue

            valid = True
            mapping = {}
            S = 0
            for a, b in zip(A, B):
                if a in skip and skip[a] == b:
                    mapping[a] = b
                    S += scores[(a, b)]
                    continue

                if strict:
                    if equality:
                        if cg.nodes[a] != o.nodes[b]:
                            valid = False
                            break
                    else:
                        if cg.nodes[a] not in o.nodes[b]:
                            valid = False
                            break
                mapping[a] = b
                S += scores[(a, b)]

            if not valid:
                continue

            if S not in map_scores:
                map_scores[S] = [mapping]
            else:
                map_scores[S].append(mapping)
    else:
        map_scores[0] = [{i: skip[i] for i in cg_primary}]

    if len(map_scores) == 0:
        dprint("returned {}", on=verbose)
        ret = mapped_type(cg_orig, o_orig, {})
        if return_all:
            # this path needs to be tested better
            return [ret]
        else:
            return ret

    dprint("map_scores1", map_scores, on=verbose)

    if mode == "high":
        best_score = max(map_scores)
    elif mode == "low":
        best_score = min(map_scores)

    # preload the cache with skip if it exists; this will prevent the permutation searches
    if skip:
        lvl = 1
        while True:
            A = graphs.structure_vertices_at_depth(cg, lvl, cg_depth_cache)
            B = graphs.structure_vertices_at_depth(o, lvl, o_depth_cache)
            if len(A) == 0 or len(B) == 0:
                break
            lvl += 1
            if all(x in skip for x in A):
                for a in A:
                    for b in B:
                        sucA = [
                            i
                            for i in graphs.structure_vertices_at_depth(
                                cg, lvl, cg_depth_cache
                            )
                            if i in nbr_cg[a]
                        ]
                        sucB = [
                            i
                            for i in graphs.structure_vertices_at_depth(
                                o, lvl, o_depth_cache
                            )
                            if i in nbr_o[b]
                        ]
                        for permA in itertools.permutations(sucA):
                            for permB in itertools.permutations(sucB):
                                score_cache[(permA, permB)] = (
                                    math.inf,
                                    {a: skip[a] for a in A},
                                )

    best_maps = map_scores.pop(best_score)

    lvl = 1
    dprint("starting map_to_descend depth", lvl, on=verbose)
    total_score, total_maps = map_to_descend(
        cg,
        o,
        best_maps,
        scores,
        lvl,
        nbr_cg,
        nbr_o,
        score_cache,
        cg_depth_cache,
        o_depth_cache,
        skip=skip,
        strict=strict,
        equality=equality,
        add_nodes=add_nodes,
        fill=fill,
        pool=pool,
    )

    if total_score < 0:
        # this means we might have a suboptimal mapping that
        # is a subset; the highest overlap at the primary
        # lead to a nonsubset mapping at further depths
        for best_score in sorted(map_scores, reverse=True):
            best_maps = map_scores[best_score]
            this_score, this_maps = map_to_descend(
                cg,
                o,
                best_maps,
                scores,
                lvl,
                nbr_cg,
                nbr_o,
                score_cache,
                cg_depth_cache,
                o_depth_cache,
                skip=skip,
                strict=strict,
                equality=equality,
                add_nodes=add_nodes,
                fill=fill,
                pool=pool,
            )
            if mode == "high" and this_score > total_score:
                total_score = this_score
                total_maps = this_maps
            elif mode == "low" and this_score < total_score:
                total_score = this_score
                total_maps = this_maps
            # if total_score >= 0:
            #     break
        if total_score < 0:
            if return_all:
                # return mapped_type(cg_orig, o_orig, [{}])
                return [
                    mapped_type(cg, o, total_map) for total_map in total_maps
                ]
            else:
                return mapped_type(cg_orig, o_orig, {})

    result = []
    for i, total_map in enumerate(total_maps):
        if total_map:
            total_maps[i] = {
                k: total_map.get(k) for k in cg.select if k not in skip
            }
        if skip:
            total_maps[i].update(skip)

        total_maps[i] = {
            k: v
            for k, v in total_map.items()
            if v is not None and k in cg.select and v in o.select
        }

        cg_orig = graphs.structure_copy(cg)
        o_orig = graphs.structure_copy(o)
        if add_nodes == 0 or add_nodes == 2:
            cg.select = tuple(
                (
                    x
                    for i, x in enumerate(cg.select)
                    if x in total_map or i < len(cg.topology.primary)
                )
            )
            cg = graphs.structure_remove_unselected(cg)
            remove = [
                x
                for x in cg.select
                if x not in total_map or total_map[x] is None
            ]
            cg = graphs.structure_remove_nodes(cg, remove)
            remove = graphs.structure_unreachable_nodes(cg)
            cg = graphs.structure_remove_nodes(cg, remove)
        if add_nodes == 0 or add_nodes == 3:
            o.select = tuple(
                (
                    x
                    for i, x in enumerate(o.select)
                    if x in total_map.values() or i < len(o.topology.primary)
                )
            )
            o = graphs.structure_remove_unselected(o)
            remove = [x for x in o.select if x not in total_map.values()]
            o = graphs.structure_remove_nodes(o, remove)
            remove = graphs.structure_unreachable_nodes(o)
            o = graphs.structure_remove_nodes(o, remove)

        for k, v in list(total_map.items()):
            if k not in cg.nodes or v not in o.nodes:
                del total_map[k]

        cg.cache.clear()
        o.cache.clear()
        ans = mapped_type(cg, o, total_map)
        cg = cg_orig
        o = o_orig
        if not return_all:
            return ans
        result.append(ans)

    return result


def map_to_descend(
    cg,
    o,
    mappings,
    scores,
    lvl,
    nbr_cg,
    nbr_o,
    score_cache,
    cg_depth_cache,
    o_depth_cache,
    skip=None,
    strict=False,
    equality=False,
    add_nodes=0,
    fill=0,
    mode="high",
    return_all=False,
    pool=None,
):
    """
    Find the mappings between a subset of nodes at a certain depth from each structure.

    Parameters
    ----------
    cg : graphs.structure
        The structure to map from
    o : graphs.structure
        The structure to map to
    scores :
        The scores of mappings already completed
    scores :
        The scores of mappings already completed
    nbr_cg :
        The adjacency map of cg
    nbr_o :
        The adjacency map of o
    score_cache :
        The score cache
    cg_depth_cache :
        The precalculated depths of the nodes of cg
    o_depth_cache :
        The precalculated depths of the nodes of o
    strict : bool
        Whether the map must satsify that cg is a subset to o
    equality : bool
        Whether the map must satsify that cg is equal to o
    skip : Dict[node_id, node_id]
        A set of mapped nodes that should be held constant. The mapper will try to map
        any remaining nodes.

    Returns
    -------
    best_s :
        The score of the best mapping
    best_mapping :
        The score of best mapping
    """

    dprint("map_to_descend", mappings, "level", lvl, on=verbose)
    best_s = 0
    best_map = {}

    best_maps = []

    new_maps = {}
    if skip is None:
        skip = {}

    for idx, mapping in enumerate(mappings):
        s = 0
        total_map = mapping.copy()

        for a, b in mapping.items():
            if b is None:
                continue

            mapped_scores = map_vertices(
                cg,
                o,
                a,
                b,
                scores,
                lvl,
                nbr_cg,
                nbr_o,
                score_cache,
                cg_depth_cache,
                o_depth_cache,
                strict=strict,
                equality=equality,
                add_nodes=add_nodes,
                fill=fill,
                mode=mode,
                pool=pool,
            )
            if mode == "high":
                this_s = max(mapped_scores)
            elif mode == "low":
                this_s = min(mapped_scores)

            if this_s < 0:
                continue

            lower_s, new_maps = map_to_descend(
                cg,
                o,
                mapped_scores[this_s],
                scores,
                lvl + 1,
                nbr_cg,
                nbr_o,
                score_cache,
                cg_depth_cache,
                o_depth_cache,
                skip=skip,
                strict=strict,
                equality=equality,
                add_nodes=add_nodes,
                fill=fill,
                return_all=return_all,
                pool=pool,
            )
            if lower_s < 0:
                s = lower_s
                continue

            if new_maps:
                for new_map in new_maps:
                    for k, v in new_map.items():
                        if k not in total_map and v not in total_map.values():
                            total_map[k] = v
                        else:
                            pass
                    break

            s += this_s

        if mode == "high":
            if s > best_s:
                best_maps.clear()
            if s >= best_s:
                best_maps.append(total_map)
                best_s = s

        elif mode == "low":
            if s < best_s:
                best_maps.clear()

            if s <= best_s:
                best_maps.append(total_map)
                best_s = s

    dprint("map_to best_map", "score:", best_s, "map", best_map, on=verbose)
    return best_s, best_maps


def map_vertices(
    cg,
    o,
    a,
    b,
    scores,
    lvl,
    neighbors_cg,
    neighbors_o,
    score_cache,
    cg_depth_cache=None,
    o_depth_cache=None,
    skip=None,
    strict=False,
    equality=False,
    add_nodes=0,
    fill=0,
    mode="high",
    pool=None,
):
    """
    Return the scores of all possible mappings of given set of nodes

    Parameters
    ----------
    cg : graphs.structure
        The structure to map from
    o : graphs.structure
        The structure to map to
    a : node_id
        The node of cg to to get neighbors from
    b : node_id
        The node of o to to get neighbors from
    scores :
        The scores of mappings already completed
    lvl : int
        The depth
    nbr_cg :
        The adjacency map of cg
    nbr_o :
        The adjacency map of o
    score_cache :
        The score cache
    cg_depth_cache :
        The precalculated depths of the nodes of cg
    o_depth_cache :
        The precalculated depths of the nodes of o
    strict : bool
        Whether the map must satsify that cg is a subset to o
    equality : bool
        Whether the map must satsify that cg is equal to o
    skip : Dict[node_id, node_id]
        A set of mapped nodes that should be held constant. The mapper will try to map
        any remaining nodes.

    Returns
    -------
    Dict[int, Dict[node_id, node_id]]
        A mapping of scores to node maps
    """

    dprint(f"map_vertices begin on lvl {lvl}",  on=verbose)
    cgl = graphs.structure_max_depth(cg)
    ol = graphs.structure_max_depth(o)
    dprint(
        f"map_vertices begin on lvl {lvl} cgl {cgl} ol {ol}", on=verbose
    )

    if lvl > cgl and lvl > ol:
        return {0: [{}]}

    cg_depth_cache = None
    o_depth_cache = None
    group_mappings = {}

    if skip is None:
        skip = {}

    sucA = [
        i
        for i in graphs.structure_vertices_at_depth(
            cg, lvl, depth_cache=cg_depth_cache
        )
        if i in neighbors_cg[a]
    ]
    preA = [i for i in neighbors_cg[a] if i not in sucA]
    if skip:
        if all([i in skip for i in sucA]):
            return {0: [{i: skip[i] for i in sucA}]}
    sucB = [
        i
        for i in graphs.structure_vertices_at_depth(
            o, lvl, depth_cache=o_depth_cache
        )
        if i in neighbors_o[b]
    ]
    preB = [i for i in neighbors_o[b] if i not in sucB]

    if len(sucA) == 0 and len(sucB) == 0:
        return {0: [{}]}

    if len(sucB) == 0 and add_nodes == 0:
        cg.select = tuple((x for x in cg.select if x not in sucA))
        cg.cache.clear()

        if neighbors_cg:
            if a in neighbors_cg:
                neighbors_cg[a] = [x for x in neighbors_cg[a] if x not in sucA]
        return {0: [{a: None for a in sucA}]}

    if len(sucA) == 0 and add_nodes == 0:
        o.select = tuple((x for x in o.select if x not in sucB))
        o.cache.clear()

        if neighbors_o:
            if b in neighbors_o:
                neighbors_o[b] = [x for x in neighbors_o[b] if x not in sucB]
    if (add_nodes == 1 or add_nodes == 3) and len(sucB) < len(sucA):
        new_b = [
            x
            for i, x in enumerate(neighbors_o[b])
            if x not in o.select and x not in sucB
        ]

        # o.select = tuple(
        #     [x for x in o.select] + new_b[: (len(sucA) - len(sucB))]
        # )
        o.cache.clear()

        dprint("0ADDING nbrs:", len(sucA), len(sucB), len(new_b), on=verbose)
        for add_idx in range(
            len(sucA) - len(sucB) - len(new_b) + len(preA) - len(preB)
        ):
            # if len(neighbors_o[b]) == 4:
            #     break
            n = o.nodes[b].copy()
            ni = max(o.nodes) + 1

            if o.edges:
                ei = next(iter(o.edges))
                e = o.edges[ei].copy()
            else:
                ei = next(iter(cg.edges))
                e = cg.edges[ei].copy()

            new_ei = graphs.edge((b, ni))

            if fill == 1 or fill == 3:
                n.fill()
                e.fill()
            else:
                n.clear()
                e.clear()

            sucB.append(ni)
            o.select = tuple(list(o.select) + [ni])
            o.cache.clear()
            o.nodes[ni] = n
            o.edges[new_ei] = e

            if b not in neighbors_o:
                neighbors_o[b] = [ni]
            else:
                neighbors_o[b].append(ni)
            neighbors_o[ni] = [b]
            if len(neighbors_o[b]) > 4:
                print(
                    "0WARNING:",
                    b,
                    "has",
                    len(neighbors_o[b]),
                    "neighbors",
                    neighbors_o[b],
                    "at add",
                    add_idx,
                    "sucA sucB newb",
                    len(sucA),
                    len(sucB),
                    len(new_b),
                )

            if o_depth_cache:
                d = graphs.structure_node_depth(o, b) + 1

                nbr = None
                if d not in o_depth_cache:
                    nbr = set([ni])
                else:
                    nbr = o_depth_cache[d]
                    nbr.add(ni)
                o_depth_cache[d] = nbr

    if (add_nodes == 1 or add_nodes == 2) and len(sucB) > len(sucA):
        new_a = [
            x
            for i, x in enumerate(neighbors_cg[a])
            if x not in cg.select and x not in sucA
        ]

        # cg.select = tuple(
        #     [x for x in cg.select] + new_a[: (len(sucB) - len(sucA))]
        # )
        cg.cache.clear()
        dprint(
            "1ADDING nbrs:",
            len(sucB),
            len(sucA),
            len(new_a),
            len(preA),
            len(preB),
            on=verbose,
        )

        for add_idx in range(
            len(sucB) - len(sucA) - len(new_a) + len(preB) - len(preA)
        ):
            # if len(neighbors_cg[a]) == 4:
            #     break

            n = cg.nodes[a].copy()
            ni = max(cg.nodes) + 1

            if o.edges:
                ei = next(iter(o.edges))
                e = o.edges[ei].copy()
            else:
                ei = next(iter(cg.edges))
                e = cg.edges[ei].copy()

            new_ei = graphs.edge((a, ni))

            if fill == 1 or fill == 3:
                n.fill()
                e.fill()
            else:
                n.clear()
                e.clear()

            sucA.append(ni)
            cg.select = tuple(list(cg.select) + [ni])
            cg.cache.clear()
            cg.nodes[ni] = n
            cg.edges[new_ei] = e
            if a not in neighbors_cg:
                neighbors_cg[a] = [ni]
            else:
                neighbors_cg[a].append(ni)
            if len(neighbors_cg[a]) > 4:
                print(
                    "1WARNING:",
                    a,
                    "has",
                    len(neighbors_cg[a]),
                    "neighbors at add",
                    add_idx,
                    "sucA sucB newb",
                    len(sucB),
                    len(sucA),
                    len(new_a),
                )
            neighbors_cg[ni] = [a]

            if cg_depth_cache:
                d = graphs.structure_node_depth(cg, a) + 1

                nbr = None
                if d not in cg_depth_cache:
                    nbr = set([ni])
                else:
                    nbr = cg_depth_cache[d]
                    nbr.add(ni)
                cg_depth_cache[d] = nbr

    H = scores

    n_cached = 0
    n_calc = 0

    if len(sucA) < len(sucB):
        pairs = itertools.product(
            enumerate([sucA]),
            enumerate(itertools.permutations(sucB, len(sucA) + 1)),
        )
    elif len(sucA) > len(sucB):
        pairs = itertools.product(
            enumerate(itertools.permutations(sucA, len(sucB))),
            enumerate([sucB]),
        )
    else:
        pairs = itertools.product(
            enumerate(tuple((sucA,))),
            enumerate(itertools.permutations(sucB)),
        )
    pairs = list(pairs)
    updated_scores = False
    for (Ai, permA), (Bi, permB) in pairs:
        for i, j in itertools.zip_longest(permA, permB):
            if i is None or j is None:
                continue
            if (i, j) not in H:
                scores.update(pairwise_overlap(cg, sucA, o, sucB))
                updated_scores = True
                break
        if updated_scores:
            break

    map_vertices_ctx.cg = cg
    map_vertices_ctx.o = o
    map_vertices_ctx.H = H
    map_vertices_ctx.strict = strict
    map_vertices_ctx.equality = equality

    pool = None
    if pool is True:
        pool = multiprocessing.pool.Pool()

    work = []
    dprint(f"Number of pairs: {len(pairs)}", on=verbose)
    for (Ai, permA), (Bi, permB) in pairs:
        permA = tuple(permA)
        permB = tuple(permB)
        S = 0
        mapping = {}

        cached = score_cache.get((permA, permB), None)
        valid = True
        if cached is None:
            n_calc += 1
            if pool:
                work.append(
                    pool.apply_async(
                        map_vertices_parallel, (permA, permB, a, b)
                    )
                )
            else:
                work.append(map_vertices_parallel(permA, permB, a, b))

        else:
            n_cached += 1
            S, mapping = cached

            if S not in group_mappings:
                group_mappings[S] = [mapping]
            elif mapping not in group_mappings[S]:
                group_mappings[S].append(mapping)

    map_vertices_ctx.cg = None
    map_vertices_ctx.o = None
    map_vertices_ctx.H = None
    map_vertices_ctx.strict = None
    map_vertices_ctx.equality = None

    for unit in work:
        if pool:
            permA, permB, S, mapping = unit.get()
        else:
            permA, permB, S, mapping = unit
        score_cache[(permA, permB)] = (S, mapping)
        if S not in group_mappings:
            group_mappings[S] = [mapping]
        elif mapping not in group_mappings[S]:
            group_mappings[S].append(mapping)

    best_score = -1
    if len(group_mappings) > 0:
        if mode == "high":
            best_score = max(group_mappings)
        elif mode == "low":
            best_score = min(group_mappings)
    else:
        if pool is not None:
            pool.terminate()
            pool.close()
        return {-1: [{}]}

    if best_score < 0:
        if pool is not None:
            pool.terminate()
            pool.close()
        return {-1: [{}]}

    if pool is not None:
        pool.terminate()
        pool.close()

    best_map = group_mappings[best_score]

    best_idx = 0
    best_s = None

    # this breaks ties so we only do this if we have multiple scores
    if len(best_map) > 1 and (lvl + 1) < len(scores):
        for idx, mapping in enumerate(best_map):
            s = 0
            for x, y in mapping.items():
                if (x, y) in H:
                    continue
                new_mappings = map_vertices(
                    cg,
                    o,
                    x,
                    y,
                    scores,
                    lvl + 1,
                    neighbors_cg,
                    neighbors_o,
                    score_cache,
                    cg_depth_cache,
                    o_depth_cache,
                    strict=strict,
                    equality=equality,
                    add_nodes=add_nodes,
                    fill=fill,
                    mode=mode,
                    pool=pool,
                )
                s += max(new_mappings)
            if (
                best_s is None
                or (mode == "high" and s > best_s)
                or (mode == "low" and s < best_s)
            ):
                best_s = s
                best_idx = idx

    return {best_score: [group_mappings[best_score][best_idx]]}


def overlap_scores(cg, o, skip=None, cg_depth_cache=None, o_depth_cache=None):
    scores = {}

    if skip is None:
        skip = {}

    seen_a = set()
    seen_b = set()
    depth = min(graphs.structure_max_depth(cg), graphs.structure_max_depth(o))
    for lvl in range(depth + 1):
        A = [
            x
            for x in graphs.structure_vertices_at_depth(cg, lvl, cg_depth_cache)
            if x in cg.select and x not in seen_a
        ]
        B = [
            x
            for x in graphs.structure_vertices_at_depth(o, lvl, o_depth_cache)
            if x in o.select and x not in seen_b
        ]
        if len(A) == 0 or len(B) == 0:
            break

        H = pairwise_overlap(cg, A, o, B)
        scores.update(H)
        if skip:
            for a in A:
                if a in skip:
                    scores[(a, skip[a])] = math.inf
    return scores


def pairwise_overlap(cg, A, o, B):
    H = {}

    dprint(f"pairwise overlap {len(A)} {len(B)}", on=verbose)
    for i in A:
        prim_i = cg.nodes[i]
        bonds_i = tuple(
            tuple(sorted((i, j))) for j in graphs.subgraph_connection(cg, i)
        )
        dprint(f"pairwise overlap bonds to permute", bonds_i, on=verbose)
        if len(bonds_i) > 4:
            breakpoint()
        for j in B:
            prim_j = o.nodes[j]
            bond_j = tuple(
                tuple(sorted((j, k))) for k in graphs.subgraph_connection(o, j)
            )
            best_score = 0
            for bi, bond_i in enumerate(itertools.permutations(bonds_i), 0):
                score = 0
                for b_i, b_j in zip(bond_i, bond_j):
                    b_i = cg.edges[b_i]
                    b_j = o.edges[b_j]
                    score += (b_i & b_j).bits(maxbits=True)
                best_score = max(best_score, score)

            H[(i, j)] = (prim_i & prim_j).bits(maxbits=True) + best_score + 1
    dprint(f"pairwise overlap {A} {B} {H[(i, j)]}", on=verbose)

    return H


def group_by_isomorphism(
    structures: List[graphs.structure], mapping_cache=None
):
    """
    Group structures that are isomorphic (equal).

    Parameters
    ----------
    structures : List[graphs.structure]
        The structures to extend.

    mapping_cache : Dict[(int, int), Dict[node_id, node_id]]
        Precalculated maps for each pair of structures. Useful when you mapped at one
        depth, then extend and need to map at the next depth.

    Returns
    -------
    List[List[graphs.structure]]
        A list of groups of structures that are isomorphic (equal).
    """

    keys = list(range(len(structures)))
    groups = [keys]
    new_groups = []

    for group in groups:
        keys = list(group)
        while keys:
            n = keys.pop()
            new_group = [n]
            ref_graph = structures[n]
            for m in keys:
                cmp_graph = structures[m]
                skip = None
                if mapping_cache is not None:
                    skip = mapping_cache.get((n, m))
                mapping = map_to(
                    ref_graph, cmp_graph, strict=True, equality=True, skip=skip
                )
                if mapping:
                    new_group.append(m)
                    for k in list(mapping):
                        if mapping[k] is None:
                            del mapping[k]
                    if mapping_cache is not None:
                        mapping_cache[(n, m)] = mapping
                        mapping_cache[(m, n)] = mapping
            new_groups.append(new_group)
            for m in new_group[1:]:
                keys.remove(m)
    return new_groups

def mapper_compose_graphs(
    cg, o, map, add_nodes=False, fill_new_nodes=False
) -> Tuple[graphs.structure, graphs.structure, dict]:
    """
    Perform a bitwise operation across the nodes and edges of pair of
    structures, adding or subtracting nodes as necessary.

    Parameters
    ----------
    cg : graphs.structure
        The input structure that defines the domain of the map
    o : graphs.structure
        The input second input structures
    config: configs.mapper_config
        The configuration for mapping new nodes
    map : Dict[node_id, node_id]
        A precalculated map to use

    Returns
    -------
    structure
        The result of the operation
    """

    g: graphs.structure = graphs.structure_remove_unselected(cg)
    o = graphs.structure_remove_unselected(o)

    M = map

    for n in g.select:
        m = M.get(n)
        if m is not None:
            pass
        elif add_nodes:
            primitive = g.nodes[n]
            empty = primitive.copy()
            empty.clear()
            if fill_new_nodes:
                empty.fill()

            oidx = max(o.nodes) + 1
            o.nodes[oidx] = empty
            M[n] = oidx
            o.select = tuple((*o.select, oidx))
            o.cache.clear()

        elif n not in (g.select[i] for i in g.topology.primary):
            g = graphs.structure_remove_nodes(g, [n])
            if n in M:
                M.pop(n)

    for edge in g.edges:
        i, j = edge
        oi, oj = M[i], M[j]
        if (i not in g.select) or (j not in g.select):
            continue
        oedge = tuple(sorted((oi, oj)))
        if oedge not in o.edges:
            primitive = cg.edges[edge]
            empty = primitive.copy()
            empty.clear()
            if fill_new_nodes:
                empty.fill()
            o.edges[oedge] = empty

    Minv = {v: k for k, v in M.items() if v is not None}

    for n in o.select:
        m = Minv.get(n)
        if m is not None:
            pass
        elif add_nodes:
            primitive = o.nodes[n]
            empty = primitive.copy()
            empty.clear()
            if fill_new_nodes:
                empty.fill()

            idx = max(g.nodes) + 1
            g.nodes[idx] = empty
            M[idx] = n
            Minv[n] = idx
            g.select = tuple((*g.select, idx))
            g.cache.clear()
        elif n not in (o.select[i] for i in o.topology.primary):
            o = graphs.structure_remove_nodes(o, [n])

    for oedge in o.edges:
        oi, oj = oedge
        if (oi not in o.select) or (oj not in o.select):
            continue
        i, j = Minv[oi], Minv[oj]
        edge = tuple(sorted((i, j)))
        if edge not in g.edges:
            primitive = o.edges[oedge]
            empty = primitive.copy()
            empty.clear()
            if fill_new_nodes:
                empty.fill()
            g.edges[edge] = empty

    return g, o, M


class filter_contains_ctx:
    bes = None
    to_check = None


def filter_contains_parallel(indices):
    bes = filter_contains_ctx.bes
    to_check = filter_contains_ctx.to_check
    work_list = tuple(
        (i for i in indices if not map_to(to_check[i], bes, strict=True).map)
    )
    return work_list


def filter_contains(bes: graphs.graph, to_check, executor=None):
    N = len(to_check)
    procs = configs.processors
    chunksize = N // procs + bool(N % procs)
    indices = list(range(N))
    chunks = arrays.batched(indices, chunksize)

    filter_contains_ctx.bes = bes
    filter_contains_ctx.to_check = to_check

    masks = [filter_contains_parallel(chunk) for chunk in chunks]
    to_check = [to_check[i] for y in masks for i in y]

    filter_contains_ctx.bes = None
    filter_contains_ctx.to_check = None

    return to_check


def align_score(G: graphs.structure, H: graphs.structure):
    """
    Return the number of overlapping bits after mapping two structures

    Parameters
    ----------
    G : graphs.structure
        The first input structure
    H : graphs.structure
        The second input structure

    Returns:
    int
        The number of overlapping bits
    """

    T = map_to(G, H, add_nodes=0, fill=False)
    g = intersection(T.G, T.H, map=T.map)
    return graphs.structure_bits(g)

def align_score_parallel(indices):
    ref = align_score_ctx.ref
    to_check = align_score_ctx.to_check
    return tuple(tuple((i, align_score(ref, to_check[i]))) for i in indices)


def ordered_align_score(i, ref, o):
    return i, align_score(ref, o)


def ordered_contains(i, ref, o):
    return i, o in ref


def intersection_list(
    A: Sequence[graphs.structure],
    config: configs.mapper_config = None,
    max_depth=None,
    reference=None,
    sort=True,
    executor=None,
    verbose=False,
) -> graphs.structure:
    if config is None:
        config = configs.mapper_config(False, False, "high")

    add_nodes = config.add_nodes
    if add_nodes is False:
        max_nodes = min((graphs.structure_max_depth(a) for a in A))
        if reference:
            max_nodes = min(max_nodes, graphs.structure_max_depth(reference))

    ref = A[0]
    ref = graphs.structure_copy(A[0])

    if max_depth is not None and max_depth >= 0:
        ref = graphs.structure_up_to_depth(ref, max_depth)

    to_check = A[1:]
    total = len(to_check)

    if reference is not None:
        reference = graphs.structure_copy(reference)

    scores = []
    i = 0

    _executor = executor

    procs = configs.processors
    while len(to_check) > 0:
        if verbose:
            print(
                f"Intersection set {total-len(to_check):5d}/{total} A={len(ref.select)} Ad={graphs.structure_max_depth(ref)}{' ':15s}",
                end="\r",
            )
        if not scores:
            if sort:
                align_score_ctx.ref = ref
                align_score_ctx.to_check = to_check
                N = len(to_check)
                ck = min(10, N // procs + bool(N % procs))
                ck = max(ck, 1)
                chunks = arrays.batched(list(range(N)), ck)
                scores = []
                for chunk in chunks:
                    work = align_score_parallel(chunk)
                    scores.extend(work)
                align_score_ctx.ref = None
                align_score_ctx.to_check = None
            else:
                scores = [(0, 0)] * len(to_check)

        s = arrays.argmax([x[1] for x in sorted(scores, key=lambda y: y[0])])

        M = None
        new_g = to_check[s]
        if max_depth is not None and max_depth >= 0:
            ref = graphs.structure_up_to_depth(ref, max_depth)
            new_g = graphs.structure_up_to_depth(new_g, max_depth)
        else:
            new_g = graphs.structure_copy(new_g)

        between_map = None
        if reference is not None:
            if verbose:
                print(
                    f"Intersection ref {total-len(to_check):5d}/{total} A={len(ref.nodes)} Ad={graphs.structure_max_depth(ref)} B={len(new_g.nodes)} Bd={graphs.structure_max_depth(new_g)}{' ':15s}",
                    end="\r",
                )

            T1 = map_to(ref, reference, strict=True, add_nodes=2, fill=True)
            _, _, M1 = T1.G, T1.H, T1.map

            T2 = map_to(new_g, reference, strict=True, add_nodes=2, fill=True)
            _, _, M2 = T2.G, T2.H, T2.map

            if M1 is None:
                print(ref.select, ref.topology.primary)
                for n in ref.select:
                    print(f"{n:3d}", ref.nodes[n])
                print(reference.nodes)
                print(reference.select, reference.topology.primary)
                for n in reference.select:
                    print(f"{n:3d}", ref.nodes[n])
                raise Exception()
            if M2 is None:
                print(new_g.nodes)
                print(reference.nodes)
                raise Exception()

            M2 = {v: k for k, v in M2.items() if v is not None}
            M = {k: M2.get(v) for k, v in M1.items() if v is not None}
            M = {k: v for k, v in M.items() if v is not None}

            if verbose:
                print(
                    f"Intersection map {total-len(to_check):5d}/{total} A={len(ref.nodes)} B={len(new_g.nodes)}{' ':15s}",
                    end="\r",
                )
            T3 = map_to(ref, new_g, skip=M, add_nodes=1, fill=False, pool=True)
            ref, new_g, between_map = T3.G, T3.H, T3.map

        if verbose:
            print(
                f"Intersection exe {total-len(to_check):5d}/{total} A={len(ref.select)} Ad={graphs.structure_max_depth(ref)} B={len(new_g.select)} Bd={graphs.structure_max_depth(new_g)}",
                end="\r",
            )

        result = intersection(ref, new_g, config, map=between_map)

        scores.pop(s)
        to_check.pop(s)

        if sort and result != ref:
            # if the union did have an effect, force a rescore
            scores = []

        i += 1

        # this seems to slow things down? evaluate every once in awhile
        if i % 20 == 0:
            to_check = filter_contains(ref, to_check, executor=None)
            scores = []

        ref = result

    if verbose:
        print()
    if _executor and executor is None:
        _executor.shutdown()

    return ref


def union_list(
    A: Sequence[graphs.structure],
    config: configs.mapper_config = None,
    max_depth=None,
    reference=None,
    sort=True,
    executor=None,
    verbose=False,
) -> graphs.structure:
    if config is None:
        config = configs.mapper_config(False, False, "high")

    
    add_nodes = config.add_nodes
    if add_nodes is False:
        max_nodes = min((graphs.structure_max_depth(a) for a in A))
        if reference:
            max_nodes = min(max_nodes, graphs.structure_max_depth(reference))

    # ref = A[0]
    ref = graphs.structure_copy(A[0])

    if max_depth is not None and max_depth >= 0:
        ref = graphs.structure_up_to_depth(ref, max_depth)

    to_check = A[1:]
    total = len(to_check)

    if reference is not None:
        reference = graphs.structure_copy(reference)

    scores = []
    i = 0

    _executor = executor

    procs = configs.processors
    while len(to_check) > 0:
        if verbose:
            print(
                f"Union set {total-len(to_check):5d}/{total} A={len(ref.select)} Ad={graphs.structure_max_depth(ref)}{' ':15s}",
                end="\r",
            )
        if not scores:
            if sort:
                align_score_ctx.ref = ref
                align_score_ctx.to_check = to_check
                N = len(to_check)
                ck = min(10, N // procs + bool(N % procs))
                ck = max(ck, 1)
                chunks = arrays.batched(list(range(N)), ck)
                scores = []
                for chunk in chunks:
                    work = align_score_parallel(chunk)
                    scores.extend(work)
                align_score_ctx.ref = None
                align_score_ctx.to_check = None
            else:
                scores = [(0, 0)] * len(to_check)

        s = arrays.argmax([x[1] for x in sorted(scores, key=lambda y: y[0])])

        M = None
        new_g = to_check[s]
        if max_depth is not None and max_depth >= 0:
            ref = graphs.structure_up_to_depth(ref, max_depth)
            new_g = graphs.structure_up_to_depth(new_g, max_depth)
        else:
            new_g = graphs.structure_copy(new_g)

        between_map = None
        if reference is not None:
            if verbose:
                print(
                    f"Union ref {total-len(to_check):5d}/{total} A={len(ref.nodes)} Ad={graphs.structure_max_depth(ref)} B={len(new_g.nodes)} Bd={graphs.structure_max_depth(new_g)}{' ':15s}",
                    end="\r",
                )

            T1 = map_to(ref, reference, strict=True, add_nodes=2, fill=True)
            _, _, M1 = T1.G, T1.H, T1.map

            T2 = map_to(new_g, reference, strict=True, add_nodes=2, fill=True)
            _, _, M2 = T2.G, T2.H, T2.map

            if M1 is None:
                print(ref.select, ref.topology.primary)
                for n in ref.select:
                    print(f"{n:3d}", ref.nodes[n])
                print(reference.nodes)
                print(reference.select, reference.topology.primary)
                for n in reference.select:
                    print(f"{n:3d}", ref.nodes[n])
                raise Exception()
            if M2 is None:
                print(new_g.nodes)
                print(reference.nodes)
                raise Exception()

            M2 = {v: k for k, v in M2.items() if v is not None}
            M = {k: M2.get(v) for k, v in M1.items() if v is not None}
            M = {k: v for k, v in M.items() if v is not None}

            if verbose:
                print(
                    f"Union map {total-len(to_check):5d}/{total} A={len(ref.nodes)} B={len(new_g.nodes)}{' ':15s}",
                    end="\r",
                )
            T3 = map_to(ref, new_g, skip=M, add_nodes=1, fill=False, pool=True)
            ref, new_g, between_map = T3.G, T3.H, T3.map

        if verbose:
            print(
                f"Union exe {total-len(to_check):5d}/{total} A={len(ref.select)} Ad={graphs.structure_max_depth(ref)} B={len(new_g.select)} Bd={graphs.structure_max_depth(new_g)}",
                end="\r",
            )

        result = union(ref, new_g, config, map=between_map)

        scores.pop(s)
        to_check.pop(s)

        if sort and result != ref:
            # if the union did have an effect, force a rescore
            scores = []

        i += 1

        # this seems to slow things down? evaluate every once in awhile
        if i % 20 == 0:
            to_check = filter_contains(ref, to_check, executor=None)
            scores = []

        ref = result

    if verbose:
        print()
    if _executor and executor is None:
        _executor.shutdown()

    return graphs.structure_copy(ref)


def intersection_list_dispatch(
    indices,
) -> graphs.structure:
    A = [union_ctx.A[i] for i in indices]
    reference = union_ctx.reference
    config = union_ctx.config
    max_depth = union_ctx.max_depth

    return intersection_list(
        A, config, max_depth, reference, sort=True, executor=None, verbose=False
    )


def union_list_dispatch_distributed(
    indices: List[int],
    shm=None
) -> graphs.structure:

    topo = shm.topology

    reference = shm.reference
    config = shm.config
    max_depth = shm.max_depth

    icd: codecs.intvec_codec = shm.icd
    G, sel = shm.A

    sel = [sel[i] for i in indices]

    work = union_list_parallel(
        G,
        sel,
        topo,
        config=config,
        max_depth=max_depth,
        reference=reference,
        icd=icd,
        verbose=False
    )
    return icd.structure_encode(work)


def union_list_dispatch(
    indices: List[int],
) -> graphs.structure:
    topo = union_ctx.topology

    reference = union_ctx.reference
    config = union_ctx.config
    max_depth = union_ctx.max_depth

    if type(union_ctx.A) is db.db_dict:
        A = union_ctx.A.read_structure_list(indices)
    else:
        icd: codecs.intvec_codec = union_ctx.icd
        if union_ctx.result is None:
            G = union_ctx.A[0]
            sel = union_ctx.A[1]
            A = [graphs.graph_to_structure(icd.graph_decode(G[sel[i][0]]), sel[i][1], topo) for i in indices]
            if max_depth is not None and max_depth > 0:
                # print(f"EXTENDING to {max_depth}")
                # for i, x in enumerate(A):
                #     print(i, x.select)
                mapper_smarts_extend(configs.smarts_extender_config(max_depth, max_depth, True), A)
                # for i, x in enumerate(A):
                #     print(i, x.select)
        else:
            A = [icd.structure_decode(union_ctx.result[i]) for i in indices]
    

    
    Q = union_list(
        A, config, max_depth, reference, sort=True, executor=None, verbose=False
    )
    A = None

    return icd.structure_encode(Q)


def intersection_list_parallel(
    A: Sequence[graphs.structure],
    config: configs.mapper_config = None,
    max_depth=None,
    reference=None,
    sort=True,
    executor=None,
) -> graphs.structure:
    procs = configs.processors

    indices = list(range(len(A)))
    procs = min(os.cpu_count(), len(indices))
    procs = min(procs, configs.processors)

    if len(indices) == 1:
        return A[0]

    union_ctx.A = A
    union_ctx.reference = reference
    union_ctx.config = config
    union_ctx.max_depth = max_depth

    work = []
    while len(indices) > 1:
        with multiprocessing.pool.Pool(processes=procs) as pool:
            chunked = arrays.batched(indices, max(1, len(indices) // procs))
            work = [
                pool.apply_async(intersection_list_dispatch, (chunk,))
                for chunk in chunked
            ]
            work = [unit.get() for unit in work]
        indices = list(range(len(work)))
        union_ctx.A = work
        if len(indices) // procs < 2:
            procs = max(1, procs // 2)
    # print()

    union_ctx.A = None
    union_ctx.reference = None
    union_ctx.config = None
    union_ctx.max_depth = None

    return work[0]


def union_list_distributed(
    G: Sequence[graphs.graph],
    selections,
    topo,
    wq,
    config: configs.mapper_config = None,
    max_depth=None,
    reference=None,
    sort=True,
    executor=None,
    icd = None,
) -> graphs.structure:
    procs = configs.processors


    if len(selections) == 1:
        # return graphs.structure_copy(A[0])
        i = selections[0][0]
        sel = selections[0][1]
        if icd:
            g = graphs.graph_to_structure(icd.graph_decode(G[i]), sel, topo)
        else:
            g = graphs.graph_to_structure(G[i], sel, topo)
        return g

    # icd = None

    # if icd:
    #     print(f"{datetime.datetime.now()} Writing structures to disk...")
    #     adb = db.db_dict(icd, "A.db")
    #     adb.write_structure({i:a for i, a in enumerate(A)})
    #     union_ctx.A = adb
    # else:
    union_ctx.A = G, selections
    union_ctx.icd = icd

    if reference is None:
        reference = graphs.graph_to_structure(icd.graph_decode(G[selections[0][0]]), selections[0][1], topo)
        # union_ctx.reference = graphs.graph_to_structure(icd.graph_decode(G[selections[0]]))
    union_ctx.reference = reference
    union_ctx.result = None
    union_ctx.topology = topo
    union_ctx.config = config
    union_ctx.max_depth = max_depth

    indices = list(range(len(selections)))
    procs = min(os.cpu_count(), len(indices))
    procs = min(procs, configs.processors)

    work = None


    while len(indices) > 1:
        # print(timestamp(), f"Initializing pool")

        if len(indices) > 1000:
            print(timestamp(), f"Distributed Union merging={len(indices)}")
            shm = compute.shm_local(0, data={
                "reference": reference,
                "result": union_ctx.result,
                "topology": topo,
                "config": config,
                "max_depth": max_depth,
                "A": (G, selections),
                "icd": icd
            })
            ws = compute.workqueue_new_workspace(
                wq,
                shm=shm
            )
            
            chunk_n = max(2, len(indices) // procs)
            chunk_n = min(chunk_n, 10000)

            iterable = {
                i: ((x,), {}) for i, x in enumerate(arrays.batched(indices, chunk_n))
            }

            chunksize = 1

            work = compute.workspace_submit_and_flush(
                ws,
                union_list_dispatch_distributed,
                iterable,
                chunksize,
                1,
                len(iterable),
                verbose=True
            )
            compute.workqueue_remove_workspace(wq, ws)
            ws.close()
            work = list(work.values())
        else:
            print(timestamp(), f"Parallel Union merging={len(indices)}")

            with multiprocessing.pool.Pool(processes=procs) as pool:
                # print(timestamp(), f"Generating batches")
                chunk_n = max(2, len(indices) // procs)
                chunk_n = min(chunk_n, 10000)

                chunked = arrays.batched(indices, chunk_n)
                # print(timestamp(), f"Submitting")
                work = [
                    pool.apply_async(union_list_dispatch, (chunk,))
                    for chunk in chunked
                ]
                # print(timestamp(), f"Collecting")
                work = [unit.get() for unit in work]

        union_ctx.result = work

        # print(timestamp(), f"Done")
        indices = list(range(len(work)))
        # print(timestamp(), f"Union merging={len(indices)}")
        if len(indices) // procs < 2:
            procs = max(1, procs // 2)
    # print()
    ans = work[0]
    union_ctx.A = None
    union_ctx.reference = None
    union_ctx.config = None
    union_ctx.max_depth = None
    union_ctx.result = None
    # if icd:
    #     adb.remove()
    return icd.structure_decode(ans)


def union_list_parallel(
    G: Sequence[graphs.graph],
    selections,
    topo,
    config: configs.mapper_config = None,
    max_depth=None,
    reference=None,
    sort=True,
    executor=None,
    icd = None,
    verbose=True
) -> graphs.structure:
    procs = configs.processors


    if len(selections) == 1:
        # return graphs.structure_copy(A[0])
        i = selections[0][0]
        sel = selections[0][1]
        if icd:
            g = graphs.graph_to_structure(icd.graph_decode(G[i]), sel, topo)
        else:
            g = graphs.graph_to_structure(G[i], sel, topo)
        return g

    # icd = None

    # if icd:
    #     print(f"{datetime.datetime.now()} Writing structures to disk...")
    #     adb = db.db_dict(icd, "A.db")
    #     adb.write_structure({i:a for i, a in enumerate(A)})
    #     union_ctx.A = adb
    # else:
    union_ctx.A = G, selections
    union_ctx.icd = icd

    if reference is None:
        reference = graphs.graph_to_structure(icd.graph_decode(G[selections[0][0]]), selections[0][1], topo)
        # union_ctx.reference = graphs.graph_to_structure(icd.graph_decode(G[selections[0]]))
    union_ctx.reference = reference
    union_ctx.result = None
    union_ctx.topology = topo
    union_ctx.config = config
    union_ctx.max_depth = max_depth

    indices = list(range(len(selections)))
    procs = min(os.cpu_count(), len(indices))
    procs = min(procs, configs.processors)

    work = []

    if verbose:
        print(timestamp(), f"Union merging={len(indices)}")
    while len(indices) > 1:
        # print(timestamp(), f"Initializing pool")
        with multiprocessing.pool.Pool(processes=procs) as pool:
            # print(timestamp(), f"Generating batches")
            chunk_n = max(2, len(indices) // procs)
            chunk_n = min(chunk_n, 10000)

            chunked = arrays.batched(indices, chunk_n)
            # print(timestamp(), f"Submitting")
            work = [
                pool.apply_async(union_list_dispatch, (chunk,))
                for chunk in chunked
            ]
            # print(timestamp(), f"Collecting")
            work = [unit.get() for unit in work]

        union_ctx.result = work

        # print(timestamp(), f"Done")
        indices = list(range(len(work)))
        if verbose:
            print(timestamp(), f"Union merging={len(indices)}")
        if len(indices) // procs < 2:
            procs = max(1, procs // 2)
    # print()
    ans = work[0]
    union_ctx.A = None
    union_ctx.reference = None
    union_ctx.config = None
    union_ctx.max_depth = None
    union_ctx.result = None
    # if icd:
    #     adb.remove()
    return icd.structure_decode(ans)


def union(
    cg: structure,
    o: structure,
    config: configs.mapper_config = None,
    map=None,
    pool=None,
) -> structure:
    """
    Calculate the union of two structures

    Parameters
    ----------
    cg : structure
        The input structure that defines the domain of the map
    o : structure
        The input second input structures
    config: configs.mapper_config
        The configuration for mapping new nodes
    map : Dict[node_id, node_id]
        A precalculated map to use

    Returns
    -------
    structure
        The result of the operation
    """
    if config is None:
        config = configs.mapper_config(False, False, "high")

    return dispatch_boolean_op(
        cg, o, chem.bechem_ior, config, map=map, pool=pool
    )


def xor(self, o, config: configs.mapper_config = None, map=None):
    """
    Calculate the exclusive or (symmetric difference) of two structures

    Parameters
    ----------
    cg : structure
        The input structure that defines the domain of the map
    o : structure
        The input second input structures
    config: configs.mapper_config
        The configuration for mapping new nodes
    map : Dict[node_id, node_id]
        A precalculated map to use

    Returns
    -------
    structure
        The result of the operation
    """
    if config is None:
        config = configs.mapper_config(False, False, "high")

    return dispatch_boolean_op(self, o, chem.bechem_ixor, config, map=map)


def neg(self):
    """
    Calculate the negation of a structure

    Parameters
    ----------
    cg : structure
        The input structure that defines the domain of the map

    Returns
    -------
    structure
        The result of the operation
    """
    g = graphs.structure_copy(self)
    for n, o in g.nodes.items():
        g.nodes[n] = ~o
    for n, o in g.edges.items():
        g.nodes[n] = ~o
    return g


def subtract(
    self: structure,
    o: structure,
    config: configs.mapper_config = None,
    map=None,
):
    """
    Calculate the exclusive or (symmetric difference) of two structures

    Parameters
    ----------
    cg : structure
        The input structure that defines the domain of the map
    o : structure
        The input second input structures
    config: configs.mapper_config
        The configuration for mapping new nodes
    map : Dict[node_id, node_id]
        A precalculated map to use

    Returns
    -------
    structure
        The result of the operation
    """
    if config is None:
        config = configs.mapper_config(True, True, "high")

    ret = dispatch_boolean_op(
        self,
        o,
        chem.bechem_isubtract,
        config,
        map=map,
    )
    return ret


def subtract_conditional_right(
    self: structure,
    o: structure,
    config: configs.mapper_config = None,
    map=None,
    pool=None,
):
    """
    Calculate the exclusive or (symmetric difference) of two structures

    Parameters
    ----------
    cg : structure
        The input structure that defines the domain of the map
    o : structure
        The input second input structures
    config: configs.mapper_config
        The configuration for mapping new nodes
    map : Dict[node_id, node_id]
        A precalculated map to use

    Returns
    -------
    structure
        The result of the operation
    """
    if config is None:
        config = configs.mapper_config(True, True, "high")

    ret = dispatch_boolean_op(
        o, self, chem.bechem_subtract_conditional, config, map=map, pool=pool
    )
    return ret


def subtract_conditional_left(
    g: structure,
    h: structure,
    config: configs.mapper_config = None,
    map=None,
    pool=None,
):
    """
    Calculate the exclusive or (symmetric difference) of two structures

    Parameters
    ----------
    g : structure
        The input structure that defines the domain of the map
    h : structure
        The input second input structures
    config: configs.mapper_config
        The configuration for mapping new nodes
    map : Dict[node_id, node_id]
        A precalculated map to use

    Returns
    -------
    structure
        The result of the operation
    """
    if config is None:
        config = configs.mapper_config(True, True, "high")

    ret = dispatch_boolean_op(
        g, h, chem.bechem_subtract_conditional, config, map=map, pool=pool
    )
    return ret


def subtract_conditional(
    g: structure,
    h: structure,
    config: configs.mapper_config = None,
    map=None,
):
    """
    Calculate the exclusive or (symmetric difference) of two structures

    Parameters
    ----------
    g : structure
        The input structure that defines the domain of the map
    h : structure
        The input second input structures
    config: configs.mapper_config
        The configuration for mapping new nodes
    map : Dict[node_id, node_id]
        A precalculated map to use

    Returns
    -------
    structure
        The result of the operation
    """
    if config is None:
        config = configs.mapper_config(True, True, "high")

    ret = dispatch_boolean_op(
        g,
        h,
        chem.bechem_subtract_conditional,
        config,
        map=map,
    )
    return ret


def dispatch_boolean_op(
    cg, o, fn, config: configs.mapper_config, map=None, pool=None
) -> structure:
    """
    Perform a bitwise operation across the nodes and edges of pair of structures,
    adding or subtracting nodes as necessary.

    Parameters
    ----------
    cg : graphs.structure
        The input structure that defines the domain of the map
    o : graphs.structure
        The input second input structures
    fn : Callable
        The bitwise operation
    config: configs.mapper_config
        The configuration for mapping new nodes
    map : Dict[node_id, node_id]
        A precalculated map to use

    Returns
    -------
    graphs.structure
        The result of the operation
    """

    add_nodes = config.add_nodes
    add_nodes = False

    fill_new_nodes = config.fill_new_nodes
    mode = config.mode

    g: graphs.structure = graphs.structure_remove_unselected(cg)
    # o = graphs.structure_copy(o)

    M = map

    if M is None:
        T = mapper(
            g,
            o,
            add_nodes=config.add_nodes,
            fill=fill_new_nodes,
            mode=mode,
            pool=pool,
        )
        g, o, M = T.G, T.H, T.map
        # print(list(g.nodes.keys()))
        # print(list(o.nodes.keys()))
        # print(M)
    else:
        M = map.copy()

    if M is None:
        breakpoint()
        M = mapper(g, o, None, pool=pool).map
        return None

    idx = max(max(g.nodes), max(o.nodes)) + 1

    for n in list(reversed(sorted(o.nodes))):
        v = o.nodes.pop(n)
        o.nodes[n + idx] = v

    for i, j in list(o.edges):
        v = o.edges.pop((i, j))
        o.edges[(i + idx, j + idx)] = v

    o.select = tuple((n + idx for n in o.select))
    o.cache.clear()

    for m, n in list(M.items()):
        M[m] = n + idx

    # performs op on the mapped nodes
    for n in list(M):
        m = M[n]
        if m is not None:
            g.nodes[n] = fn(g.nodes[n], o.nodes[m])

        elif n not in (g.select[i] for i in g.topology.primary):
            graphs.structure_remove_nodes(g, [n])

    for edge in g.edges:
        i, j = edge
        if i not in M or j not in M:
            continue
        if M[i] is None or M[j] is None:
            continue

        mapped_edge = tuple(sorted((M[i], M[j])))
        edge_exists = mapped_edge in o.edges  # edges()
        if edge_exists:
            g.edges[edge] = fn(g.edges[edge], o.edges[mapped_edge])

    return g

def intersection(
    cg: structure,
    o: structure,
    config: configs.mapper_config = None,
    map=None,
    pool=None,
):
    """
    Calculate the intersection of two structures

    Parameters
    ----------
    cg : structure
        The input structure that defines the domain of the map
    o : structure
        The input second input structures
    config: configs.mapper_config
        The configuration for mapping new nodes
    map : Dict[node_id, node_id]
        A precalculated map to use

    Returns
    -------
    structure
        The result of the operation
    """
    if config is None:
        config = configs.mapper_config(True, True, "high")

    return dispatch_boolean_op(
        cg, o, chem.bechem_iand, config, map=map, pool=pool
    )


def intersection_conditional(
    cg: structure,
    o: structure,
    config: configs.mapper_config = None,
    map=None,
    pool=None,
):
    """
    Calculate the intersection of two structures

    Parameters
    ----------
    cg : structure
        The input structure that defines the domain of the map
    o : structure
        The input second input structures
    config: configs.mapper_config
        The configuration for mapping new nodes
    map : Dict[node_id, node_id]
        A precalculated map to use

    Returns
    -------
    structure
        The result of the operation
    """
    if config is None:
        config = configs.mapper_config(True, True, "high")

    return dispatch_boolean_op(
        cg, o, chem.bechem_and_conditional, config, map=map, pool=pool
    )


def difference(
    cg: graphs.subgraph,
    o: graphs.subgraph,
    config: configs.mapper_config = None,
    map=None,
    pool=None,
):
    """
    Calculate the difference of two structures

    Parameters
    ----------
    cg : structure
        The input structure that defines the domain of the map
    o : structure
        The input second input structures
    config: configs.mapper_config
        The configuration for mapping new nodes
    map : Dict[node_id, node_id]
        A precalculated map to use

    Returns
    -------
    structure
        The result of the operation
    """
    if config is None:
        config = configs.mapper_config(True, True, "high")

    return dispatch_boolean_op(
        cg, o, chem.bechem_subtract, config, map=map, pool=pool
    )

