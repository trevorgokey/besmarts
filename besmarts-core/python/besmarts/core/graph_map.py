
"""
besmarts.core.graph_map
"""

class mapped_type:
    """
    A pair of graphs with an associated node mapping between them
    """

    def __init__(self, G: graphs.subgraph, H: graphs.subgraph, map: mapping):
        self.G: graphs.subgraph = G
        self.H: graphs.subgraph = H
        self.map = map

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
        dprint("returned {}")
        ret = mapped_type(cg_orig, o_orig, {})
        if return_all:
            # this path needs to be tested better
            return [ret]
        else:
            return ret

    dprint("map_scores1", map_scores)

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
    dprint("starting map_to_descend depth", lvl)
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

    dprint("map_to_descend", mappings, "level", lvl, on=add_nodes == 2)
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

            mapped_scores = map_nodes(
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

    dprint("map_to best_map", "score:", best_s, "map", best_map)
    return best_s, best_maps

def map_nodes(
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

    dprint(f"map_vertices begin on lvl {lvl}", on=add_nodes == 2)
    cgl = graphs.structure_max_depth(cg)
    ol = graphs.structure_max_depth(o)
    dprint(
        f"map_vertices begin on lvl {lvl} cgl {cgl} ol {ol}", on=add_nodes == 2
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

        dprint("0ADDING nbrs:", len(sucA), len(sucB), len(new_b))
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
            on=True,
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
    dprint(f"Number of pairs: {len(pairs)}", on=add_nodes == 2)
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
                        map_nodes_parallel, (permA, permB, a, b)
                    )
                )
            else:
                work.append(map_nodes_parallel(permA, permB, a, b))

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
                new_mappings = map_nodes(
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

    dprint(f"pairwise overlap {len(A)} {len(B)}")
    for i in A:
        prim_i = cg.nodes[i]
        bonds_i = tuple(
            tuple(sorted((i, j))) for j in graphs.subgraph_connection(cg, i)
        )
        dprint(f"pairwise overlap bonds to permute", bonds_i)
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
    dprint(f"pairwise overlap {A} {B} {H[(i, j)]}")

    return H

def map_nodes_parallel(permA, permB, a, b):
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
