"""
besmarts.core.graph_bitwise

Bitwise operations on BESMARTS graphs. All functions require a bijection.
Additionally, we are not allowed to

"""

from besmarts.core import graphs
from besmarts.core import graph_bijections
from besmarts.core import graph_multijections

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

def intersection_multijection(MF):
    pass

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


def union_multijection(FF):
    pass

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
) -> graphs.structure:
    procs = configs.processors


    if len(selections) == 1:
        # return graphs.structure_copy(A[0])
        i = selections[0][0][0]
        sel = selections[0][1]
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
    cg: graphs.structure,
    o: graphs.structure,
    config: configs.mapper_config = None,
    map=None,
    pool=None,
) -> graphs.structure:
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

def intersection_bijection():
    pass

def intersection_bijection(f: graph_bijections.bijection):
    pass

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

