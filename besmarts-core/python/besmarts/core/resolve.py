"""
besmarts.core.enumerate

Enumerate a SMARTS pattern into one or more SMILES
"""


from typing import Sequence, List, Tuple
import itertools

from besmarts.core import graphs
from besmarts.core import compute
from besmarts.core.rulesets import visitor_ruleset


def resolve_smiles_recurse_iter(
    beg: graphs.graph,
    beg_adj,
    visitor_rulesets: List[visitor_ruleset],
    frag_cache,
    nodes,
    edges,
    adj,
    seen,
):
    for idx, node in beg.nodes.items():
        if idx in nodes:
            continue

        edge_frags = get_edge_frags(beg, edges, beg_adj, adj, idx)

        frags = frag_cache.get(hash(node), None)
        if frags is None:
            frags = [
                frag
                for frag in node.to_fragments()
                if all(
                    (ruleset.on_node(idx, frag) for ruleset in visitor_rulesets)
                )
            ]

        for frag in frags:
            nodes[idx] = frag

            # molecule-so-far based rules
            if not all(ruleset.on_nodes(nodes) for ruleset in visitor_rulesets):
                nodes.pop(idx)
                continue

            nodes.pop(idx)

            for edge_set in itertools.product(*edge_frags):
                # atom and bonds rules
                if not all(
                    ruleset.on_node_edges(idx, frag, edge_set)
                    for ruleset in visitor_rulesets
                ):
                    continue

                for eidx, e in edge_set:
                    edges[eidx] = e

                # print("ACCEPTED")
                nodes[idx] = frag

                frag_hash = hash_fragment(nodes, edges)

                if frag_hash in seen:
                    continue

                seen.add(frag_hash)
                if len(nodes) == len(beg.nodes) and len(edges) == len(
                    beg.edges
                ):
                    if all(
                        (
                            ruleset.on_graph(nodes, edges)
                            for ruleset in visitor_rulesets
                        )
                    ):
                        P = graphs.graph_copy(graphs.graph(nodes, edges))
                        # print("-----------")
                        # print(P.nodes)
                        # print(P.edges)
                        # print("***********")
                        yield P
                P = graphs.graph_copy(graphs.graph(nodes, edges))
                yield from resolve_smiles_recurse_iter(
                    beg,
                    beg_adj,
                    visitor_rulesets,
                    frag_cache,
                    P.nodes,
                    P.edges,
                    {i: j.copy() for i, j in adj.items()},
                    seen,
                )


def resolve_process_fragment(
    beg,
    beg_adj,
    visitor_rulesets,
    frag_cache,
    edge_frags,
    nodes,
    edges,
    adj,
    seen,
    idx,
    gcd,
    frag,
    shm=None,
):

    # print("Enter RPF")
    results = []

    nodes[idx] = frag

    # molecule-so-far based rules
    if not all(ruleset.on_nodes(nodes) for ruleset in visitor_rulesets):
        nodes.pop(idx)
        return []

    nodes.pop(idx)

    prods = list(itertools.product(*edge_frags))
    # nprod = list(itertools.product(range(0, len(edge_frags))))
    for i, edge_set in enumerate(prods, 0):
        # atom and bonds rules
        # print(" "*len(nodes), f"Enumerate atom {idx} cnd={frag} N={len(nodes)} edgeset {i+1}/{len(prods)} {nprod[0]} {edge_set}")
        if not all(
            ruleset.on_node_edges(idx, frag, edge_set)
            for ruleset in visitor_rulesets
        ):
            continue

        for eidx, e in edge_set:
            edges[eidx] = e

        # print("ACCEPTED")
        nodes[idx] = frag

        frag_hash = hash_fragment(nodes, edges)

        if frag_hash in seen:
            continue

        seen.add(frag_hash)
        if len(nodes) == len(beg.nodes) and len(edges) == len(
            beg.edges
        ):
            if all(
                (
                    ruleset.on_graph(nodes, edges)
                    for ruleset in visitor_rulesets
                )
            ):
                P = graphs.graph(nodes, edges)
                # print("-----------")
                # print(P.nodes)
                # print(P.edges)
                # print("***********")
                results.append(gcd.smiles_encode(P))
        P = graphs.graph_copy(graphs.graph(nodes, edges))
        results.extend((
            gcd.smiles_encode(x)
            for x in resolve_smiles_recurse_iter(
                beg,
                beg_adj,
                visitor_rulesets,
                frag_cache,
                P.nodes,
                P.edges,
                {i: j.copy() for i, j in adj.items()},
                seen,
            )
        ))
    return results


def resolve_smiles_recurse_distributed(
    beg: graphs.graph,
    beg_adj,
    visitor_rulesets: List[visitor_ruleset],
    frag_cache,
    nodes,
    edges,
    adj,
    seen,
    ws,
    gcd,
    library_atoms=None
):
    for idx, node in beg.nodes.items():
        if idx in nodes:
            continue

        # print(f"Visit atom {idx}/{len(beg.nodes)}")
        edge_frags = get_edge_frags(beg, edges, beg_adj, adj, idx)

        frags = frag_cache.get(hash(node), None)
        if frags is None:
            if library_atoms is None:
                frags = [
                    frag
                    for frag in node.to_fragments()
                    if all(
                        (ruleset.on_node(idx, frag) for ruleset in visitor_rulesets)
                    )
                ]
            else:
                frags = [
                    frag
                    for frag in gcd.smarts_decode(lib).nodes[1].to_fragments()
                    for lib in library_atoms
                    if all(
                        (ruleset.on_node(idx, frag) for ruleset in visitor_rulesets)
                    )
                ]
            # frags = [
            #     frag
            #     for frag in node.to_fragments()
            #     if all(
            #         (ruleset.on_node(idx, frag) for ruleset in visitor_rulesets)
            #     )
            # ]

        args = (
            beg,
            beg_adj,
            visitor_rulesets,
            frag_cache,
            edge_frags,
            nodes,
            edges,
            adj,
            seen,
            idx,
            gcd
        )

        iterable = {
            i: [(*args, frag), {}]
            for i, frag in enumerate(frags)
        }
        if ws:
            ws.reset()
            results = compute.workspace_submit_and_flush(
                ws,
                resolve_process_fragment,
                iterable,
                verbose=True
            )
        else:
            results = {
                i: resolve_process_fragment(*args[0], **args[1])
                for i, args in iterable.items()
            }

        for k, v in results.items():
            if v:
                yield from v


def get_edge_frags(beg, edges, beg_adj, adj, idx):
    if not beg_adj:
        return []
    adj_edges: Sequence[Tuple[int, int]] = adj.get(idx)
    if adj_edges is None:
        adj_edges = []
        for j in beg_adj[idx]:
            new_edge = (idx, j) if idx < j else (j, idx)
            if new_edge not in adj_edges:
                adj_edges.append(new_edge)
        adj[idx] = adj_edges

    new_edges: Sequence[Tuple[int, int]] = [
        e for e in adj_edges if e not in edges
    ]
    old_edges: Sequence[Tuple[int, int]] = [e for e in adj_edges if e in edges]
    edge_frags = []
    edge: Tuple[int, int]
    for edge in new_edges:
        edge_frags.append(
            list(
                zip(
                    itertools.repeat(edge),
                    [frag for frag in beg.edges[edge].to_fragments()],
                )
            )
        )
    old_edge_frags = []
    for edge in old_edges:
        old_edge_frags.append([(edge, edges[edge])])
    if old_edge_frags:
        edge_frags += old_edge_frags
    # print(f"Returned {len(edge_frags)} fragments")
    return edge_frags


def hash_fragment(nodes, edges):
    h = hash(
        (
            tuple(sorted(nodes.keys())),
            tuple(sorted((hash(v) for v in nodes.values()))),
            tuple(sorted(edges.keys())),
            # tuple(nodes.keys()),
            tuple(sorted((hash(v) for v in edges.values()))),
        )
    )
    return h


def resolve_smiles_iter(
    beg: graphs.graph, visitor_rulesets: List[visitor_ruleset]
):
    for ruleset in visitor_rulesets:
        ruleset.on_start(beg)

    beg_adj = graphs.graph_connections(beg)

    # pre generate the atom fragments and filter them
    # this allows reuse of atom hashes that have the
    # same smarts, such as * atoms
    frag_cache = {}
    for idx, chem in beg.nodes.items():
        h = hash(chem)
        if h not in frag_cache:
            frags = [
                frag
                for frag in chem.to_fragments()
                if all(
                    (ruleset.on_node(idx, frag) for ruleset in visitor_rulesets)
                )
            ]
            # print(len(chem.to_fragments()))
            frag_cache[hash(chem)] = frags

    for chem in beg.edges.values():
        h = hash(chem)
        if h not in frag_cache:
            frag_cache[hash(chem)] = chem.to_fragments()

    yield from resolve_smiles_recurse(
        beg, beg_adj, visitor_rulesets, frag_cache, {}, {}, {}, set()
    )

    for ruleset in visitor_rulesets:
        ruleset.on_stop(beg, None)


def resolve_smiles(
    beg: graphs.graph, visitor_rulesets: List[visitor_ruleset]
):
    ret = list(resolve_smiles_iter(beg, visitor_rulesets))
    return ret


def resolve_smiles_distributed(
    beg: graphs.graph, visitor_rulesets: List[visitor_ruleset], ws, gcd, library_atoms=None
):
    for ruleset in visitor_rulesets:
        ruleset.on_start(beg)

    beg_adj = graphs.graph_connections(beg)

    # pre generate the atom fragments and filter them
    # this allows reuse of atom hashes that have the
    # same smarts, such as * atoms
    frag_cache = {}
    for idx, chem in beg.nodes.items():
        h = hash(chem)
        if h not in frag_cache:
            if library_atoms is None:
                frags = [
                    frag
                    for frag in chem.to_fragments()
                    if all(
                        (ruleset.on_node(idx, frag) for ruleset in visitor_rulesets)
                    )
                ]
            else:
                frags = [
                    frag
                    for lib in library_atoms
                    for frag in gcd.smarts_decode(lib).nodes[1].to_fragments()
                    if all(
                        (ruleset.on_node(idx, frag) for ruleset in visitor_rulesets)
                    )
                ]
                # print(frags)
            # print(len(chem.to_fragments()))
            frag_cache[hash(chem)] = frags

    for chem in beg.edges.values():
        h = hash(chem)
        if h not in frag_cache:
            frag_cache[hash(chem)] = chem.to_fragments()

    yield from resolve_smiles_recurse_distributed(
        beg, beg_adj, visitor_rulesets, frag_cache, {}, {}, {}, set(), ws, gcd, library_atoms=library_atoms
    )

    for ruleset in visitor_rulesets:
        ruleset.on_stop(beg, None)
