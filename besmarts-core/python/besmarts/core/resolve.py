"""
besmarts.core.resolve

Resolve a SMARTS pattern into one or more SMILES
"""


from typing import Sequence, List, Tuple
import itertools

from besmarts.core import graphs
from besmarts.core.rulesets import visitor_ruleset


def resolve_smiles_recurse(
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
                        yield P
                P = graphs.graph_copy(graphs.graph(nodes, edges))
                yield from resolve_smiles_recurse(
                    beg,
                    beg_adj,
                    visitor_rulesets,
                    frag_cache,
                    P.nodes,
                    P.edges,
                    {i: j.copy() for i, j in adj.items()},
                    seen,
                )


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
    return edge_frags


def hash_fragment(nodes, edges):
    h = hash(
        (
            tuple(nodes.keys()),
            tuple((hash(v) for v in nodes.values())),
            tuple(edges.keys()),
            tuple(nodes.keys()),
            tuple((hash(v) for v in edges.values())),
        )
    )
    return h


def resolve_smiles(
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
            print(len(chem.to_fragments()))
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
