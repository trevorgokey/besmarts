"""
besmarts.core.hierarchy_merge
"""

import pprint
import functools
import multiprocessing
import multiprocessing.pool
from typing import List

from besmarts.core import hierarchies
from besmarts.core import trees
from besmarts.core import codecs
from besmarts.core import assignments

from besmarts.core import tree_iterators
from besmarts.core import mapper
from besmarts.core import configs

from besmarts.core import graphs

from besmarts.core.mapper import intersection, union, mapper_match


def structure_hierarchy_prune_unused(
    hA: hierarchies.structure_hierarchy,
    gcd,
    labeler: assignments.smarts_hierarchy_assignment,
    smi_list: List,
) -> hierarchies.structure_hierarchy:
    topo = hA.topology
    lA: assignments.smiles_assignment_group = labeler.assign(
        hA, gcd, smi_list, topo
    )
    seen = set()
    for assignment in lA.assignments:
        seen.update(set(assignment.selections.values()))

    roots = [hA.index.nodes[i] for i, x in hA.index.above.items() if x is None]
    for root in roots:
        nodes = list(tree_iterators.tree_iter_dive_reverse(hA.index, root))
        for n in nodes:
            if n.index in hA.index.nodes and n.name not in seen:
                hA.index.node_remove(n.index)
                hA.smarts.pop(n.index)
                hA.subgraphs.pop(n.index)
    return hA


def structure_hierarchy_fit(
    hA: hierarchies.structure_hierarchy,
    gcd,
    labeler: assignments.smarts_hierarchy_assignment,
    smi_list: List,
    min_depth=0,
) -> hierarchies.structure_hierarchy:
    topo = hA.topology
    lA: assignments.smiles_assignment_group = labeler.assign(
        hA, gcd, smi_list, topo
    )
    groups = {}
    for assignment in lA.assignments:
        g = gcd.smiles_decode(assignment.smiles)
        for sel, lbl in assignment.selections.items():
            if lbl not in groups:
                groups[lbl] = []
            groups[lbl].append(graphs.graph_to_structure(g, sel, topo))
    for idx, node in hA.index.nodes.items():
        lbl = node.name
        A = groups.get(lbl)
        if A:
            S0 = graphs.subgraph_to_structure(hA.subgraphs[idx], topo)
            depth = graphs.structure_max_depth(S0)
            depth = max(depth, min_depth)
            mapper.mapper_smarts_extend(
                configs.smarts_extender_config(depth, depth, True), A
            )
            Q = mapper.union_list_parallel(A, reference=S0, max_depth=depth)
            Q = mapper.intersection(Q, S0, configs.mapper_config(0, 1, "high"))
            relabel = {x: i for i, x in enumerate(Q.select, 1)}
            for i, x in enumerate(Q.nodes, len(Q.select) + 1):
                if x not in Q.select:
                    relabel[x] = i
            Q = graphs.structure_relabel_nodes(Q, relabel)
            hA.subgraphs[idx] = Q
            hA.smarts[idx] = gcd.smarts_encode(Q)
    return hA


def return_g(x):
    return x


def take_left(a, b):
    return a


def take_right(a, b):
    return b


def structure_hierarchy_merge(
    hA: hierarchies.structure_hierarchy,
    hB: hierarchies.structure_hierarchy,
    gcd,
):
    topo = hA.topology

    assert hA.topology == hB.topology

    # force_counts = hindex_assign_counts(hA, hB, gcd, smi_list)

    # PRUNE = (smi_list is not None) and smi_list
    # force_counts = {'Bonds': force_counts['Bonds']}

    # The marginals inform whether there is a parameter in either
    # force that is not covered by intersections, A > A&B or B > A&B
    # marginals = {}
    # for force, counts in force_counts.items():
    #     marginals[force] = {}
    #     for (lblA, lblB), n in counts.items():
    #         if lblA is None:
    #             if lblB is None:
    #                 marginals[force][(None, None)] = (
    #                     marginals[force].get((None, None), 0) + n
    #                 )
    #             else:
    #                 marginals[force][(None, lblB)] = (
    #                     marginals[force].get((None, lblB), 0) + n
    #                 )
    #         else:
    #             if lblB is None:
    #                 marginals[force][(lblA, None)] = (
    #                     marginals[force].get((lblA, None), 0) + n
    #                 )
    #             else:
    #                 marginals[force][(lblA, None)] = (
    #                     marginals[force].get((lblA, None), 0) - n
    #                 )
    #                 marginals[force][(None, lblB)] = (
    #                     marginals[force].get((None, lblB), 0) - n
    #                 )
    # print("Marginals")
    # pprint.pprint(marginals)

    # print("Counts")
    # pprint.pprint(force_counts)

    # hidx = hierarchies.structure_hierarchy_copy(hA)
    hidx = hierarchies.structure_hierarchy(trees.tree_index(), {}, {}, topo)
    # hidx.db = {"SMIRNOFF": {"version": "0.3", "aromaticity_model": "OEAroModel_MDL"}}
    # for idx in hA.root.down:
    #     if idx in hA.force_entry_idx.values():
    #         continue
    #     hent = hidx.hentry_add(0, hentry())
    #     hent.key = hA[idx].key
    #     attrib = hA.db[hent.key]
    #     hidx.db[hent.key] = dict(attrib)

    # opers = {
    #     "|": mapper.union,
    #     "l": take_left,
    #     "r": take_right,
    #     ":-": mapper.subtract_conditional_left,
    #     "-:": mapper.subtract_conditional_right,
    #     "^": mapper.intersection,
    # }
    opers = {
        "|": functools.partial(
            mapper.union, config=configs.mapper_config(1, False, "high")
        ),
        "l": take_left,
        "r": take_right,
        ":-": functools.partial(
            mapper.subtract_conditional_left,
            config=configs.mapper_config(3, False, "high"),
        ),
        "-:": functools.partial(
            mapper.subtract_conditional_right,
            config=configs.mapper_config(2, False, "high"),
        ),
        "^": functools.partial(
            mapper.intersection, config=configs.mapper_config(1, True, "high")
        ),
    }
    # mergers = {
    #     "|": parameter_terms_merge,
    #     "l": parameter_terms_merge_left,
    #     "r": parameter_terms_merge_right,
    #     ":-": parameter_terms_merge_left,
    #     "-:": parameter_terms_merge_right,
    #     "^": parameter_terms_merge,
    # }

    rootsA = [hA.index.nodes[i] for i, x in hA.index.above.items() if x is None]
    # next((x for x in hA.entries.values() if x.key == force))
    rootsB = [hB.index.nodes[i] for i, x in hB.index.above.items() if x is None]
    # rootB = hB.index.nodes[
    #     0
    # ]  # next((x for x in hB.entries.values() if x.key == force))

    print(f"Forming product on <{len(rootsA)}|{len(rootsB)}>")
    # hidx.force_entry_idx[force] = fhent.index
    # hidx.db[force] = dict(hA.db[rootA.key])

    hentA = None

    # fhent = trees.tree_node(0, "parameter", "", "NULL")
    # fhent = hidx.index.node_add(None, fhent)
    # hidx.smarts[fhent.index] = None
    # hidx.subgraphs[fhent.index] = None

    total = len(list(tree_iterators.tree_iter_dive(hA.index, rootsA)))
    for sym, oper in opers.items():
        seen = set()
        count = 1
        for rootA in rootsA:
            fhent = hidx.index.node_add(None, trees.tree_node_copy(rootA))
            hidx.smarts[fhent.index] = hA.smarts[rootA.index]
            hidx.subgraphs[fhent.index] = hA.subgraphs[rootA.index]
            entriesA = list(tree_iterators.tree_iter_dive(hA.index, rootA))
            # hentB = next(tree_iter_breadth_first(hB, rootB))

            # mapping = {rootA.index: fhent.index}
            for ei, hentA in enumerate(entriesA, 1):
                print(
                    f"{count:4d}/{total:4d}",
                    "Param",
                    hentA.name,
                    "Operation",
                    sym,
                )
                count += 1
                # fhent = hidx[mapping[hentA.up]]

                # breakpoint()
                for rootB in rootsB:
                    # for bidx in hB.index.below[rootB.index]:
                    # hentB = hB.index.nodes[bidx]
                    # if hentB.index == 132:
                    #     breakpoint()

                    hidx = structure_hierarchy_add_hierarchy(
                        hidx, fhent, hB, rootB
                    )
                    # hidx = tree_merge(hidx, fhent, hB, hentB)

                new_desc = list(
                    tree_iterators.tree_iter_breadth_first(hidx.index, fhent)
                )

                # this will take the first node of the just-added index from hB

                # now go into this node
                # is is important not to start from fhent, as fhent has
                # descendents we already processed, and want to skip. This
                # is why we take the last descendent above

                # print("Composing keys for", new_desc)
                # print("Seen", seen)
                work = {}
                if sym in "lr" or configs.processors == 1:
                    _pool = multiprocessing.pool.ThreadPool
                else:
                    _pool = multiprocessing.Pool
                # with multiprocessing.pool.ThreadPool(configs.processors) as pool:
                with _pool(configs.processors) as pool:
                    for hent in new_desc:
                        if hent.name in seen:
                            # print("Continue", hent.name)
                            continue
                        # print("Consider", hent.name)
                        g1 = hA.subgraphs[hentA.index]  # .db[lblA].get("graph")
                        g2 = hidx.subgraphs[
                            hent.index
                        ]  # .db[lblA].get("graph")

                        g = None
                        if g1:
                            # print("    A:", gcd.smarts_encode(g1))
                            if g2:
                                work[hent.index] = (
                                    hentA.index,
                                    pool.apply_async(
                                        oper,
                                        (
                                            graphs.subgraph_to_structure(
                                                g1, topo
                                            ),
                                            graphs.subgraph_to_structure(
                                                g2, topo
                                            ),
                                        ),
                                    ),
                                )
                                # g = graphs.structure_to_subgraph(g)
                                # print("    B:", gcd.smarts_encode(g2))
                            else:
                                work[hent.index] = (
                                    hentA.index,
                                    pool.apply_async(return_g, (g1,)),
                                )
                                g = g1
                        elif g2:
                            # print("    B:", gcd.smarts_encode(g2))
                            work[hent.index] = (
                                hentA.index,
                                pool.apply_async(return_g, (g2,)),
                            )
                            # g = g2
                        # print(hent.key, counts.get(hent.key, 0), gcd.smarts_encode(g))
                    for idx, (idxb, unit) in work.items():
                        if idx in seen:
                            continue
                        g = unit.get()
                        hent = hidx.index.nodes[idx]
                        if g is not None and graphs.graph_is_null(g):
                            # print("Removing dead parameter", hent.name)
                            hidx.index.node_remove(idx)
                            hidx.subgraphs.pop(idx)
                            hidx.smarts.pop(idx)
                        else:
                            hentA = hA.index.nodes[idxb]
                            hent.name = (str(hentA.name), sym, str(hent.name))
                            seen.add(hent.name)
                            hidx.subgraphs[hent.index] = g
                            hidx.smarts[hent.index] = gcd.smarts_encode(g)

            hidx.smarts.pop(fhent.index)
            hidx.subgraphs.pop(fhent.index)
            hidx.index.node_remove(fhent.index)

    for idx in hidx.index.nodes:
        hidx.smarts[idx] = gcd.smarts_encode(hidx.subgraphs[idx])

    hierarchies.smarts_hierarchy_print(hidx)

    roots = [
        hidx.index.nodes[i] for i, x in hidx.index.above.items() if x is None
    ]
    removed = set()
    print("Removing local occlusions")
    work = {}
    nodes = []
    for root in roots:
        nodes += list(tree_iterators.tree_iter_dive_reverse(hidx.index, root))

    with multiprocessing.Pool(configs.processors) as pool:
        for ai, a in enumerate(nodes, 1):
            break
            # if a.index in removed:
            #     continue

            if hidx.index.above[a.index] is None:
                continue

            ga = hidx.subgraphs[a.index]
            if not ga:
                continue
            ga = graphs.subgraph_to_structure(ga, topo)

            gb = hidx.subgraphs[hidx.index.above[a.index]]
            # gb = hidx.db[hidx[a.up].key].get("graph")
            if not gb:
                continue

            gb = graphs.subgraph_to_structure(gb, topo)
            # print("Consider", a, hidx.db[a.key].get('smirks'), gcd.smarts_encode(ga))
            work[a.index] = pool.apply_async(mapper_match, (gb, ga))

        for idx, matched in work.items():
            print(
                f"{idx+len(removed):4d}/{len(nodes):4d}",
                "Visiting pruner",
                idx,
                hidx.index.nodes[idx].name,
            )

            if idx in hidx.index.nodes:
                if matched.get():
                    # print("Removing due to duplicate:", a, hidx.db[a.key]['smirks'], "b:", b, hidx.db[hidx[a.up].key]['smirks'])
                    # print("Removing local occlusion due to duplicate:", a, a.key)
                    # hidx.tree.hentry_remove(a.index)
                    hidx.index.node_remove(idx)
                    hidx.subgraphs.pop(idx)
                    hidx.smarts.pop(idx)
                    # removed.add(idx)
        pool.terminate()

    removed = set()
    roots = [
        hidx.index.nodes[i] for i, x in hidx.index.above.items() if x is None
    ]
    nodes = list(tree_iterators.tree_iter_dive(hidx.index, roots))
    ordering = {x.index: i for i, x in enumerate(nodes)}
    # print("Removing nonlocal occlusions")
    for ai, a in enumerate(nodes, 1):
        break
        print(
            f"{ai+len(removed):4d}/{len(nodes):4d}",
            "Visiting pruner",
            a.index,
            a.name,
        )
        if a.index in removed:
            continue
        # ga = hidx.db[a.key].get("graph")
        ga = hidx.subgraphs[a.index]
        if not ga:
            continue
        ga = graphs.subgraph_to_structure(ga, topo)
        # print("Consider", a, hidx.db[a.key].get('smirks'), gcd.smarts_encode(ga))
        work = {}
        with multiprocessing.Pool(configs.processors) as pool:
            rootsB = [
                hidx.index.nodes[i]
                for i, x in hidx.index.above.items()
                if x is None
            ]
            for b in list(
                tree_iterators.tree_iter_dive_reverse(hidx.index, rootsB)
            ):
                if a.index == b.index:
                    continue
                if ordering[a.index] > ordering[b.index]:
                    continue
                # if b.name not in hidx.db:
                #     hidx.hentry_remove(b.index)
                #     continue
                gb = hidx.subgraphs[b.index]
                if not gb:
                    continue
                gb = graphs.subgraph_to_structure(gb, topo)
                # if (hidx.db[a.key]['smirks'] == hidx.db[b.key]['smirks']):
                work[b.index] = (
                    a.index,
                    pool.apply_async(mapper_match, (ga, gb)),
                )
            for idx, (idxa, matched) in work.items():
                if idx in hidx.index.nodes:
                    # print("Consider is", idx, hidx.smarts[idx], gcd.smarts_encode(hidx.subgraphs[idx]), "subset of", idxa, hidx.smarts[idxa], gcd.smarts_encode(hidx.subgraphs[idxa]))
                    # print("###")
                    # graphs.subgraph_print(hidx.subgraphs[idx])
                    # print("###")
                    # graphs.subgraph_print(hidx.subgraphs[idxa])
                    # if idx == 511 and idxa == 427:
                    #     breakpoint()
                    if matched.get():
                        # if b.key == ('n3', '^', 'n3'):
                        #     breakpoint()
                        # print("Yes, Removing due to duplicate")
                        # print("Removing nonlocal occlusion due to duplicate:", b, b.key, "kept", a.key)
                        hidx.index.node_remove(idx)
                        hidx.subgraphs.pop(idx)
                        hidx.smarts.pop(idx)
                        removed.add(idx)
            pool.terminate()

    return hidx


def structure_hierarchy_add_hierarchy(
    sA: hierarchies.structure_hierarchy,
    rootA: trees.tree_node,
    sB: hierarchies.structure_hierarchy,
    rootB: trees.tree_node,
    index=None,
):
    node = trees.tree_node(0, "", "", rootB.name)
    # hent.key = rootB.key
    node = sA.index.node_add(rootA.index, node, index=index)
    sA.subgraphs[node.index] = graphs.subgraph_copy(sB.subgraphs[rootB.index])
    sA.smarts[node.index] = str(sB.smarts[rootB.index])

    up = node.index
    mapping = {rootB.index: up}
    for eb in tree_iterators.tree_iter_breadth_first(sB.index, rootB):
        ei = eb.index
        up = mapping[sB.index.above[ei]]

        node = trees.tree_index_node_add(
            sA.index, up, trees.tree_node(None, "", "", eb.name)
        )
        ni = node.index
        sA.subgraphs[node.index] = graphs.subgraph_copy(sB.subgraphs[eb.index])
        # sA.subgraphs[ni] = sB.subgraphs[ei]
        sA.smarts[ni] = str(sB.smarts[ei])
        mapping[eb.index] = ni

    return sA
