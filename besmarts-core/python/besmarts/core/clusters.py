"""
besmarts.core.clusters

Associates a SMARTS hierarchy to a dataset group (of assignments)
"""

import os
import sys
import pickle
import datetime
import collections
import multiprocessing.pool
import threading
import time
from typing import Dict, Sequence, Tuple, List
import heapq
import gc
from pprint import pprint

from besmarts.core import (
    codecs,
    configs,
    mapper,
    graphs,
    topology,
    hierarchies,
    assignments,
    optimization,
    trees,
    tree_iterators,
    splits,
    compute,
    arrays,
)
from besmarts.cluster import cluster_assignment


class smarts_clustering:
    __slots__ = "hierarchy", "group", "mappings", "group_prefix_str"

    def __init__(self, structure_hierarchy, assign_group, mappings):
        # this is the SMARTS hierarchy
        self.hierarchy: hierarchies.structure_hierarchy = structure_hierarchy

        # this will be the labeling of the structures
        self.group: assignments.smiles_assignment_group = assign_group

        # this is the structures grouped by the labels
        self.mappings: assignments.assignment_mapping = mappings

        # when new nodes are created, use this prefix
        self.group_prefix_str = "p"


class clustering_objective:
    def split(self, A, B) -> float:
        raise NotImplementedError()

    def merge(self, A, B) -> float:
        raise NotImplementedError()

    def single(self, A) -> float:
        raise NotImplementedError()

    def report(self, A) -> str:
        raise NotImplementedError()

    def is_discrete(self) -> bool:
        raise NotImplementedError()

    def sum(self) -> bool:
        raise NotImplementedError()


def objective_total(hidx, groups, objective):
    X = 0.0
    for s in hidx.index.nodes.values():
        X += objective.single(groups[s.name])
    return X


def get_objective(cst, assn, objfn, edits, splitting=True):
    hidx = cst.hierarchy
    new_match = cst.mappings
    keep = True
    obj = 0.0
    X = 0.0
    for n in tree_iterators.tree_iter_dive(
        hidx.index, trees.tree_index_roots(hidx.index)
    ):
        m = hidx.index.above.get(n.index, None)
        if m is not None:
            m = hidx.index.nodes[m]
            if m.type != "parameter" or n.type != "parameter":
                continue
            n_match = new_match[n.name]
            m_match = new_match[m.name]
            # if not (n_match and m_match):
            #     keep = False
            #     print(f"RETURNING FALSE1 because {n_match} {m_match}")
            #     continue
            n_group = tuple(((assn[i] for i in n_match)))
            m_group = tuple(((assn[i] for i in m_match)))
            # print(f"Objective for Sj: {n.name} {m.name} ->")
            obj = objfn(n_group, m_group, overlap=edits)
            if False and splitting:
                if obj >= 0.0:
                    keep = False
                # print(f"Object increased to {obj} for {n.name} parent {m.name}")
                # continue
            X += obj
    # print("Total objective:", X)
    return keep, X



def find_successful_candidates_distributed(S, Sj, operation, edits, shm=None):

    sag = shm.sag
    cst = shm.cst
    hidx = cst.hierarchy.copy()
    labeler = shm.labeler
    gcd = shm.gcd

    objective = shm.objective

    smiles = [a.smiles for a in sag.assignments]
    topo = hidx.topology
    assn = shm.assn # get_assns(sag.assignments, topo)

    obj = edits
    X = edits

    # (S, Sj, step, _, _, _, _) = candidates[key]
    # (edits, _, p_j) = key
    param_name = "pX"
    sma = ""
    added = False
    dX = 0.0

    if operation == optimization.optimization_strategy.SPLIT:
        # param_name = "p" + str(group_number)

        # print(datetime.datetime.now(), '*** 2')
        hent = hidx.index.node_add(
            S.index,
            trees.tree_node(None, S.category, S.type, param_name),
            index=0,
        )
        # print(datetime.datetime.now(), '*** 3')
        Sj = graphs.subgraph_relabel_nodes(
            Sj, {n: i for i, n in enumerate(Sj.select, 1)}
        )
        Sj = graphs.subgraph_to_structure(Sj, topo)

        # sma = Sj_sma[cnd_i] #gcd.smarts_encode(Sj)

        hidx.subgraphs[hent.index] = Sj
        hidx.smarts[hent.index] = shm.gcd.smarts_encode(graphs.subgraph_as_structure(Sj, topo))

        # print(datetime.datetime.now(), '*** 4')
        new_assignments = labeler.assign(hidx, gcd, smiles, topo)

        # print(datetime.datetime.now(), '*** 5')
        new_match = clustering_build_assignment_mappings(hidx, new_assignments)

        cst = smarts_clustering(hidx, new_assignments, new_match)

        # print(datetime.datetime.now(), '*** 6')
        groups = clustering_build_ordinal_mappings(cst, sag, [S.name, hent.name])

        if not (groups[S.name] and groups[hent.name]):
            keep = False
            _, X = get_objective(
                cst, assn, objective.split, edits, splitting=True
            )
        else:
            obj = objective.split(groups[S.name], groups[hent.name], overlap=edits)
            # print(datetime.datetime.now(), '*** 7')
            keep, X = get_objective(
                cst, assn, objective.split, edits, splitting=True
            )
            # dX = X - X
            # if dX > 0:
            #     keep = False
            if obj >= 0.0:
                keep = False

            if not (cst.mappings[S.name] and cst.mappings[hent.name]):
                keep = False
            # print(datetime.datetime.now(), '*** 8', keep, X, obj)

            # keep the splits that match the most (general)
            # match_len = len(cst.mappings[S.name])

            # keep the splits that match the least (specific)
            # this is better since it leaves more in S, so more can be split
            # at a time
        match_len = len(cst.mappings[hent.name])

        return keep, X, obj, match_len

        # print(datetime.datetime.now(), '*** 8')
    elif operation == optimization.optimization_strategy.MERGE:

        hent = Sj
        # print(datetime.datetime.now(), '*** 8')
        groups = clustering_build_ordinal_mappings(cst, sag, [S.name, hent.name])
        obj = 0.0
        if S.name in groups and hent.name in groups:
            obj = objective.merge(groups[S.name], groups[hent.name], overlap=edits)
        trees.tree_index_node_remove(hidx.index, Sj.index)
        # print(datetime.datetime.now(), '*** 9')
        new_assignments = labeler.assign(hidx, gcd, smiles, topo)
        # print(datetime.datetime.now(), '*** 10')
        new_match = clustering_build_assignment_mappings(hidx, new_assignments)
        cst = smarts_clustering(hidx, new_assignments, new_match)
        # print(datetime.datetime.now(), '*** 11')
        _, X = get_objective(cst, assn, objective.split, edits, splitting=False)
        # Sj = hidx.subgraphs[Sj.index]
        # cst.hierarchy.subgraphs.pop(hent.index)
        # cst.hierarchy.smarts.pop(hent.index)
        # sma = Sj_sma[cnd_i-1]
        keep = (
            obj < 0
            or hent.name in groups and len(groups[hent.name]) == 0
        )
        if S.name in cst.mappings:
            match_len = len(cst.mappings[S.name])
        else:
            keep = False

        return keep, X, obj, match_len
    else:
        return False, 0.0, 0.0, 0

def perform_operations(
        hidx: hierarchies.structure_hierarchy,
        candidates,
        keys,
        group_number,
        Sj_sma,
        strategy,
        prefix="p"
    ):

    topo = hidx.topology

    nodes = {}
    for cnd_i, key in keys.items():
        (S, Sj, step, _, _, _, _) = candidates[key]
        (edits, _, p_j, oper) = key
        param_name = "p."
        sma = ""
        added = False

        if oper == strategy.SPLIT:
            # param_name = prefix + str(group_number)
            i = max(hidx.index.nodes) + 1
            name = f"{prefix}{i}"
            existing_names = [x.name for x in hidx.index.nodes.values()]
            while name in existing_names:
                i += 1
                name = f"{prefix}{i}"

            # print(datetime.datetime.now(), '*** 2')
            hent = hidx.index.node_add(
                S.index,
                trees.tree_node(None, S.category, S.type, name),
                index=0,
            )

            # print(datetime.datetime.now(), '*** 3')
            Sj = graphs.subgraph_relabel_nodes(
                Sj, {n: i for i, n in enumerate(Sj.select, 1)}
            )
            Sj = graphs.subgraph_to_structure(Sj, topo)

            sma = Sj_sma[cnd_i-1]  # gcd.smarts_encode(Sj)

            hidx.subgraphs[hent.index] = Sj
            hidx.smarts[hent.index] = sma
            nodes[cnd_i] = hent

            group_number += 1

            # print(datetime.datetime.now(), '*** 4')

            #####
        elif oper == strategy.MERGE:
            # if (S.name not in groups) or (Sj.name not in groups):
            #     continue

            hent = Sj
            nodes[cnd_i] = hent
            # obj += objective.merge(groups[S.name], groups[hent.name], overlap=edits)
            trees.tree_index_node_remove(hidx.index, Sj.index)

            hidx.subgraphs.pop(hent.index)
            hidx.smarts.pop(hent.index)

    return hidx, nodes


def check_lbls_data_selections_equal(
    lbls: assignments.smiles_assignment_group,
    data: assignments.smiles_assignment_group,
):
    # check if there as many labeled ICs to data points.
    warning_max = 10
    warnings = 0
    for idx, (lbl_assn, data_assn) in enumerate(
        zip(lbls.assignments, data.assignments)
    ):
        assert lbl_assn.smiles == data_assn.smiles
        for lbl_ic in lbl_assn.selections:
            if lbl_ic not in data_assn.selections:
                if warnings < warning_max:
                    print(
                        f"WARNING: mol {idx} atoms {lbl_ic} does not have data! This will likely fail."
                    )
                    print(f"WARNING:     SMILES: {lbl_assn.smiles}")
                warnings += 1
    if warnings > warning_max:
        print(
            f"WARNING: suppressed {warnings - warning_max} additional warnings"
        )


def smarts_clustering_optimize(
    gcd: codecs.graph_codec,
    labeler: assignments.smarts_hierarchy_assignment,
    sag: assignments.smiles_assignment_group,
    objective: clustering_objective,
    strategy: optimization.optimization_strategy,
    initial_conditions: smarts_clustering,
) -> smarts_clustering:

    # gc.disable()
    started = datetime.datetime.now()

    smiles = [a.smiles for a in sag.assignments]

    topo = sag.topology

    group_prefix_str = initial_conditions.group_prefix_str

    hidx = initial_conditions.hierarchy.copy()

    groups = clustering_build_ordinal_mappings(initial_conditions, sag)
    # print(groups)

    # match = clustering_build_assignment_mappings(initial_conditions, initial_conditions.group)
    match = initial_conditions.mappings
    assn = get_assns(sag.assignments, topo)

    icd = codecs.intvec_codec(gcd.primitive_codecs, gcd.atom_primitives, gcd.bond_primitives)

    wq = compute.workqueue_local('0.0.0.0', configs.workqueue_port)

    check_lbls_data_selections_equal(initial_conditions.group, sag)


    """
    nomenclature

    match is node_name:data_idx (mol, ic)
    assn is data_idx:data
    group is node_name:data
    mapping node_name:data_idx
    """

    group_number = max(
        [0] + [*map(int, [x for x in 
            ["".join(a for a in x.name[1:] if a.isdigit())
            for x in hidx.index.nodes.values()] if x])]
    )

    group_number += 1

    # gc.collect()
    if len(sag.assignments) > 100000 and (configs.remote_compute_enable or configs.processors > 1):
        batch_size = 10000
        print(f"{datetime.datetime.now()} Large number of graphs detected... using a workspace")
        # wq = compute.workqueue_local('', configs.workqueue_port)
        ws = compute.workqueue_new_workspace(wq, address=('127.0.0.1', 0), shm={"gcd": gcd})

        work = compute.workspace_submit_and_flush(
            ws,
            codecs.smiles_decode_list_distributed,
            {i: ((list(e),), {}) for i, e in enumerate(arrays.batched((s.smiles for s in sag.assignments), batch_size))},
            chunksize=10
        )
        
        G0 = {}
        n_ics = 0
        for ii in sorted(work):
            i = 0
            for i, ig in enumerate(work[ii], ii*batch_size):
                G0[i] = ig
                sels = sag.assignments[i].selections
                n_ics += len(sels)
            print(f"\r{datetime.datetime.now()} graphs= {i+1:8d}/{len(sag.assignments)} subgraphs= {n_ics:8d}", end="")
        print()
            
            
        compute.workqueue_remove_workspace(wq, ws)
        ws.close()
        #threading.Thread(target=ws.close).start()
        ws = None

        # wq.close()
        # wq = None  
        # gc.collect()

    else:

        G0 = {i: icd.graph_encode(gcd.smiles_decode(a.smiles)) for i, a in enumerate(sag.assignments)}
        n_ics = sum((len(s.selections) for s in sag.assignments))

    N = len(hidx.index.nodes)
    try:
        N = len(set(assn.values()))
    except Exception:
        pass

    repeat = set()
    visited = set()
    iteration = 1
    N_str = "{:" + str(len(str(n_ics))) + "d}"

    success = False

    roots = trees.tree_index_roots(hidx.index)

    print(f"{datetime.datetime.now()} Labeling subgraphs")
    assignments = labeler.assign(hidx, gcd, smiles, topo)
    cst = smarts_clustering(
        hidx,
        assignments,
        clustering_build_assignment_mappings(hidx, assignments),
    )
    cst.group_prefix_str = group_prefix_str
    print(f"{datetime.datetime.now()} Checking consistency...")
    check_lbls_data_selections_equal(cst.group, sag)
    _, X0 = get_objective(
        cst, assn, objective.split, strategy.overlaps[0], splitting=True
    )

    if not strategy.steps:
        print("Optimization strategy is building steps...")
        strategy.build_steps()
    print(f"{datetime.datetime.now()} The optimization strategy has the following iterations:")
    for ma_i, macro in enumerate(strategy.steps, 1):
        cur = "  "
        if ma_i == strategy.cursor + 1:
            cur = "->"
        for mi_i, micro in enumerate(macro.steps):
            s = micro.pcp.splitter
            a = micro.overlap
            b0 = s.bit_search_min
            b1 = s.bit_search_limit
            d0 = s.branch_depth_min
            d1 = s.branch_depth_limit
            n0 = s.branch_min
            n1 = s.branch_limit
            print(
                f"{cur} {ma_i:3d}. op={micro.operation:2d} a={a} b={b0}->{b1} d={d0}->{d1} n={n0}->{n1}"
            )

    union_cache = {}

    step_tracker = strategy.step_tracker

    while True:
        if success:
            print("Restarting optimization search")
            strategy = optimization.optimization_strategy_restart(strategy)
            success = False

        elif optimization.optimization_strategy_is_done(strategy):
            print("Nothing found. Done.")
            break

        groups = clustering_build_ordinal_mappings(cst, sag)

        roots = trees.tree_index_roots(cst.hierarchy.index)
        nodes = [
            x
            for x in strategy.tree_iterator(cst.hierarchy.index, roots)
            if cst.hierarchy.smarts.get(x.index)
            and cst.hierarchy.subgraphs.get(x.index) is not None
            and (strategy.cursor >= step_tracker.get(x.name, -1))
        ]
        for n in nodes:
            tkey = n.name
            if tkey not in step_tracker:
                step_tracker[tkey] = 0

        print(f"Targets for this macro step {strategy.cursor+1}:")
        for nidx, n in enumerate(nodes, 1):
            print(nidx, n.name)
        print(f"N Targets: {len(nodes)}")

        print(f"Step tracker for current macro step {strategy.cursor+1}")
        for n, v in step_tracker.items():
            print(n, v + 1)

        macro: optimization.optimization_iteration = strategy.macro_iteration(
            nodes
        )

        candidates = {}
        pq = []
        n_added = 0
        n_macro = len(strategy.steps)

        t = datetime.datetime.now()
        config = strategy.bounds
        spg = "Y" if config.splitter.split_general else "N"
        sps = "Y" if config.splitter.split_specific else "N"
        print(
            f"\n\n*******************\n {t}"
            f" iteration={iteration:4d}"
            f" macro={strategy.cursor:3d}/{n_macro}"
            f" X={X0:9.5g}"
            f" params=({len(cst.mappings)}|{N})"
            f" G={spg}"
            f" S={sps}"
            f" bits={config.splitter.bit_search_min}->{config.splitter.bit_search_limit}"
            f" depth={config.splitter.branch_depth_min}->{config.splitter.branch_depth_limit}"
            f" branch={config.splitter.branch_min}->{config.splitter.branch_limit}"
        )
        print(f"*******************")
        print()
        print("Tree:")
        for ei, e in enumerate(
            tree_iterators.tree_iter_dive(
                cst.hierarchy.index, trees.tree_index_roots(cst.hierarchy.index)
            )
        ):
            s = trees.tree_index_node_depth(cst.hierarchy.index, e)
            obj_repo = ""
            # if groups[e.name]:
            obj_repo = objective.report(groups[e.name])
            print(
                f"** {s:2d} {ei:3d} {e.name:4s}",
                obj_repo,
                cst.hierarchy.smarts.get(e.index),
            )
        print("=====\n")

        print(f"{datetime.datetime.now()} Saving checkpoint to chk.cst.p")
        pickle.dump([sag, cst, strategy], open("chk.cst.p", "wb"))

        step = None

        while not optimization.optimization_iteration_is_done(macro):
            # t = datetime.datetime.now()
            # print(f"{t} Initializing new loop on macro {strategy.cursor}")
            step: optimization.optimization_step = (
                optimization.optimization_iteration_next(macro)
            )

            n_micro = len(macro.steps)
            config: configs.smarts_perception_config = step.pcp
            S = step.cluster

            # step_tracker[S.name] = strategy.cursor

            if type(cst.hierarchy.subgraphs[S.index]) is str:
                step_tracker[S.name] = strategy.cursor
                continue

            S0 = graphs.subgraph_to_structure(
                cst.hierarchy.subgraphs[S.index], topo
            )

            cfg = config.extender.copy()

            S0_depth = graphs.structure_max_depth(S0)
            d = max(S0_depth, config.splitter.branch_depth_limit)
            # print("Depth is", d)
            cfg.depth_max = d
            cfg.depth_min = S0_depth

            t = datetime.datetime.now()
            print(
                f"{t} Collecting SMARTS for {S.name} N={len(cst.mappings[S.name])}/{n_ics} and setting to depth={S0_depth}"
            )

            # find only unique graphs!
            aa = cst.mappings[S.name]
            selected_graphs = set((x[0] for x in aa))
            G = {k:v for k, v in G0.items() if k in selected_graphs}
            del selected_graphs

            assn_s = {i: assn[i] for i in cst.mappings[S.name]}

            iteration += 1

            t = datetime.datetime.now()
            print(
                f" =="
                f" iteration={iteration:4d}"
                f" macro={strategy.cursor:3d}/{n_macro}"
                f" micro={macro.cursor:3d}/{n_micro}"
                # f" overlap={strategy.overlaps:3d}"
                f" operation={step.operation}"
                f" params=({len(cst.mappings)}|{N})"
                f" cluster={S.name:4s}"
                f" N= " + N_str.format(len(aa)) + ""
                f" overlap={step.overlap}"
                f" bits={config.splitter.bit_search_min}->{config.splitter.bit_search_limit}"
                f" depth={config.splitter.branch_depth_min}->{config.splitter.branch_depth_limit}"
                f" branch={config.splitter.branch_min}->{config.splitter.branch_limit}"
            )
            print()

            if step.operation == strategy.SPLIT:
                new_pq = []
                new_candidates = {}
                new_candidates_direct = {}
                direct_success = False

                print(f"Attempting to split {S.name}:")
                print("S0:", gcd.smarts_encode(S0))



                if not assn_s:
                    print("No matches.")
                    step_tracker[S.name] = strategy.cursor
                    continue

                print(f"Matched N={len(assn_s)}")
                seen = set()
                # depth = graphs.structure_max_depth(S0)
                extend_config = config.extender.copy()
                extend_config.depth_max = config.splitter.branch_depth_limit
                extend_config.depth_min = config.splitter.branch_depth_min
                for seen_i, (i, x) in enumerate(assn_s.items(), 1):
                    g = graphs.graph_to_structure(icd.graph_decode(G[i[0]]), i[1], topo) 
                    graphs.structure_extend(extend_config, [g])
                    g = graphs.structure_remove_unselected(g)
                    # if g not in seen:
                    # seen_g.add(g)
                    if seen_i < (configs.match_print_limit or len(assn_s)+1):
                        print(
                            f"{seen_i:06d} {str(i):24s}",
                            objective.report([x]),
                            gcd.smarts_encode(g),
                        )
                        seen.add(g)
                print()
                if len(seen) < 2 and len(assn_s) < (configs.match_print_limit or len(assn_s)+1):
                    print(f"Skipping {S.name} since all graphs are the same")
                    step_tracker[S.name] = strategy.cursor
                    continue

                if len(set(map(tuple, groups[S.name]))) == 1:
                    print(f"Skipping {S.name} since all data are the same")
                    step_tracker[S.name] = strategy.cursor
                    continue

                if len(seen) < 2 and len(assn_s) < (configs.match_print_limit or len(assn_s)+1):
                    print(f"Skipping {S.name} since all graphs are the same")
                    step_tracker[S.name] = strategy.cursor
                    continue

                if objective.single(assn_s.values()) == 0.0:
                    print(f"Skipping {S.name} due to no objective")
                    step_tracker[S.name] = strategy.cursor
                    continue

                if graphs.structure_max_depth(S0) > config.splitter.branch_depth_min:
                    print("This parameter exceeds current depth. Skipping")
                    step_tracker[S.name] = strategy.cursor
                    continue

                if step.direct_enable and config.splitter.bit_search_limit < len(assn_s):
                    assn_i = []
                    if len(assn_s) < step.direct_limit:
                        # or form matches based on unique smarts
                        a = []
                        extend_config = config.extender.copy()
                        extend_config.depth_max = config.splitter.branch_depth_limit
                        extend_config.depth_min = config.splitter.branch_depth_min
                        for seen_i, (i, x) in enumerate(assn_s.items(), 1):
                            g = graphs.graph_to_structure(icd.graph_decode(G[i[0]]), i[1], topo) 
                            graphs.structure_extend(extend_config, [g])
                            g = graphs.structure_remove_unselected(g)
                            a.append(g)
                        # if len(a) < step.direct_limit:
                        #     lbls = dict()
                            # for (i, idx), sg in a.items():
                            #     lbl_i = lbls.get(a, len(lbls))
                            #     lbls[a] = lbl_i
                            #     assn_i.append(lbl_i)

                        if objective.is_discrete():
                            assn_i.extend(groups[S.name])
                        else:
                            assn_i = list(range(len(a)))

                        pcp = step.pcp.copy()
                        pcp.extender = cfg
                        print("Direct splitting....")
                        # pcp.extender.depth_max = pcp.splitter
                        # pcp.extender.depth_min = strategy.bounds.extender.depth_min
                        ret = splits.split_all_partitions(
                            topo,
                            pcp,
                            a,
                            assn_i,
                            gcd=gcd,
                            maxmoves=0,
                        )

                        for p_j, (Sj, Sj0, matches, unmatches) in enumerate(ret.value, 0):
                            print(f"Found {p_j+1} {gcd.smarts_encode(Sj)}")

                            edits = 0
                            matches = [
                                y
                                for x, y in enumerate(aa)
                                if x in matches
                            ]
                            unmatches = [
                                y
                                for x, y in enumerate(aa)
                                if x in unmatches
                            ]
                            matched_assn = tuple((assn[i] for i in matches))
                            unmatch_assn = tuple((assn[i] for i in unmatches))

                            new_candidates_direct[(step.overlap[0], None, p_j)] = (
                                S,
                                Sj,
                                step,
                                matches,
                                unmatches,
                                matched_assn,
                                unmatch_assn,
                            )
                        if len(new_candidates_direct):
                            direct_success = True
                        else:
                            print("Direct found nothing")


                if step.iterative_enable and not direct_success:
                    # Q = mapper.union_list_parallel(
                    #     list(a.values()),
                    #     reference=S0,
                    #     max_depth=graphs.structure_max_depth(S0),
                    #     icd=icd if len(a) > 100000 else None
                    # )

                    (Q, new_candidates) = union_cache.get((S.index, S.name, strategy.cursor), (None, None))
                    if Q is None:
                        if len(aa) < 1000:
                            Q = mapper.union_list_parallel(
                                G, aa, topo,
                                # list(a.values()),
                                reference=S0,
                                max_depth=graphs.structure_max_depth(S0),
                                icd=icd
                                # icd=icd if len(a) > 100000 else None
                            )
                        else:
                            Q = mapper.union_list_distributed(
                                G, aa, topo, wq,
                                # list(a.values()),
                                reference=S0,
                                max_depth=graphs.structure_max_depth(S0),
                                icd=icd
                                # icd=icd if len(a) > 100000 else None
                            )

                        union_cache[(S.index, S.name, strategy.cursor)] = (Q, new_candidates)
                    else:
                        print(
                            f"{datetime.datetime.now()} Candidates retreived from cache for node {S.index}:{S.name}"
                        )

                    t = datetime.datetime.now()
                    print(f"{t} Union is {gcd.smarts_encode(Q)}")

                    if new_candidates is None:
                        return_matches = config.splitter.return_matches
                        config.splitter.return_matches = True

                        print("PRIMITIVES IS", config.splitter.primitives)
                        ret = splits.split_structures_distributed(
                            config.splitter,
                            S0,
                            G,
                            aa,
                            wq,
                            icd,
                            Q=Q,
                        )
                        config.splitter.return_matches = return_matches

                        backmap = {i: j for i, j in enumerate(cst.mappings[S.name])}
                        print(
                            f"{datetime.datetime.now()} Collecting new candidates"
                        )

                        new_candidates = clustering_collect_split_candidates_serial(
                            S, ret, step, strategy.SPLIT
                        )
                        union_cache[(S.index, S.name, strategy.cursor)] = (Q, new_candidates)


                p_j_max = -1
                if candidates:
                    p_j_max = max(x[2] for x in candidates) + 1
                for k, v in new_candidates.items():
                    k = (k[0], k[1], k[2]+p_j_max, step.operation)
                    candidates[k] = v
                new_candidates = None

                p_j_max = -1
                if candidates:
                    p_j_max = max(x[2] for x in candidates) + 1
                for k, v in new_candidates_direct.items():
                    k = (k[0], k[1], k[2]+p_j_max, step.operation)
                    candidates[k] = v
                new_candidates_direct = None

            elif step.operation == strategy.MERGE:

                for p_j, jidx in enumerate(cst.hierarchy.index.below[S.index]):
                    J = cst.hierarchy.index.nodes[jidx]

                    for overlap in step.overlap:
                        key = (overlap, macro.cursor, p_j, strategy.MERGE)
                        cnd = (S, J, step, None, None, None, None)
                        candidates[key] = cnd

        if step is None:
            print(
                f"{datetime.datetime.now()} Warning, this macro step had no micro steps"
            )
            continue

        print(f"{datetime.datetime.now()} Scanning done.")

        print(datetime.datetime.now())
        print(f"\n\nGenerating SMARTS on {len(candidates)}")
        Sj_sma = []

        with multiprocessing.pool.Pool(configs.processors) as pool:
            Sj_lst = []
            for (_,_,_, oper), x in candidates.items():
                if oper == strategy.SPLIT:
                    Sj_lst.append(graphs.subgraph_as_structure(x[1], topo))
                elif oper == strategy.MERGE:
                    Sj_lst.append(graphs.subgraph_as_structure(cst.hierarchy.subgraphs[x[1].index], topo))
                elif oper == strategy.MODIFY:
                    Sj_lst.append(x[0])
            Sj_sma = pool.map_async(gcd.smarts_encode, Sj_lst).get()
            del Sj_lst

        print(f"{datetime.datetime.now()} Labeling")
        cur_assignments = labeler.assign(cst.hierarchy, gcd, smiles, topo)
        print(f"{datetime.datetime.now()} Rebuilding assignments")
        cur_mappings = clustering_build_assignment_mappings(
            cst.hierarchy, cur_assignments
        )
        cur_cst = smarts_clustering(
            cst.hierarchy.copy(), cur_assignments, cur_mappings
        )
        cur_cst.group_prefix_str = group_prefix_str
        print(f"{datetime.datetime.now()} Rebuilding mappings")
        groups = clustering_build_ordinal_mappings(cur_cst, sag)
        check_lbls_data_selections_equal(cst.group, sag)

        cnd_n = len(candidates)

        t = datetime.datetime.now()

        print("Tree:")
        for ei, e in enumerate(
            tree_iterators.tree_iter_dive(
                cur_cst.hierarchy.index,
                trees.tree_index_roots(cur_cst.hierarchy.index),
            )
        ):
            s = trees.tree_index_node_depth(cur_cst.hierarchy.index, e)
            obj_repo = ""
            # if groups[e.name]:
            obj_repo = objective.report(groups[e.name])
            print(
                f"** {s:2d} {ei:3d} {e.name:4s}",
                obj_repo,
                cur_cst.hierarchy.smarts.get(e.index),
            )
        print("=====\n")

        visited.clear()
        repeat.clear()


        pq_idx = 0
        procs = (
            os.cpu_count() if configs.processors is None else configs.processors
        )

        n_keep = None

        print(f"Scanning {len(candidates)} candidates for operation={step.operation}")

        macroamt = strategy.macro_accept_max_total
        macroampc = strategy.macro_accept_max_per_cluster
        microamt = strategy.micro_accept_max_total
        microampc = strategy.micro_accept_max_per_cluster

        cnd_n = len(candidates)
        n_added = 0
        added = True
        kept = set()
        macro_count = collections.Counter()
        ignore = set()
        reuse = {}
        # wq = compute.workqueue_local("", configs.workqueue_port)
        # print(f"{datetime.datetime.now()} workqueue started on {wq.mgr.address}")
        n_nano = 0
        while added:

            case1 = macroamt == 0 or n_added < macroamt
            case2 = macroampc == 0 or all([x < macroampc for x in macro_count.values()])
            if not (case1 and case2):
                break

            n_nano += 1

            added = False
            best = {}

            cout = {}
            cout_sorted_keys = []

            shm = compute.shm_local(1, data={
                "cst": cur_cst,
                "sag": sag,
                "gcd": gcd,
                "labeler": labeler,
                "objective": objective,
                "assn": assn
            })

            iterable = {
                i: ((S, Sj, oper, edits), {})
                for i, (
                    (edits, _, p_j, oper),
                    (S, Sj, step, _, _, _, _),
                ) in enumerate(candidates.items(), 1)
            }

            chunksize = 10

            if n_ics > 100000000:
                procs = max(1, procs // 10)
            elif n_ics > 50000000:
                procs = max(1, procs // 5)
            elif n_ics > 10000000:
                procs = max(1, procs // 3)
            elif n_ics > 5000000:
                procs = max(1, procs // 2)
            if n_ics > len(candidates)*10:
                shm.procs_per_task = 0
                chunksize = 1

            addr = ("", 0)
            if len(iterable)*(shm.procs_per_task or procs) <= procs:
                addr = ('127.0.0.1', 0)
                procs = len(iterable)

            cnd_keys = {i: k for i, k in enumerate(candidates, 1)}

            for k in kept:
                if k in iterable:
                    iterable.pop(k)

            for k in ignore:
                if k in iterable:
                    iterable.pop(k)

            for k in reuse:
                if k in iterable:
                    iterable.pop(k)

            for k, v in list(candidates.items()):
                S = v[0]
                if S.name in cur_cst.mappings:
                    if objective.single([
                        assn[i] for i in cur_cst.mappings[S.name]
                    ]) == 0.0:
                        if k in iterable:
                            iterable.pop(k)
                elif k in iterable:
                    iterable.pop(k)

            # if step.operation != strategy.MERGE and (macroamt + macroampc + microamt + microampc == 0):
            #     # use X0 so later dX will be 0 and kept
            #     # if we do this for merges, every merge will be taken..
            #     work = {i: (1, X0, 0.0, 1) for i in iterable}

            # else:
            ws = compute.workqueue_new_workspace(
                wq,
                address=addr,
                nproc=procs,
                shm=shm
            )
            # this modifies the cst, relabels and computes objective
            work = compute.workspace_submit_and_flush(
                ws,
                find_successful_candidates_distributed,
                iterable,
                chunksize,
                0.0,
                len(iterable),
                verbose=True
            )
            compute.workqueue_remove_workspace(wq, ws)
            ws.close()
            ws = None

            print(f"The unfiltered results of the candidate scan N={len(work)} total={len(iterable)}:")

            max_line = 0

            best_reuse = None
            if reuse:
                # just append the best to work and let the loop figure it out
                best_reuse = sorted(reuse.items(), key=lambda y: (-y[1][0], y[1][1], y[1][2], y[1][3]))[0]
                work[best_reuse[0]] = best_reuse[1]
            for j, cnd_i in enumerate(sorted(work), 1):
                (keep, X, obj, match_len) = work[cnd_i]
                # cnd_i, key, unit = unit
                oper = cnd_keys[cnd_i][3]
                (S, Sj, step, _, _, _, _) = candidates[cnd_keys[cnd_i]]

                if oper == strategy.SPLIT:
                    visited.add((S.name, oper))
                elif oper == strategy.MERGE:
                    visited.add((Sj.name, oper))
                elif oper == strategy.MODIFY:
                    visited.add((S.name, oper))

                dX = X - X0
                reused_line = ""
                if best_reuse is not None and cnd_i == best_reuse[0]:
                    reused_line="*"
                C = "Y" if keep else "N"
                cout_line = (
                    f"Cnd. {cnd_i:4d}/{len(work)}"
                    f" OP={oper}"
                    f" {S.name:6s} {reused_line}"
                    f" X= {X:10.5f}"
                    f" dX= {dX:10.5f} N= {match_len:6d} C= {C} {Sj_sma[cnd_i-1]}"
                )
                max_line = max(len(cout_line), max_line)
                # print(datetime.datetime.now())
                print(cout_line, end=" " * (max_line - len(cout_line)) + '\n')
                sys.stdout.flush()

                if match_len == 0:
                    if oper == strategy.SPLIT:
                        keep = False
                        ignore.add(cnd_i)
                        continue

                if not keep:
                    ignore.add(cnd_i)
                    continue


                if cnd_i in kept:
                    ignore.add(cnd_i)
                    continue

                # We prefer to add in this order
                cout_key = None

                # print sorted at the end but only for new
                # this is to speed things up
                cout_key = (-int(keep), X, match_len, cnd_i, S.name)
                cout[cout_key] = cout_line

                # use these below to determine the best ones to keep
                heapq.heappush(cout_sorted_keys, cout_key)

            print("\r" + " " * max_line)

            # if best:
            #     for name, v in best.items():
            #         kept.add(v[0])
            #     added = True
            # else:
            #     break

            # print sorted at the end
            print(f"Nanostep {n_nano}: The filtered results of the candidate scan N={len(cout)} total={len(iterable)}:")
            if len(cout) == 0:
                continue
            ck_i = 1

            cnd_keep = []
            # best_params = [x[0] for x in best.values()]
            macroamt = strategy.macro_accept_max_total
            macroampc = strategy.macro_accept_max_per_cluster
            microamt = strategy.micro_accept_max_total
            microampc = strategy.micro_accept_max_per_cluster
            micro_added = 0
            micro_count = collections.Counter()

            while len(cout_sorted_keys):
                ck = heapq.heappop(cout_sorted_keys)

                keeping = "  "

                dX = ck[1] - X0
                case0 = not (strategy.filter_above is not None and strategy.filter_above <= dX)

                if case0:
                    ignore.add(ck[0])

                case1 = macroamt == 0 or n_added < macroamt
                case2 = microamt == 0 or micro_added < microamt
                if case0 and case1 and case2:
                    sname = ck[4]
                    case3 = macroampc == 0 or macro_count[sname] < macroampc
                    case4 = microampc == 0 or micro_count[sname] < microampc
                    if case3 and case4:
                        cnd_keep.append(ck)
                        micro_count[sname] += 1
                        macro_count[sname] += 1
                        micro_added += 1
                        n_added += 1
                # if ck[3] in best_params:
                        keeping = "->"
                        kept.add(ck[0])
                    else:
                        # if we only accept 1 but this is still a good param
                        # put it on the repeat list or else it will get skipped
                        # as we increment nonrepeating params to the next step
                        repeat.add(ck[4])
                else:
                    repeat.add(ck[4])
                print(f"{keeping} {ck_i:4d}", cout[ck])
                ck_i += 1
            ck = None

            # keys = {x[0]: cnd_keys[x[0]] for x in best.values()}
            keys = {x[3]: cnd_keys[x[3]] for x in cnd_keep}

            group_number = max(cur_cst.hierarchy.index.nodes)+1
            print(f"Performing {len(keys)} operations")
            hidx, nodes = perform_operations(
                cur_cst.hierarchy,
                candidates,
                keys,
                group_number,
                Sj_sma,
                strategy,
                prefix=cur_cst.group_prefix_str
            )

            print(f"There are {len(nodes)} nodes returned")

            print("Operations per parameter for this micro:")
            print(micro_count)
            print(f"Micro total: {sum(micro_count.values())} should be {micro_added}")

            print("Operations per parameter for this macro:")
            print(macro_count)
            print(f"Macro total: {sum(macro_count.values())} should be {n_added}")

            if len(nodes) == 0:
                added = False
                continue

            new_assignments = labeler.assign(hidx, gcd, smiles, topo)
            # print(datetime.datetime.now(), '*** 5')
            new_match = clustering_build_assignment_mappings(hidx, new_assignments)

            cst = smarts_clustering(hidx, new_assignments, new_match)

            groups = clustering_build_ordinal_mappings(cst, sag)

            if strategy.prune_empty:
                prune_count = 0
                pruned = []
                new_lbls = [n.name for n in nodes.values()]
                for n  in list(nodes):
                    # only consider new nodes since we have other state like
                    # the tracker that needs to be managed and would be more
                    # complicated to handle
                    n = nodes[n]
                    lbl = n.name

                    # also consider the case where we occlude the old pattern
                    m = cst.hierarchy.index.above.get(n.index, None)
                    if m is None:
                        continue
                    m = cst.hierarchy.index.nodes[m]
                    
                    if (len(groups.get(m.name, [])) == 0 or len(groups[lbl]) == 0) and lbl in new_lbls:
                        n_list = trees.tree_index_node_remove_by_name(cst.hierarchy.index, lbl)
                        for n in n_list:
                            cst.hierarchy.subgraphs.pop(n.index)
                            cst.hierarchy.smarts.pop(n.index)
                            prune_count += 1
                            pruned.append(n.index)
                            micro_count[m.name] -= 1
                            micro_count[m.name] = max(micro_count[m.name], 0)
                            macro_count[m.name] -= 1
                            macro_count[m.name] = max(macro_count[m.name], 0)
                            print(f"Pruned {n.name} parent len {len(groups.get(m.name, []))} sj len {len(groups[lbl])}")
                        # del groups[lbl]
                        # if lbl in cst.mappings:
                        #     del cst.mappings[lbl]
                        
                for cnd_i in list(nodes):
                    hent = nodes[cnd_i]
                    if hent.index in pruned:
                        del nodes[cnd_i]
                        ignore.add(cnd_i)
                        del keys[cnd_i]
                if prune_count:
                #     groups = clustering_build_ordinal_mappings(cst, sag)
                    new_assignments = labeler.assign(hidx, gcd, smiles, topo)
                    # print(datetime.datetime.now(), '*** 5')
                    new_match = clustering_build_assignment_mappings(hidx, new_assignments)

                    cst = smarts_clustering(hidx, new_assignments, new_match)

                    groups = clustering_build_ordinal_mappings(cst, sag)

                print(f"Pruned {prune_count} empty nodes; candidates now {len(keys)}/{len(cnd_keep)}")
                print(pruned)
                del pruned
                del prune_count
                del new_lbls

            


            if not nodes:
                success = False
                added = False
                continue
            success = False
            added = False
            
            #group_number += len(keys)
            _, X = get_objective(cst, assn, objective.split, step.overlap[0], splitting=False)
            dX = X-X0

            recalc = False
            for cnd_i, hent in nodes.items():
                key = keys[cnd_i]

                oper = key[3]
                (S, Sj, step, _, _, _, _) = candidates[key]
                repeat.add(S.name)
                sma = Sj_sma[cnd_i-1]
                kept.add(cnd_i)
                # cnd_i = best[S.name][0] 
                # hent = Sj
                visited.add((hent.name, oper))
                edits = step.overlap[0]
                for union_idx in [(S.index, S.name), (hent.index, hent.name)]:
                    for k in list(union_cache):
                        if union_idx == tuple(k[:2]):
                            union_cache.pop(k)


                if oper == strategy.SPLIT:

                    obj = edits
                    if groups[S.name] and groups[hent.name]:
                        obj = objective.split(
                            groups[S.name], groups[hent.name], overlap=edits
                        )

                    if obj >= 0.0:
                        # we get here if we lazily add many params, so now 
                        # some will no longer satisfy the constraint
                        kept.remove(cnd_i)
                        repeat.remove(S.name)

                        n_list = trees.tree_index_node_remove_by_name(cst.hierarchy.index, hent.name)
                        for n in n_list:
                            recalc = True
                            cst.hierarchy.subgraphs.pop(n.index)
                            cst.hierarchy.smarts.pop(n.index)
                            micro_count[m.name] -= 1
                            micro_count[m.name] = max(micro_count[m.name], 0)
                            macro_count[m.name] -= 1
                            macro_count[m.name] = max(macro_count[m.name], 0)
                        print(
                            f"\n>>>>> Skipping parameter {cnd_i:4d}/{cnd_n}",
                            hent.name,
                            "parent",
                            S.name,
                            "Objective",
                            f"{X:10.5f}",
                            "Delta",
                            f"{dX:10.5f}",
                            f"Partition {len(cst.mappings[S.name])}|{len(cst.mappings[hent.name])}",
                        )
                        print(" >>>>>", key, f"Local dObj {obj:10.5f}", sma, end="\n\n")
                    else:
                        success = True
                        added = True
                        print(
                            f"\n>>>>> New parameter {cnd_i:4d}/{cnd_n}",
                            hent.name,
                            "parent",
                            S.name,
                            "Objective",
                            f"{X:10.5f}",
                            "Delta",
                            f"{dX:10.5f}",
                            f"Partition {len(cst.mappings[S.name])}|{len(cst.mappings[hent.name])}",
                        )
                        print(" >>>>>", key, f"Local dObj {obj:10.5f}", sma, end="\n\n")

                        repeat.add(hent.name)
                        step_tracker[hent.name] = 0

                elif oper == strategy.MERGE:

                    if hent.name in step_tracker:
                        step_tracker.pop(hent.name)
                    else:
                        print("WARNING", hent.name, "missing from the tracker")

                    visited.add((S.name, oper))
                    visited.remove(hent.name, oper)

                    above = cst.hierarchy.index.above.get(S.index)
                    if above is not None:
                        repeat.add(cst.hierarchy.index.nodes[above].name)

                    success = True
                    added = True
                    print(
                        f">>>>> Delete parameter {cnd_i:4d}/{cnd_n}",
                        hent.name,
                        "parent",
                        S.name,
                        "Objective",
                        f"{X:10.5f}",
                        "Delta",
                        f"{dX:10.5f}",
                    )
                    print(" >>>>>", key, f"Local dObj {obj:10.5f}", sma, end="\n\n")

            if recalc:
                print("Detected change in result")
                print("Operations per parameter for this micro:")
                print(micro_count)
                print(f"Micro total: {sum(micro_count.values())}")

                print("Operations per parameter for this macro:")
                print(macro_count)
                print(f"Macro total: {sum(macro_count.values())}")

                new_assignments = labeler.assign(hidx, gcd, smiles, topo)
                # print(datetime.datetime.now(), '*** 5')
                new_match = clustering_build_assignment_mappings(hidx, new_assignments)

                cst = smarts_clustering(hidx, new_assignments, new_match)

                groups = clustering_build_ordinal_mappings(cst, sag)


            for ei, e in enumerate(
                tree_iterators.tree_iter_dive(
                    cst.hierarchy.index,
                    trees.tree_index_roots(cst.hierarchy.index),
                )
            ):
                s = trees.tree_index_node_depth(cst.hierarchy.index, e)
                obj_repo = ""
                # if groups[e.name]:
                obj_repo = objective.report(groups[e.name])
                print(
                    f"** {s:2d} {ei:3d} {e.name:4s}",
                    obj_repo,
                    cst.hierarchy.smarts.get(e.index),
                )

            mod_lbls = cluster_assignment.smiles_assignment_str_modified(
                cur_cst.group.assignments, cst.group.assignments
            )
            repeat.update(mod_lbls)
            cur_cst = cst
            cst = None
            X0 = X

        # wq.close()
        # wq = None

        if strategy.macro_accept_max_total > 0 and n_added > 0:
            strategy.repeat_step()

        print(f"There were {n_added} successful operations")
        cst = cur_cst
        for ei, e in enumerate(
            tree_iterators.tree_iter_dive(
                cst.hierarchy.index,
                trees.tree_index_roots(cst.hierarchy.index),
            )
        ):
            s = trees.tree_index_node_depth(cst.hierarchy.index, e)
            obj_repo = ""
            # if groups[e.name]:
            obj_repo = objective.report(groups[e.name])
            print(
                f"** {s:2d} {ei:3d} {e.name:4s}",
                obj_repo,
                cst.hierarchy.smarts.get(e.index),
            )

        print(f"{datetime.datetime.now()} Visited", visited)
        for name in (node.name for node in cst.hierarchy.index.nodes.values()):
            if name not in step_tracker:
                continue

            if name not in repeat:
                step_tracker[name] = max(strategy.cursor, step_tracker[name])
            else:
                print(f"Assignments changed for {name}, will retarget")
                step_tracker[name] = 0

        pickle.dump([sag, cst, strategy], open("chk.cst.p", "wb"))

    new_assignments = labeler.assign(cst.hierarchy, gcd, smiles, topo)
    mappings = clustering_build_assignment_mappings(
        cst.hierarchy, new_assignments
    )
    cst = smarts_clustering(cst.hierarchy, new_assignments, mappings)
    pickle.dump([sag, cst, strategy], open("chk.cst.p", "wb"))

    ended = datetime.datetime.now()

    print(f"Start time: {started}")
    print(f"End   time: {ended}")

    wq.close()
    wq = None

    # gc.enable()
    return cst


def smarts_clustering_find_max_depth(
    group: assignments.structure_assignment_group, maxdepth, gcd=None
) -> int:
    prev_problems = 99999999
    N = 0
    struct_lbls = {}
    max_depth = 0

    for i in range(0, maxdepth+1):
        topo = group.topology
        assignments = {}
        struct_lbls = {}
        for moli, mol in enumerate(group.assignments, 0):
            print(
                f"Assigning molecule {moli+1:5d}/{len(group.assignments):d} at depth {i}",
                end="\r",
            )
            N += len(mol.selections)
            for molj, (selection, lbl) in enumerate(mol.selections.items()):
                idx = tuple((selection[j] for j in topo.primary))
                atom = graphs.graph_to_structure(mol.graph, idx, topo)
                suc = mapper.mapper_smarts_extend(
                    configs.smarts_extender_config(i, i, False), [atom]
                )
                atom = graphs.structure_remove_unselected(atom)
                atom_lbls = struct_lbls.get(atom, dict())
                assignments[(moli, molj)] = lbl

                if lbl not in atom_lbls:
                    atom_lbls[lbl] = []

                atom_lbls[lbl].append((moli, molj, idx, hash(atom)))
                struct_lbls[atom] = atom_lbls

                # print((moli, molj, idx), gcd.smarts_encode(atom))
        print()
        print("Labels per unique structure that need more depth")
        problems = set()
        for j, (h, lbl) in enumerate(struct_lbls.items()):
            if len(lbl.values()) > 1:
                problems.add(tuple((h, *sorted(set(lbl.keys())))))

        print(
            f"There are {len(set(struct_lbls))}/{N} unique structures at depth",
            i,
        )
        print(f"There are {len(problems)} problems:")
        for problem in problems:
            print("SMARTS", gcd.smarts_encode(problem[0]))
            print("DATA")
            pprint(problem[1:])
            print("ASSN:")
            pprint(struct_lbls[problem[0]])
        # print(problems)
        if len(problems) < prev_problems:
            max_depth = i
            prev_problems = len(problems)

        if len(problems) == 0:
            break
        bad_mols = set()
        if gcd:
            for moli, mol in enumerate(group.assignments, 0):
                # print(mol.smiles)
                for molj, (selection, lbl) in enumerate(mol.selections.items()):
                    # idx = tuple((selection[i] for i in topo.primary))
                    # atom = graphs.graph_to_structure(mol.graph, idx, topo)
                    # suc = mapper.mapper_smarts_extend(
                    #     configs.smarts_extender_config(i, i, True), [atom]
                    # )
                    these_probs = []
                    for problem in problems:
                        if lbl in problem[1:]:
                            bad_mols.add(moli)
                            # these_probs.append(problem)
                    # print(
                    #     "   ",
                    #     hash(problem[0]), gcd.smarts_encode(problem[0]),
                    #     these_probs,
                    #     selection,
                    #     lbl,
                    #     hash(atom),
                    #     gcd.smarts_encode(atom),
                    # )

        print(f"There were {len(bad_mols)}/{len(group.assignments)} that could not be distinguished:")
        for mid in bad_mols:
            print(mid)

        if all(len(x.values()) == 1 for x in struct_lbls.values()):
            break

    if any(len(x) > 1 for x in struct_lbls.values()):
        print(
            "WARNING",
            "there are environments with multiple labels. will be"
            "impossible to split environment",
        )
    print("Max depth is set to", max_depth)
    return max_depth


def clustering_collect_structures(A, matches, topo, extend):
    a = [graphs.subgraph_as_structure(A[i], topo) for i in matches]

    for ai in a:
        ai.select = tuple(ai.select[x] for x in ai.topology.primary)

    if extend.depth_max > 0:
        mapper.mapper_smarts_extend(extend, a)

    a = {i: x for i, x in zip(matches, a)}
    return a


def clustering_update_assignments(
    group: assignments.structure_assignment_group, match
) -> assignments.structure_assignment_group:
    new_group = group.copy()
    inverted_match = {x: lbl for lbl, y in match.items() for x in y}
    topo = group.topology

    for i, assns_i in enumerate(new_group.assignments):
        for idx in assns_i.selections:
            primary_idx = tuple((idx[j] for j in topo.primary))
            assns_i.selections[idx] = inverted_match[(i, primary_idx)]

    return new_group


def clustering_initial_conditions(
    gcd, sag: assignments.smiles_assignment_group, hidx=None, labeler=None, prefix="p"
):
    topo = sag.topology
    group_prefix_str = prefix

    if hidx is None:
        hidx = hierarchies.structure_hierarchy(trees.tree_index(), {}, {}, topo)

        hidx.index.node_add(None, trees.tree_node(0, tuple((0,0,0)), "parameter", "p0"))

        if topo == topology.atom:
            S0 = gcd.smarts_decode("[*:1]")
        elif topo == topology.bond:
            S0 = gcd.smarts_decode("[*:1]~[*:2]")
        elif topo == topology.angle:
            S0 = gcd.smarts_decode("[*:1]~[*:2]~[*:3]")
        elif topo == topology.torsion:
            S0 = gcd.smarts_decode("[*:1]~[*:2]~[*:3]~[*:4]")
        elif topo == topology.outofplane:
            S0 = gcd.smarts_decode("[*:1]~[*:2](~[*:3])~[*:4]")


        hidx.subgraphs[0] = graphs.structure_to_subgraph(S0)
        if gcd:
            hidx.smarts[0] = gcd.smarts_encode(S0)
        group_name = group_prefix_str + "0"

        assn = []
        for sag_i in sag.assignments:
            sels = {}
            for idx in sag_i.selections:
                sels[idx] = group_name
            assn.append(
                cluster_assignment.smiles_assignment_str(sag_i.smiles, sels)
            )

        groups: assignments.assignment_mapping = {
            group_name: list(
                x for y in sag.assignments for x in list(y.selections.values())
            )
        }
    else:
        assn = []
        new_assignments = labeler.assign(hidx, gcd, [x.smiles for x in sag.assignments], topo)
        all_lbls = set()
        for sag_i, lbls in zip(sag.assignments, new_assignments.assignments):
            sels = {}
            for idx in sag_i.selections:
                sels[idx] = lbls.selections[idx]
                all_lbls.add(sels[idx])

            assn.append(
                cluster_assignment.smiles_assignment_str(sag_i.smiles, sels)
            )
        groups = {lbl: [] for lbl in all_lbls}
        for sag_i, lbls in zip(sag.assignments, new_assignments.assignments):
            for idx, vals in sag_i.selections.items():
                groups[lbls.selections[idx]].append(vals)

        if gcd:
            for idx, smarts in hidx.smarts.items():
                if smarts:
                    hidx.subgraphs[idx] = gcd.smarts_decode(smarts)
                else:
                    hidx.subgraphs[idx] = None


    new_assn_group = assignments.smiles_assignment_group(assn, sag.topology)

    initial_conditions = smarts_clustering(hidx, new_assn_group, groups)
    initial_conditions.group_prefix_str = group_prefix_str

    return initial_conditions


class clustering_collect_split_candidates_ctx:
    ret = None
    assn = None
    backmap = None


def clustering_collect_split_candidates_single_distributed(
    j, matched, shm=None
):

    unmatch = [v for k, v in shm.backmap.items() if k not in matched]
    matched = [shm.backmap[i] for i in matched]

    matched_assn = tuple((shm.assn[i] for i in matched))
    unmatch_assn = tuple((shm.assn[i] for i in unmatch))

    return j, unmatch, matched, matched_assn, unmatch_assn


def clustering_collect_split_candidates_single(j):
    ret = clustering_collect_split_candidates_ctx.ret
    assn = clustering_collect_split_candidates_ctx.assn
    backmap = clustering_collect_split_candidates_ctx.backmap

    Sj = ret.splits[j]
    # bj = ret.shards[i]
    matched = ret.matched_idx[j]

    unmatch = [v for k, v in backmap.items() if k not in matched]
    matched = [backmap[i] for i in matched]

    matched_assn = tuple((assn[i] for i in matched))
    unmatch_assn = tuple((assn[i] for i in unmatch))

    return j, unmatch, matched, matched_assn, unmatch_assn


def clustering_collect_split_candidates_serial(S, ret, step, operation):
    candidates = {}

    for p_j, Sj in enumerate(ret.splits):
        overlaps = step.overlap

        if overlaps is None:
            overlaps = [0]

        for edits in step.overlap:
            # if edits not in candidates:
            # candidates[edits] = []
            key = (edits, None, p_j, operation)
            candidates[key] = (S, Sj, step, None, None, None, None)

    return candidates


def clustering_node_remove_by_name(
    ph: smarts_clustering,
    gcd: codecs.graph_codec,
    assign: assignments.smarts_hierarchy_assignment,
    name: str,
):
    """
    removes the nodes and reassigns tree
    """
    n = trees.tree_index_node_remove_by_name(ph.hierarchy.index, name)
    clustering_assign(ph, gcd, assign)
    ph.mappings.clear()
    return n


def clustering_assign(
    ph: smarts_clustering,
    gcd: codecs.graph_codec,
    assign: assignments.smarts_hierarchy_assignment,
):
    hierarchy = hierarchies.structure_hierarchy_to_smarts_hierarchy(
        ph.hierarchy, gcd
    )
    for assignment in ph.group.assignments:
        smiles = assignment.smiles
        selections = list(assignment.selections)
        matches = assign.assign(hierarchy, gcd, smiles, selections)
        assignment.selections.update(matches)


def clustering_build_assignment_mappings(
    hierarchy: hierarchies.smarts_hierarchy,
    assns: assignments.smiles_assignment_group,
) -> assignments.assignment_mapping:

    mappings = {n.name: [] for n in hierarchy.index.nodes.values()}
    # allow a None if there is no parameter for diagostic purposes
    mappings[None] = []
    for i, ag in enumerate(assns.assignments):
        for sel, lbl in ag.selections.items():
            mappings[lbl].append((i, sel))

    return mappings


def clustering_build_ordinal_mappings(
    initial_conditions: smarts_clustering, stuag, select=None
):
    """
    parameter:data mapping
    """
    mapping = {n.name: [] for n in initial_conditions.hierarchy.index.nodes.values()}
    for a, b in zip(initial_conditions.group.assignments, stuag.assignments):

        assert a.smiles == b.smiles

        for sel, x in a.selections.items():
            assert x is not None
            if select is None or x in select:
                y = b.selections.get(sel)
                mapping[x].append(y)

    return mapping


def clustering_build_label_mappings(
    initial_conditions: smarts_clustering, stuag
):
    mapping = {}
    for a, b in zip(initial_conditions.group.assignments, stuag.assignments):
        sa = a.selections
        sb = b.selections
        for x, y in zip(sa.values(), sb.values()):
            if x not in mapping:
                mapping[x] = set()
            mapping[x].add(y)
    return mapping


def match_group_assignments(
    assignments, topo
) -> Dict[str, List[Tuple[int, Sequence[int]]]]:
    match = {}
    for i, assns_i in enumerate(assignments):
        for idx in assns_i.selections:
            lbl = assns_i.compute(idx)
            idx = tuple((idx[j] for j in topo.primary))
            assns = match.get(lbl, list())
            assns.append((i, idx))
            match[lbl] = assns
    return match


def get_assns(assignments, topo):
    assn = {}
    for i, assns_i in enumerate(assignments):
        for idx, lbl in assns_i.selections.items():
            idx = tuple((idx[j] for j in topo.primary))
            assn[(i, idx)] = lbl
    return assn


def smarts_filter_data(
    gcd: codecs.graph_codec,
    labeler: assignments.smarts_hierarchy_assignment,
    sag: assignments.smiles_assignment_group,
    hierarchy: hierarchies.smarts_hierarchy,
    bounds: Dict[str, Tuple[float, float]]
) -> List[int]:
    smiles = [a.smiles for a in sag.assignments]

    lbl_assn = labeler.assign(hierarchy, gcd, smiles, sag.topology)

    keep = []
    for i, (mol, lbls) in enumerate(zip(sag.assignments, lbl_assn.assignments)):
        valid = True
        for ic, r in mol.selections.items():
            lbl = lbls.selections.get(ic)
            if lbl is None:
                print(f"Warning, bond {ic} was not assigned a label")
                continue
            R = bounds.get(lbl)
            if R is None:
                print(f"Warning, cutoff for {lbl} was not provided")
                continue
            sma = gcd.smarts_encode(
                graphs.structure_remove_unselected(graphs.graph_to_structure(
                    mol.graph, ic, sag.topology
                ))
            )

            lR, rR = R
            status = "filter:"
            if (
                (rR is not None and r[0] >= rR) or
                (lR is not None and r[0] <= lR)
            ):
                valid = False
            else:
                status = "keep:  "

            print(
                status, i, ic, r, R, lbl, sma
            )
        if valid:
            keep.append(i)

    return keep

def smarts_filter_bond_lengths(
    gcd: codecs.graph_codec,
    labeler: assignments.smarts_hierarchy_assignment,
    sag: assignments.smiles_assignment_group,
    hierarchy: hierarchies.smarts_hierarchy,
    cutoffs: Dict[str, float]
) -> List[int]:

    bounds = {k: [0, v] for k, v in cutoffs.items()}
    return smarts_filter_data(gcd, labeler, sag, hierarchy, bounds)
