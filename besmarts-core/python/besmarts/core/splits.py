"""
besmarts.core.splits

Chemical perception in SMARTS. Functions to find numerical splits based on
bit iteration, and analytical splits based on a predefined partitioning.
"""

import array
import os
import itertools
import math
from typing import List, Set, Tuple, Sequence
import multiprocessing.pool

import datetime

from besmarts.core.topology import structure_topology
from besmarts.core.graphs import (
    subgraph,
    structure,
)
from besmarts.core import (
    mapper,
    configs,
    graphs,
    graph_visitors,
    codecs,
    compute,
    arrays,
)
from besmarts.codecs import codec_native

from besmarts.core.returns import return_value, success

from besmarts.core.configs import (
    smarts_perception_config,
    smarts_splitter_config,
)


class split_return_type:
    __slots__ = (
        "splits",
        "shards",
        "matched_idx",
        "unmatch_idx",
        "subgraphs",
        "topology",
    )

    def __init__(
        self,
        splits: Sequence[subgraph],
        shards: Sequence[subgraph],
        matched_idx: Sequence[Sequence[int]],
        unmatch_idx: Sequence[Sequence[int]],
        subgraphs: Sequence[subgraph],
        topology: structure_topology,
    ):
        self.splits: Sequence[subgraph] = splits
        self.shards: Sequence[subgraph] = shards
        self.matched_idx: Sequence[Sequence[int]] = matched_idx
        self.unmatch_idx: Sequence[Sequence[int]] = unmatch_idx
        self.subgraphs: Sequence[subgraph] = subgraphs
        self.topology: structure_topology = topology


class process_split_ctx:
    S0 = None
    splitter = None
    A = None
    pool = None


def process_split_matches_distributed(Sj, indices, shm=None):

    # try to send back the minimum amount of data

    icd = shm.icd
    G = shm.G
    s = shm.selections
    values = (j for j in indices if mapper.mapper_match(
            graphs.graph_to_structure(
                icd.graph_decode(G[s[j][0]]),
                s[j][1],
                Sj.topology),
            Sj
            )
    )
    ret = array.array(indices.typecode, values)
    # t = time.perf_counter() - t0
    # print(f"Time taken t={t:.6f} answer is {ret}")
    return ret


def assert_bonds(g, n=5, line="None"):
    for n in g.select:
        con = graphs.subgraph_connection(g, n)
        assert len(con) < 5


def dprint(*args, **kwargs):
    if kwargs.pop("on", False):
        print(*args, **kwargs)


def split_all_partitions(
    topology: structure_topology,
    perception: smarts_perception_config,
    fragments: List[subgraph],
    assignments: List[str],
    gcd=None,
    maxmoves=0,
) -> return_value[List[Tuple[structure, List[int]]]]:
    """
    Find a shard that induces the given partition in a sequence of fragments.
    Up to two shards can be found; one will match everything in the specified
    partition, and the other will match the fragments not in the partition.

    Parameters
    ----------
    topology: structure_topology
        The topology of the structures to split
    perception: smarts_perception_config
        The settings that control the search depth
    fragments: List[subgraph]
        The list of subgraphs that are the target of splitting
    assignments: List[str]
        The label for each subgraph. The function will try to create partitions
        that group these labels together
    gcd: graph_codec
        A graph codec
    maxmoves: int
        The number of subgraphs that we can allow to move around to arrive at
        a solution

    Returns
    -------
    A return_value containing pairs of structures and indices of the subgraphs
    that belong to the structures. The structures are the SMARTS that cluster
    the subgraphs..
    """

    lbls = []
    for x in assignments:
        if x not in lbls:
            lbls.append(x)

    bitmin = perception.splitter.bit_search_min
    bitmax = perception.splitter.bit_search_limit
    # bitmax = max(1, len(lbls) // 2)

    results = []
    seen = set()
    for b in range(bitmin, bitmax + 1):
        if b > len(lbls):
            continue
        if b == len(lbls) and maxmoves == 0:
            continue
        for combo in itertools.combinations(lbls, b):
            deselect = tuple(sorted([x for x in lbls if x not in combo]))
            if tuple(sorted(combo)) in seen or deselect in seen:
                continue
            seen.add(tuple(sorted(combo)))
            seen.add(deselect)
            print(
                "Direct on",
                b,
                "combo",
                tuple(sorted(combo)),
                "depth",
                perception.extender.depth_min,
                perception.extender.depth_max,
            )
            matched = set(
                [i for i, lbl in enumerate(assignments) if lbl in combo]
            )
            unmatch = set(
                [i for i, lbl in enumerate(assignments) if i not in matched]
            )
            if not (matched and unmatch):
                print("Occluding. Skipped")
                continue
            ret = split_partition(
                topology,
                perception,
                fragments,
                matched,
                gcd=gcd,
                maxmoves=maxmoves,
            )
            (
                lhs,
                rhs,
                lhs_removeA,
                lhs_removeB,
                rhs_removeA,
                rhs_removeB,
            ) = ret.value
            matched.difference_update(lhs_removeA)
            matched.update(lhs_removeB)
            unmatch.difference_update(rhs_removeB)
            unmatch.update(rhs_removeA)

            matched.difference_update(rhs_removeA)
            matched.update(rhs_removeB)
            unmatch.difference_update(lhs_removeB)
            unmatch.update(lhs_removeA)

            # if lhs or rhs:
            #     results.append((lhs, rhs, matched, unmatch))
            if lhs:
                results.append((lhs, rhs, matched, unmatch))
            elif rhs:
                results.append((rhs, lhs, unmatch, matched))

    if not results and maxmoves > 0:
        print("Brute force on", len(fragments))
        indices = range(len(fragments))
        for b in range(1, len(indices) - 1):
            combos = list(itertools.combinations(indices, b))
            for i, matched in enumerate(combos, 1):
                print(
                    "Brute force",
                    b,
                    i,
                    len(combos),
                    "depth",
                    perception.extender.depth_min,
                    perception.extender.depth_max,
                )
                matched = set(matched)
                unmatch = set(indices).difference(matched)

                ret = split_partition(
                    topology,
                    perception,
                    fragments,
                    matched,
                    gcd=gcd,
                    maxmoves=maxmoves,
                )
                (
                    lhs,
                    rhs,
                    lhs_removeA,
                    lhs_removeB,
                    rhs_removeA,
                    rhs_removeB,
                ) = ret.value
                matched.difference_update(lhs_removeA)
                matched.update(lhs_removeB)
                unmatch.difference_update(rhs_removeB)
                unmatch.update(rhs_removeA)
                matched.difference_update(rhs_removeA)
                matched.update(rhs_removeB)
                unmatch.difference_update(lhs_removeB)
                unmatch.update(lhs_removeA)
                if lhs:
                    print("Hit")
                    results.append((lhs, rhs, matched, unmatch))
                    # return success(results)
                elif rhs:
                    print("Hit")
                    results.append((rhs, lhs, unmatch, matched))
                    # return success(results)

    return success(results)


def split_partition(
    topology: structure_topology,
    perception: smarts_perception_config,
    fragments: List[subgraph],
    partition: Set[int],
    gcd=None,
    maxmoves=0,
) -> return_value[Tuple[structure, structure, List, List]]:
    """
    Find a shard that induces the given partition in a sequence of fragments.
    Up to two shards can be found; one will match everything in the specified
    partition, and the other will match the fragments not in the partition.
    """

    add_nodes = True
    lhs, rhs = None, None
    partition = set(partition)
    bestlhs = None
    bestrhs = None
    best_lhs_removeA = set()
    best_lhs_removeB = set()
    best_rhs_removeA = set()
    best_rhs_removeB = set()
    lhs_nummoves = 0
    rhs_nummoves = 0
    Ta = None
    Tb = None

    union_config = configs.mapper_config(0, False, "high")
    config = configs.mapper_config(3, False, "high")
    intr_config = configs.mapper_config(0, True, "high")

    # gcd = None

    for depth in range(
        perception.extender.depth_min, perception.extender.depth_max + 1
    ):
        A = [graphs.structure_copy(a) for a in fragments]
        for a in A:
            a.topology = topology
            a.select = tuple((a.select[i] for i in topology.primary))
        lhs_removeA = set()
        lhs_removeB = set()

        rhs_removeA = set()
        rhs_removeB = set()

        extender = perception.extender.copy()
        extender.depth_max = depth
        extender.depth_min = depth

        suc = mapper.mapper_smarts_extend(extender, A)
        A = [graphs.structure_remove_unselected(a) for a in A]
        if not suc and depth > perception.extender.depth_min:
            break

        Ai = list([A[i] for i in partition])
        # if gcd:
        #     for ai in Ai:
        #         print("LHS_MATCH: ", gcd.smarts_encode(ai))

        # a = Ai[0]

        # if Ta is None:
        #     Ta = mapper.structure_mapper(a, add_nodes=)

        # for i, ai in enumerate(Ai):
        #     Ta.add(ai)

        ref_a = mapper.union_list(Ai, union_config, max_depth=depth)
        a = mapper.intersection_list(
            Ai, config, max_depth=depth, reference=ref_a
        )
        # a = ref_a
        # graphs.structure_print(a)
        # Tb = None
        if gcd:
            print("LUN: ", gcd.smarts_encode(ref_a))
            print("LHS: ", gcd.smarts_encode(a))
        Ta = None

        part_diff = [i for i in range(len(A)) if i not in partition]
        Bi = list([A[i] for i in part_diff])
        b = Bi[0]

        # if gcd:
        #     for ai in Bi:
        #         print("RHS_MATCH: ", gcd.smarts_encode(ai))
        # if Tb is None:
        #     Tb = mapper.structure_mapper(b, add_nodes=True)
        # for i, bi in enumerate(Bi):
        #     Tb.add(bi)

        ref_b = mapper.union_list(Bi, union_config, max_depth=depth)
        b = mapper.intersection_list(
            Bi, config, max_depth=depth, reference=ref_b
        )
        # b = ref_b
        # Tb = None
        # print("bu: ", gcd.smarts_encode(ref))
        if gcd:
            print("RUN: ", gcd.smarts_encode(ref_b))
            print("RHS: ", gcd.smarts_encode(b))
        # print("b: ", gcd.smarts_encode(b))

        lhs: structure = mapper.difference(a, b, config)
        if gcd:
            print("LHS_DIFF: ", gcd.smarts_encode(lhs))

        if graphs.graph_any(lhs):
            lhs = a
            graphs.subgraph_invert_null(lhs)
            if gcd:
                print("LHS_INVE: ", gcd.smarts_encode(lhs))

            # lhs = mapper.intersection_list_parallel([lhs] + Ai, intr_config, max_depth=depth, reference=lhs)
            # lhs = graphs.subgraph_invert_null(lhs)
            # if gcd:
            #     print("LHS_UNIO: ", gcd.smarts_encode(lhs))

            valid = True
            for i in partition:
                ai = A[i]
                if not mapper.mapper_match(ai, lhs):
                    lhs_removeA.add(i)
                    # print("LHS Didn't match but should:", i, gcd.smarts_encode(ai))

            for i in part_diff:
                ai = A[i]
                if mapper.mapper_match(ai, lhs):
                    lhs_removeB.add(i)
                    # print("LHS Did match but shouldn't:", i, gcd.smarts_encode(ai))

            if len(lhs_removeA) + len(lhs_removeB) > maxmoves:
                valid = False

            if not valid:
                lhs = None

        else:
            lhs = None

        rhs: structure = mapper.difference(b, a, config)
        if gcd:
            print("RHS_DIFF: ", gcd.smarts_encode(rhs))
        if graphs.graph_any(rhs):
            rhs = b
            graphs.subgraph_invert_null(rhs)
            if gcd:
                print("RHS_INVE: ", gcd.smarts_encode(rhs))

            # rhs = mapper.intersection_list_parallel([rhs] + Bi, intr_config, max_depth=depth, reference=rhs)
            # rhs = graphs.subgraph_invert_null(rhs)
            if gcd:
                print("RHS_INTR: ", gcd.smarts_encode(rhs))
            valid = True
            for i in partition:
                if len(rhs_removeA) > maxmoves:
                    break
                ai = A[i]
                if mapper.mapper_match(ai, rhs):
                    rhs_removeA.add(i)
                    # print("RHS Did match but shouldn't:", i, gcd.smarts_encode(ai))

            for i in part_diff:
                if len(rhs_removeA) + len(rhs_removeB) > maxmoves:
                    break
                ai = A[i]
                if not mapper.mapper_match(ai, rhs):
                    rhs_removeB.add(i)
                    # print("RHS Didn't match but should:", i, gcd.smarts_encode(ai))
            if len(rhs_removeA) + len(rhs_removeB) > maxmoves:
                valid = False

            if not valid:
                rhs = None
        else:
            rhs = None

        winner = False
        this_lhs_nummoves = len(lhs_removeA) + len(lhs_removeB)
        if (this_lhs_nummoves == 0) and lhs:
            relabel = {x: i for i, x in enumerate(lhs.select, 1)}
            for i, x in enumerate(lhs.nodes, len(lhs.select) + 1):
                if x not in lhs.select:
                    relabel[x] = i
            lhs = graphs.structure_relabel_nodes(lhs, relabel)
            if gcd:
                print("BESTLHS: ", gcd.smarts_encode(lhs))
            bestlhs = graphs.structure_copy(lhs)
            best_lhs_removeA = lhs_removeA
            best_lhs_removeB = lhs_removeB
            winner = True
            lhs_nummoves = this_lhs_nummoves
        this_rhs_nummoves = len(rhs_removeA) + len(rhs_removeB)
        if (this_rhs_nummoves == 0) and rhs:
            relabel = {x: i for i, x in enumerate(rhs.select, 1)}
            for i, x in enumerate(rhs.nodes, len(rhs.select) + 1):
                if x not in rhs.select:
                    relabel[x] = i
            rhs = graphs.structure_relabel_nodes(rhs, relabel)
            if gcd:
                print("BESTRHS: ", gcd.smarts_encode(rhs))
            bestrhs = graphs.structure_copy(rhs)
            best_rhs_removeA = rhs_removeA
            best_rhs_removeB = rhs_removeB
            winner = True
            rhs_nummoves = this_rhs_nummoves
        if winner:
            break

        lhs_better = bool(lhs and this_lhs_nummoves <= lhs_nummoves)
        rhs_better = bool(rhs and this_rhs_nummoves <= rhs_nummoves)
        if lhs_better and rhs_better:
            if lhs_nummoves <= rhs_nummoves:
                rhs_better = False
            else:
                lhs_better = False
        first = bool((bestlhs is None and bestrhs is None) and (lhs or rhs))
        if lhs and (first or lhs_better):
            relabel = {x: i for i, x in enumerate(lhs.select, 1)}
            for i, x in enumerate(lhs.nodes, len(lhs.select) + 1):
                if x not in lhs.select:
                    relabel[x] = i
            lhs = graphs.structure_relabel_nodes(lhs, relabel)
            if gcd:
                print("BESTLHS: ", gcd.smarts_encode(lhs))
            bestlhs = graphs.structure_copy(lhs)

            best_lhs_removeA = lhs_removeA
            best_lhs_removeB = lhs_removeB

            best_rhs_removeA = rhs_removeA
            best_rhs_removeB = rhs_removeB
            lhs_nummoves = this_lhs_nummoves
            rhs_nummoves = this_rhs_nummoves

        if rhs and (first or rhs_better):
            relabel = {x: i for i, x in enumerate(rhs.select, 1)}
            for i, x in enumerate(rhs.nodes, len(rhs.select) + 1):
                if x not in rhs.select:
                    relabel[x] = i
            rhs = graphs.structure_relabel_nodes(rhs, relabel)
            if gcd:
                print("BESTRHS: ", gcd.smarts_encode(rhs))
            bestrhs = graphs.structure_copy(rhs)

            best_lhs_removeA = lhs_removeA
            best_lhs_removeB = lhs_removeB
            best_rhs_removeA = rhs_removeA
            best_rhs_removeB = rhs_removeB
            lhs_nummoves = this_lhs_nummoves
            rhs_nummoves = this_rhs_nummoves

    return success(
        (
            bestlhs,
            bestrhs,
            best_lhs_removeA,
            best_lhs_removeB,
            best_rhs_removeA,
            best_rhs_removeB,
        )
    )


def make_branches(
    ref: structure,
    G: List[graphs.graph],
    selections,
    max_branch_depth,
    max_branches=6,
    gcd=None,
    icd=None):

    if not G or not selections:
        return []

    if gcd is None:
        cdcs = codec_native.primitive_codecs_get()
        atom_p = list(codec_native.primitive_codecs_get_atom())
        bond_p = list(codec_native.primitive_codecs_get_bond())
        gcd = codec_native.graph_codec_native(cdcs, atom_p, bond_p)

    if icd is None:
        icd = codecs.intvec_codec(
            gcd.primitive_codecs,
            gcd.atom_primitives,
            gcd.bond_primitives
        )

    config = configs.mapper_config(1, False, "high")
    branch_length = max_branch_depth
    branches = []
    new_branches = []
    trunks = None
    branch_count = 0

    # print(f"ref")
    # print(gcd.smarts_encode(ref))

    _ref = ref
    topo = ref.topology
    trunks = [(None, _ref)]

    for depth in range(max_branch_depth, max_branch_depth + 1):
    # for depth in range(1, max_branch_depth + 1):
        print(
            datetime.datetime.now(),
            f"Generating branched for depth={depth} trunks={len(trunks)}",
        )

        for _, ref in list(trunks):
            print(
                datetime.datetime.now(),
                f"Unioning {len(selections)} unique branches",
            )
            branched_group_ref = mapper.union_list_parallel(
                G, selections, topo,
                # prims_pruned,
                config,
                max_depth=depth,
                reference=None,
                icd=icd
            )
            # print(f"branched_group_ref")
            # print(gcd.smarts_encode(branched_group_ref))

            Mref = mapper.map_to(branched_group_ref, ref, add_nodes=1)
            branched_group, ref_mapped, Mref = (
                Mref.G,
                Mref.H,
                Mref.map,
            )

            Mparam = {k: v for k, v in Mref.items() if v in ref.select}

            trunk = ref
            if trunks is None:
                trunks = [(Mref, ref)]

            branch_depth = 0
            added = True

            print(datetime.datetime.now(), f"Remapping branch to reference")
            T = mapper.map_to(branched_group_ref, trunk, add_nodes=1, pool=True)

            branched_group, trunk, M = T.G, T.H, T.map
            M = {v: k for k, v in M.items() if v in trunk.select}
            trunk = mapper.union(trunk, branched_group, map=M)
            m = {k: k for k in _ref.nodes}

            i = 0
            print(
                datetime.datetime.now(),
                f"Generating branched for depth={depth} trunks={len(trunks)}",
            )
            # print(f"Trunk")
            # print(gcd.smarts_encode(trunk))
            # this takes the trunk and starts with only the S0 nodes
            all_branches = list(
                graphs.structure_branch(
                    trunk, m, max_branches, max_branch_depth
                )
            )
            print(
                datetime.datetime.now(),
                f"Generating branched for depth={depth} trunks={len(trunks)} branch={len(all_branches)}",
            )

            for bi, branch in enumerate(all_branches, 1):
                # parallelize this
                # reduce_branch
                i += 1
                configurations = mapper.map_to(
                    graphs.structure_clear(branch),
                    graphs.structure_clear(trunk),
                    add_nodes=0,
                    return_all=True,
                )

                if not configurations:
                    continue

                # print(f"BRANCH {bi}")
                # print(gcd.smarts_encode(branch))
                new_branch = branch

                for Tmult in configurations:
                    new_branch = mapper.union(
                        new_branch, Tmult.H, map=Tmult.map
                    )
                    for n in _ref.nodes:
                        if n in new_branch.nodes:
                            new_branch.nodes[n].clear()
                    for n in _ref.edges:
                        if n in new_branch.edges:
                            new_branch.edges[n].clear()

                skip = False
                for i, c in new_branch.nodes.items():
                    for n, p in c.primitives.items():
                        if p.all():
                            skip = True
                            break
                    if skip:
                        break
                if skip:
                    continue

                if not graphs.graph_any(new_branch):
                    continue

                result = (m, graphs.structure_copy(new_branch))
                combined = False
                skip = False


                for i, existing_branch in enumerate(new_branches):
                    if len(result[1].nodes) != len(existing_branch[1].nodes):
                        continue

                    T = mapper.map_to(
                        result[1], existing_branch[1], add_nodes=0, pool=True
                    )
                    # if set(result.nodes).difference(T.G.nodes) or set(T.map).symmetric_difference(result[1].nodes):
                    #     continue

                    if len(T.G.nodes) == len(result[1].nodes) and len(T.H.nodes) == len(existing_branch[1].nodes):
                        # print("COMBINED A + B = C:")
                        # print(gcd.smarts_encode(existing_branch[1]))
                        # print(gcd.smarts_encode(result[1]))
                        # print("U: ", gcd.smarts_encode(T.G))
                        # print("U: ", gcd.smarts_encode(T.H))
                        union_result = mapper.union(T.G, T.H, map=T.map)
                        if hash(union_result) != hash(new_branches[i][1]):
                            print(
                                datetime.datetime.now(),
                                f"Update {i:4d}: trunks={len(trunks)} d={depth} branches={bi}/{len(all_branches)} new={len(new_branches)} {gcd.smarts_encode(new_branches[i][1])}",
                            )
                        new_branches[i] = (m, union_result)
                        combined = True
                        # print(gcd.smarts_encode(new_branches[i][1]))
                        break

                if not combined:
                    new_branches.append(result)

                if branch_count != len(new_branches):
                    print(
                        datetime.datetime.now(),
                        f"New    {branch_count:4d}: trunks={len(trunks)} d={depth} branches={bi}/{len(all_branches)} new={len(new_branches)} {gcd.smarts_encode(new_branch)}",
                    )
                    branch_count = len(new_branches)
                continue

    return new_branches

def split_branches(
    ref: structure,
    G: List[graphs.graph],
    selections,
    max_branch_depth,
    max_branches=6,
    primitives=None,
    icd=None
) -> List[structure]:

    branches = []
    new_branches = make_branches(ref, G, selections, max_branch_depth, max_branches=max_branches, icd=icd)
    print(
        datetime.datetime.now(),
        f"Generating branched bits over branches={len(new_branches)}",
    )

    if len(new_branches) > 0:
        branched_bits = []
        for mapping, branch in new_branches:
            new_bits = list(
                (bit, mapping)
                for bit in graph_visitors.structure_iter_bits(
                    branch, iter_inverse=True, skip_ones=False, primitives=primitives
                )
            )
            branched_bits.extend(new_bits)

        branches.extend(branched_bits)

    return branches


def process_split(pack):
    S0 = process_split_ctx.S0
    pool = process_split_ctx.pool

    try:
        bit, mapping = pack[0]
    except:
        pack = (pack,)

    bj = None
    mj = None
    shard = None

    for bit, mapping in pack:
        if bj is None:
            bj = bit
            mj = mapping
            # invert = {v: k for k, v in mj.items() if v is not None}
        else:
            composed_map = None
            bj = mapper.union(
                bj,
                bit,
                configs.mapper_config(1, False, "low"),
                map=composed_map,
                pool=pool,
            )
    shard = bj

    if not mj:
        mj = None

    T1 = mapper.map_to(shard, S0, add_nodes=1, fill=True, skip=mj, pool=pool)
    shard = T1.G

    M = T1.map
    T2 = mapper.mapped_type(
        graphs.structure_copy(shard),
        graphs.structure_copy(T1.H),
        dict(M),
    )

    general = process_split_intersect(pack, T1)
    specific = process_split_difference(pack, T2)

    # print("RETURNING", general[2], general[4], specific[2], specific[4])

    return general, specific


def split_shm_load(shm):
    """
    this is the process init, so called by pool init, or maybe even during wsr creation
    """
    process_split_ctx.A = shm.get_A()
    process_split_ctx.splitter = shm.get_splitter()
    process_split_ctx.S0 = shm.get_S0()



def process_split_general_distributed(pack, shm=None):

    S0 = shm.S0

    _, T1 = process_split_distributed(pack, shm=shm)

    general = process_split_difference_distributed(pack, T1, shm.S0, shm.splitter, shm.G, shm.selections, shm.icd)
    # print("RETURNING", general[2], general[4], specific[2], specific[4])

    return (general,)


def process_split_specific_distributed(pack, shm=None):

    shard, T1 = process_split_distributed(pack, shm=shm)

    specific = process_split_intersect_distributed(
        pack, T1, shm.S0, shm.splitter, shm.G, shm.selections, shm.icd
    )

    # print("RETURNING", general[2], general[4], specific[2], specific[4])

    return (specific,)


def process_split_distributed(pack, shm=None):

    splitter = shm.splitter
    S0 = shm.S0
    pool = None

    try:
        bit, mapping = pack[0]
    except Exception:
        pack = (pack,)

    bj = None
    mj = None
    shard = None

    for bit, mapping in pack:
        if bj is None:
            bj = bit
            mj = mapping
            # invert = {v: k for k, v in mj.items() if v is not None}
        else:
            composed_map = None
            bj = mapper.union(
                bj,
                bit,
                configs.mapper_config(1, False, "low"),
                map=composed_map,
                pool=pool,
            )
    shard = bj

    if not mj:
        mj = None

    T1 = mapper.map_to(shard, S0, add_nodes=1, fill=True, skip=mj, pool=None)
    return shard, T1


def process_split_difference(pack, T):
    S0 = process_split_ctx.S0
    splitter = process_split_ctx.splitter
    A = process_split_ctx.A
    pool = process_split_ctx.pool

    try:
        bit, mapping = pack[0]
    except:
        pack = (pack,)

    makes_split = False
    shard = T.G
    matches = None
    shard.cache["hash"] = None

    mj = {v: k for k, v in T.map.items() if v is not None}

    Sj = mapper.subtract_conditional_left(T.H, T.G, map=None, pool=pool)
    Sj = graphs.structure_remove_full_leaves(Sj)
    Sj.cache["hash"] = None

    mj = {k: v for k, v in mj.items() if v in Sj.nodes}

    Tsj = mapper.mapped_type(S0, Sj, mj)

    if len(Sj.nodes) - len(S0.nodes) > splitter.branch_limit:
        # print(f"2WORKER DIFFERENCE FOR A {hash(Sj)} SPLITS {makes_split}")
        return Tsj, shard, hash(shard), matches, makes_split

    if len(Sj.nodes) - len(S0.nodes) < splitter.branch_min:
        # print(f"3WORKER DIFFERENCE FOR A {hash(Sj)} SPLITS {makes_split}")
        return Tsj, shard, hash(shard), matches, makes_split

    if graphs.structure_max_depth(Sj) < splitter.branch_depth_min:
        # print(f"4WORKER DIFFERENCE FOR A {hash(Sj)} SPLITS {makes_split}")
        return Tsj, shard, hash(shard), matches, makes_split

    if splitter.return_matches:
        matches = tuple(
            (
                j
                for j, ai in enumerate(A)
                if mapper.mapper_match(ai, Sj, pool=pool)
            )
        )

        makes_split = len(matches) > 0 and len(matches) < len(A)
    else:
        yes = 0
        no = 0
        matches = tuple()

        for ai in A:
            if mapper.mapper_match(ai, Sj, pool=pool):
                yes = 1
            else:
                no = 1
            if yes and no:
                makes_split = True
                break

    Sj.cache["hash"] = None
    h = hash(Sj)
    # print(f"WORKER DIFFERENCE FOR A {h} SPLITS {makes_split}")
    return (Tsj, shard, hash(Sj), matches, makes_split)


def process_split_intersect(pack, T):
    S0 = process_split_ctx.S0
    splitter = process_split_ctx.splitter
    A = process_split_ctx.A
    pool = process_split_ctx.pool

    try:
        bit, mapping = pack[0]
    except:
        pack = (pack,)

    makes_split = False
    shard = T.G
    shard.cache["hash"] = None

    matches = tuple()

    mj = {v: k for k, v in T.map.items() if v is not None}

    Sj = mapper.intersection_conditional(T.H, T.G, map=mj, pool=pool)
    Sj.cache["hash"] = None

    Tsj = mapper.mapped_type(S0, Sj, mj)

    if len(Sj.nodes) - len(S0.nodes) > splitter.branch_limit:
        # print(f"2WORKER INTERSECT FOR A {hash(Sj)} SPLITS {makes_split}")
        return Tsj, shard, hash(shard), matches, makes_split

    if len(Sj.nodes) - len(S0.nodes) < splitter.branch_min:
        # print(f"3WORKER INTERSECT FOR A {hash(Sj)} SPLITS {makes_split}")
        return Tsj, shard, hash(shard), matches, makes_split

    if graphs.structure_max_depth(Sj) < splitter.branch_depth_min:
        # print(f"4WORKER INTERSECT FOR A {hash(Sj)} SPLITS {makes_split}")
        return Tsj, shard, hash(shard), matches, makes_split

    if splitter.return_matches:
        matches = tuple(
            (
                j
                for j, ai in enumerate(A)
                if mapper.mapper_match(ai, Sj, pool=pool)
            )
        )

        makes_split = len(matches) > 0 and len(matches) < len(A)
    else:
        yes = 0
        no = 0
        matches = tuple()

        for ai in A:
            # T = mapper.map_to(
            #     ai, Sj, strict=True, equality=False, add_nodes=1, fill=True
            # )
            if mapper.mapper_match(ai, Sj, pool=pool):
                yes = 1
            else:
                no = 1
            if yes and no:
                makes_split = True
                break

    h = hash(Sj)
    # print(f"WORKER INTERSECT FOR A {h} SPLITS {makes_split}")
    return (Tsj, shard, h, matches, makes_split)


def process_split_difference_distributed(pack, T, S0, splitter, A, selections, icd):

    try:
        bit, mapping = pack[0]
    except Exception:
        pack = (pack,)

    makes_split = False
    shard = T.G
    matches = tuple()
    shard.cache["hash"] = None

    mj = {v: k for k, v in T.map.items() if v is not None}

    Sj = mapper.subtract_conditional_left(T.H, T.G, map=None)
    Sj = graphs.structure_remove_full_leaves(Sj)
    Sj.cache["hash"] = None

    mj = {k: v for k, v in mj.items() if v in Sj.nodes}

    Tsj = mapper.mapped_type(S0, Sj, mj)

    Sj.cache["hash"] = None
    h = hash(Sj)

    if len(Sj.nodes) - len(S0.nodes) > splitter.branch_limit:
        # print(f"2WORKER DIFFERENCE FOR A {hash(Sj)} SPLITS {makes_split}")
        return Tsj, shard, h, matches, makes_split

    if len(Sj.nodes) - len(S0.nodes) < splitter.branch_min:
        # print(f"3WORKER DIFFERENCE FOR A {hash(Sj)} SPLITS {makes_split}")
        return Tsj, shard, h, matches, makes_split

    if not graphs.graph_is_valid(Sj):
        return None, shard, h, matches, makes_split

    if graphs.structure_max_depth(Sj) < splitter.branch_depth_min:
        # print(f"4WORKER DIFFERENCE FOR A {hash(Sj)} SPLITS {makes_split}")
        return Tsj, shard, h, matches, makes_split

    if graphs.graph_same(S0, Sj):
        # print(f"4WORKER DIFFERENCE FOR A {hash(Sj)} SPLITS {makes_split}")
        return None, shard, h, matches, makes_split

    matches, makes_split = process_split_matches(Sj, A, selections, icd, False)

    # print(f"WORKER DIFFERENCE FOR A {h} SPLITS {makes_split}")
    return (Tsj, shard, h, matches, makes_split)


def process_split_matches(Sj, A, selections, icd: codecs.intvec_codec, return_matches=True):

    yes = 0
    no = 0
    # matches: List[bool] = list([None] * len(A))
    matches = [None, len(A)]
    makes_split = False

    for i, (idx, sel) in enumerate(selections):
        ai = graphs.graph_as_structure(icd.graph_decode(A[idx]), sel, Sj.topology)

        if mapper.mapper_match(ai, Sj):
            yes = 1
            # matches[i] = True
            if matches[0] is None:
                matches[0] = True
        else:
            no = 1
            # matches[i] = False
            if matches[0] is None:
                matches[0] = False
        if yes and no:
            makes_split = True
            matches[1] = i
            break

    return matches, makes_split


def process_split_intersect_distributed(pack, T, S0, splitter, A, selections, icd: codecs.intvec_codec):

    try:
        bit, mapping = pack[0]
    except Exception:
        pack = (pack,)

    makes_split = False
    shard = T.G
    shard.cache["hash"] = None

    matches = tuple()

    mj = {v: k for k, v in T.map.items() if v is not None}

    Sj = mapper.intersection_conditional(T.H, T.G, map=mj)
    Sj = graphs.structure_remove_full_leaves(Sj)
    Sj.cache["hash"] = None
    h = hash(Sj)

    Tsj = mapper.mapped_type(S0, Sj, mj)

    if len(Sj.nodes) - len(S0.nodes) > splitter.branch_limit:
        # print(f"2WORKER INTERSECT FOR A {hash(Sj)} SPLITS {makes_split}")
        return Tsj, shard, h, matches, makes_split

    if len(Sj.nodes) - len(S0.nodes) < splitter.branch_min:
        # print(f"3WORKER INTERSECT FOR A {hash(Sj)} SPLITS {makes_split}")
        return Tsj, shard, h, matches, makes_split

    if not graphs.graph_is_valid(Sj):
        return None, shard, hash(shard), matches, makes_split

    if graphs.structure_max_depth(Sj) < splitter.branch_depth_min:
        # print(f"4WORKER INTERSECT FOR A {hash(Sj)} SPLITS {makes_split}")
        return Tsj, shard, h, matches, makes_split

    if graphs.graph_same(S0, Sj):
        return None, shard, h, matches, makes_split

    matches, makes_split = process_split_matches(Sj, A, selections, icd, False)

    # print(f"WORKER INTERSECT FOR A {h} SPLITS {makes_split}")
    return (Tsj, shard, h, matches, makes_split)


def split_subgraphs(
    topology: structure_topology,
    splitter: smarts_splitter_config,
    S0: subgraph,
    G: Sequence[graphs.graph],
    selections,
    Q=None,
    verbose=False,
    debug=False,
) -> List[mapper.mapped_type]:
    # verbose = True
    # debug = True

    cdcs = codec_native.primitive_codecs_get()
    gcd = codec_native.graph_codec_native(
        cdcs,
        list(codec_native.primitive_codecs_get_atom()),
        list(codec_native.primitive_codecs_get_bond())
    )
    icd = codecs.intvec_codec(
        gcd.primitive_codecs,
        gcd.atom_primitives,
        gcd.bond_primitives
    )

    # A = tuple((structure(ai.nodes, ai.edges, ai.select, topology) for ai in A))
    S0 = structure(S0.nodes, S0.edges, S0.select, topology)

    procs = configs.processors

    matched = []
    shards = []
    S = []

    print(datetime.datetime.now(), "Generating splits")
    single_bits = split_single_bits(topology, splitter, S0, G, selections, icd, Q=Q)
    print(datetime.datetime.now(), f"Generated {len(single_bits)} splits")
    # single_bits.extend([(, m.copy()) for (b,m) in single_bits])

    max_bits = splitter.bit_search_limit
    min_bits = min(splitter.bit_search_min, len(single_bits))
    uptobits = min(max_bits + 1, len(single_bits))
    hits = 0

    if len(single_bits) == 0:
        return S, shards, matched

    output = False

    single_bits_red = list()
    single_bits_gra = list()
    single_bits_sma = list()
    if True:
        seen = set()

        smarts = [gcd.smarts_encode(ai) for ai in A]

        for b, m in single_bits:
            g = gcd.smarts_encode(b)
            # print("BIT", g)
            if b not in single_bits_gra:  # and g not in single_bits_sma:
                single_bits_red.append((b, m))
                single_bits_sma.append(g)
                single_bits_gra.append(b)

        for b, sma in zip(single_bits_gra, single_bits_sma):
            # Should be all unique bits
            print("BIT", sma)

        single_bits = single_bits_red

    process_split_ctx.A = A
    process_split_ctx.splitter = splitter
    process_split_ctx.S0 = S0

    for i in range(min_bits, uptobits):
        visited = set()
        Bn = math.factorial(len(single_bits)) // (
            math.factorial(len(single_bits) - i) * math.factorial(i)
        )
        procs = min(Bn, procs)
        clen = max(procs, procs * ((Bn // procs) + bool(Bn % procs)) // 10)
        chunksize = min(1, clen // procs * 1000 // len(A))
        if len(A) * (clen // procs) < 10:
            chunksize = clen // procs
        else:
            chunksize = max(1, (clen // procs // len(A)))
            if chunksize * procs > clen:
                chunksize = max(clen // procs, 1)
        if Bn // clen < 10:
            clen = Bn
            chunksize = clen // procs
        # chunksize = 1
        print(
            f"{datetime.datetime.now()} Splits N: {Bn} Chunks N: {clen} Chunk Sz: {chunksize}"
        )

        all_completed = 0
        with multiprocessing.pool.Pool(
            processes=procs, maxtasksperchild=None
        ) as pool:
            print()
            for ci, chunk in enumerate(
                arrays.batched(itertools.combinations(single_bits, i), clen)
            ):
                work = None
                process_split_ctx.pool = pool
                # print(datetime.datetime.now(), "Processing batch", ci)

                work = [
                    pool.apply_async(process_split, (unit,)) for unit in chunk
                ]
                print(f"{datetime.datetime.now()} Submitted {len(work)}")
                completed = set()
                updates = set()
                while len(work) != len(completed):
                    for jj, splits in enumerate(work, clen * ci + 1):
                        if jj in completed:
                            continue
                        if splits.ready():
                            completed.add(jj)
                        else:
                            continue

                        splits = splits.get()
                        all_completed += 1

                        for j, unit in enumerate(splits):
                            matches = None
                            shard = None
                            Tsj = None

                            if unit is not None:
                                (
                                    Tsj,
                                    shard,
                                    hashshard,
                                    matches,
                                    makes_split,
                                ) = unit
                            else:
                                continue

                            Sj = Tsj.H

                            # if jj + j == clen * ci + len(work) and j == len(splits):
                            progress = int(((all_completed) / Bn * 10))
                            report_number = int(
                                progress / 10 * Bn * len(splits)
                            )
                            if (verbose and debug) or (
                                report_number not in updates
                            ):
                                updates.add(report_number)
                                print(
                                    datetime.datetime.now(),
                                    f"Searching atoms={len(shard.nodes)}"
                                    f" data={len(A)}"
                                    f" bit_depth={i}/{uptobits-1}"
                                    f" b_j={report_number}/{Bn*len(splits)}"
                                    f" hits={hits}            ",
                                    end="\n",
                                )
                                output = True

                            if verbose and debug:
                                print(
                                    "S0 =>",
                                    gcd.smarts_encode(Tsj.G),
                                    "\nSj =>",
                                    gcd.smarts_encode(Sj),
                                    "\nbj =>",
                                    gcd.smarts_encode(shard),
                                    makes_split,
                                    hashshard,
                                )
                            # Sj.cache['hash'] = None
                            # S0.cache['hash'] = None
                            # hashshard = hash(Sj)
                            if hashshard in visited:
                                # print(f"HASH DUP: {hashshard}")
                                continue
                            else:
                                # print(f"HASH NEW: {hashshard}")
                                visited.add(hashshard)
                            if S0 == Sj:
                                # print("COPY")
                                continue
                            if (
                                len(Sj.nodes) - len(S0.nodes)
                                > splitter.branch_limit
                            ):
                                # print("TOO LARGE")
                                continue

                            if verbose and debug:
                                for i, sma in enumerate(smarts):
                                    if splitter.return_matches and matches:
                                        is_match = i in matches
                                    else:
                                        is_match = mapper.mapper_match(A[i], Sj)
                                    print("     ", f"{str(is_match):6s}", sma)
                                print()

                            unique_split = matches not in matched
                            if makes_split and (
                                (not splitter.unique) or unique_split
                            ):
                                hits += 1
                                S.append(Tsj)
                                shards.append(shard)
                                matched.append(matches)
                                if (
                                    splitter.max_splits > 0
                                    and hits > splitter.max_splits
                                ):
                                    break
                            else:
                                pass
                                # print(f"SPLITS? {makes_split}")
                                # print(f"UNIQUE? {unique_split}")

        if splitter.max_splits > 0 and hits > splitter.max_splits:
            break

    if not splitter.return_matches and hits > 0:
        print(
            datetime.datetime.now(), f"Calculating partitions for hits={hits}"
        )
        with multiprocessing.pool.Pool(
            processes=min(procs, len(S)), maxtasksperchild=None
        ) as pool:
            work = [pool.apply_async(process_split_matches, (T.H,)) for T in S]
            for i, unit in enumerate(work):
                matches = unit.get()
                matched[i] = matches

    print(
        datetime.datetime.now(),
        f"Searching atoms done; data={len(A)} hits={hits}",
        end="\n",
    )
    output = True
    if output:
        print()

    process_split_ctx.A = None
    process_split_ctx.splitter = None
    process_split_ctx.S0 = None
    process_split_ctx.pool = None

    return S, shards, matched


def split_subgraphs_distributed(
    topology: structure_topology,
    splitter: smarts_splitter_config,
    S0: subgraph,
    G: Sequence[graphs.graph],
    selections,
    wq: compute.workqueue_local,
    icd: codecs.intvec_codec,
    Q=None,
    verbose=False,
    debug=False,
) -> List[mapper.mapped_type]:
    # verbose = True
    # debug = True

    # A = tuple((structure(ai.nodes, ai.edges, ai.select, topology) for ai in A))
    S0 = structure(S0.nodes, S0.edges, S0.select, topology)

    procs = configs.processors

    matched = []
    shards = []
    S = []

    print(datetime.datetime.now(), "Generating splits")
    single_bits = split_single_bits(topology, splitter, S0, G, selections, icd, Q=Q)
    print(datetime.datetime.now(), f"Generated {len(single_bits)} splits")

    max_bits = splitter.bit_search_limit
    min_bits = min(splitter.bit_search_min, len(single_bits))
    uptobits = min(max_bits + 1, len(single_bits))
    hits = 0

    if len(single_bits) == 0:
        return S, shards, matched

    closer = None
    output = False

    single_bits_red = list()
    single_bits_gra = list()
    single_bits_sma = list()

    if True:

        cdcs = codec_native.primitive_codecs_get()
        gcd = codec_native.graph_codec_native(
            cdcs,
            list(codec_native.primitive_codecs_get_atom()),
            list(codec_native.primitive_codecs_get_bond())
        )

        if verbose and debug:
            print(f"{datetime.datetime.now()} Decoding graphs")
            decoded = {}
            for (i, sel) in selections:
                if i not in decoded:
                    decoded[i] = icd.graph_decode(G[i])

            smarts = [gcd.smarts_encode(graphs.graph_as_subgraph(decoded[idx], sel)) for idx, sel in selections]
            del decoded

        for b, m in single_bits:
            g = gcd.smarts_encode(b)
            # print("BIT", g)
            if b not in single_bits_gra:  # and g not in single_bits_sma:
                single_bits_red.append((b, m))
                single_bits_sma.append(g)
                single_bits_gra.append(b)

        for b, sma in zip(single_bits_gra, single_bits_sma):
            # Should be all unique bits
            print("BIT", sma)

        single_bits = single_bits_red

    # this should be encoded and whatnot before start
    # shm = shm_split_subgraphs(splitter, S0, A)

    shm = {
        "splitter": splitter,
        "S0": S0,
        "G": G,
        "selections": selections,
        "icd": icd
    }

    # we need 1 for this main process, and the other for the workspace server
    nproc = max(1, configs.processors - 1)


    offsets = {}
    iterable = {}
    idx_offset = 0
    n_ops = int(splitter.split_specific) + int(splitter.split_general)
    print(f"{datetime.datetime.now()} Building tasks")
    for i in range(min_bits, uptobits):
        iterable.update(
            {
                n_ops * idx: unit
                for idx, unit in enumerate(
                    itertools.combinations(single_bits, i),
                    idx_offset,
                )
            }
        )
        if n_ops == 2:
            iterable.update({idx + 1: unit for idx, unit in iterable.items()})
        idx_offset += len(iterable)
        offsets[i] = idx_offset

    addr = ("", 0)
    if len(iterable) <= nproc:
        addr = ('127.0.0.1', 0)
        nproc = len(iterable)

    ws = compute.workqueue_new_workspace(wq, address=addr, nproc=nproc, shm=shm)

    Bn = 0
    completed = set()
    all_completed = 0
    results = {}
    updates = set()
    visited = set()
    sma_visited = set()
    k = 0
    n = 0
    chunksize = nproc
    while len(iterable):
        Bn = n_ops * (
            math.factorial(len(single_bits))
            // (math.factorial(len(single_bits) - i) * math.factorial(i))
        )

        # unfinished = {
        #     idx: unit
        #     for idx, unit in iterable.items()
        #     if idx not in completed
        # }
        # clear the queue
        compute.workspace_submit_and_flush(ws, None, {})
        for batch in arrays.batched(list(iterable.items()), 1000000):
            ids = set()
            for chunk in arrays.batched(batch, 10):
                tasks = {}
                for idx, unit in chunk:
                    ids.add(idx)
                    if idx % n_ops:
                        if splitter.split_specific:
                            tasks[idx] = (
                                process_split_specific_distributed,
                                (unit,),
                                {},
                            )
                        elif splitter.split_general:
                            tasks[idx] = (
                                process_split_general_distributed,
                                (unit,),
                                {},
                            )
                    else:
                        if splitter.split_general:
                            tasks[idx] = (
                                process_split_general_distributed,
                                (unit,),
                                {},
                            )
                        elif splitter.split_specific:
                            tasks[idx] = (
                                process_split_specific_distributed,
                                (unit,),
                                {},
                            )

                compute.workspace_local_submit(ws, tasks)

            k = 0
            these_results = compute.workspace_flush(
                ws, ids, timeout=.1
            )
            while these_results:
                for i in range(min_bits, uptobits):
                    n = offsets[i]
                    n0 = 0
                    if i > min_bits:
                        n0 = offsets[i-1]
                    for idx, splits in sorted(
                        these_results.items(), key=lambda x: x[0]
                    ):

                        if idx >= n or idx < n0:
                            continue

                        if idx in iterable:
                            iterable.pop(idx)
                            these_results.pop(idx)
                            all_completed += 1

                        for j, unit in enumerate(splits):
                            matches = None
                            shard = None
                            Tsj = None

                            if unit is not None:
                                Tsj, shard, hashshard, matches, makes_split = unit
                            else:
                                # print(f"unit IS NONE")
                                continue

                            # if jj + j == clen * ci + len(work) and j == len(splits):
                            progress = int(len([x for x in completed if x < n]) / Bn * 10)
                            report_number = progress
                            if (verbose and debug) and (report_number not in updates):
                                updates.add(report_number)
                                print(
                                    datetime.datetime.now(),
                                    f"Searching atoms={len(shard.nodes)}"
                                    f" data={len(selections)}"
                                    f" bit_depth={i}/{uptobits-1}"
                                    f" b_j={all_completed}/{Bn}"
                                    f" hits={hits}            ",
                                    end="\n",
                                )

                            if Tsj is None:
                                # print(f"Tsj IS NONE")
                                continue

                            Sj = Tsj.H

                            if verbose and debug:
                                print(
                                    "S0 =>",
                                    gcd.smarts_encode(Tsj.G),
                                    "\nSj =>",
                                    gcd.smarts_encode(Sj),
                                    "\nbj =>",
                                    gcd.smarts_encode(shard),
                                    makes_split,
                                    hashshard,
                                )

                            if matches is None or len(matches) == 0:
                                # print(f"MATCHES IS {matches}")
                                # print(f"CND SPLITS=O {sma}")
                                continue

                            if hashshard in visited:
                                # print(f"HASH DUP: {hashshard}")
                                # print(f"{idx+1:5d} CND DUPLICATE {sma}")
                                continue
                            else:
                                # print(f"HASH NEW: {hashshard}")
                                visited.add(hashshard)

                            sma = gcd.smarts_encode(Sj)

                            if sma in sma_visited:
                                # print(f"{idx+1:5d} CND DUPLICATE {sma}")
                                continue
                            else:
                                sma_visited.add(sma)

                            _matches = list([matches[0]] * matches[1])

                            if matches[1] < len(selections):
                                _matches.append(not matches[0])

                            matches = _matches

                            if len(matches) < len(selections):
                                matches.extend([None] * (len(selections) - len(matches)))

                            matches = tuple(matches)
                            if False and verbose and debug:
                                for ii, sma in enumerate(smarts):
                                    is_match = mapper.mapper_match(A[ii], Sj)
                                    print("     ", f"{str(is_match):6s}", sma)
                                print()

                            makes_split_str = "Y" if makes_split else "N"
                            print(f"{idx+1:5d} CND SPLITS={makes_split_str}  {sma}")
                            if makes_split: #and ((not splitter.unique) or unique_split):
                                hits += 1
                                S.append(Tsj)
                                shards.append(shard)
                                matched.append(matches)
                                if (
                                    splitter.max_splits > 0
                                    and hits > splitter.max_splits
                                ):
                                    break
                    if splitter.max_splits > 0 and hits > splitter.max_splits:
                        break
                k = all_completed
                # print(f"Progress: {k/n*100:5.2f}%  {k:8d}/{n}")
                if splitter.max_splits > 0 and hits > splitter.max_splits:
                    break
    if n > 0:
        print(f"Finished: {k/n*100:5.2f}%  {k:8d}/{n}")
        if len(completed) < len(iterable) and not (splitter.max_splits > 0 and hits > splitter.max_splits):
            breakpoint()
            print("something wrong...")

    compute.workqueue_remove_workspace(wq, ws)
    print("Closing workspace")
    ws.close()
    ws = None

    n = 0
    if splitter.return_matches and hits > 0:
        print(
            datetime.datetime.now(), f"Calculating partitions for hits={hits}"
        )

        results.clear()
        matched.clear()
        completed = set()

        code = arrays.find_unsigned_typecode_min(len(selections))

        chunksize = 1000  # should be about 10 seconds for torsions at d=0

        for idx in enumerate(S):
            matched.append(array.array(code))

        iterable = {
            (idx, j): ((tsj.H, array.array(code, chunk)), {})
            for idx, tsj in enumerate(S, 0)
            for j, chunk in enumerate(arrays.batched(range(len(selections)), chunksize))
        }

        ip = ""
        port = 0
        if len(iterable) <= nproc:
            ip = "127.0.0.1"
            nproc = len(iterable)

        ws = compute.workqueue_new_workspace(
            wq, (ip, port), shm=shm, nproc=nproc
        )

        j = 0
        n = len(iterable)
        print(f"Submitting {n} packets of work")

        results = compute.workspace_submit_and_flush(
            ws,
            process_split_matches_distributed,
            iterable,
            chunksize=100,
            batchsize=50000,
            verbose=True
        )

        for (idx, i), matches in results.items():
            if (idx, i) not in completed:
                completed.add((idx, i))
                j = len(completed)
                matched[idx].extend(array.array(code, matches))
        # while len(completed) < n:
        #     unfinished_all = [
        #         (key, unit)
        #         for key, unit in iterable.items()
        #         if key not in completed
        #     ]
        #     # break into large chunks or else memory becomes an issue since
        #     # the jobs are sections of matches. The sections are condensed here
        #     # so we need to take a breath every once in awhile
        #     compute_chunksize = min(len(unfinished_all), 50000)
        #     unfinished_chunks = arrays.batched(
        #         unfinished_all, compute_chunksize
        #     )

        #     for unfinished in unfinished_chunks:
        #         # this chunk is for how many are packed into a single task for
        #         # a worker. Try to make this about the number of processors per
        #         # remote compute worker. This will make a single queue get to
        #         # saturate the worker. However, if there is a lot of work and
        #         # the workers finish too fast it came gum up the queues.
        #         for chunk in arrays.batched(unfinished, 100):
        #             tasks = {}
        #             for (idx, i), (T, x) in chunk:
        #                 tasks[(idx, i)] = (
        #                     process_split_matches_distributed,
        #                     (T.H, x),
        #                     {},
        #                 )
        #             compute.workspace_local_submit(ws, tasks)
        #         # print("Flushing")
        #         new_results = compute.workspace_flush(
        #             ws, set((key for key, _ in unfinished)), timeout=0.0
        #         )
        #         new_results = compute.workspace_submit_and_flush(
        #             ws,
        #             process_split_matches_distributed,
        #             unfinished_all,
        #             chunksize=100,
        #             batchsize=50000,
        #             verbose=True
        #         )
        #         # print("Flushing Done")
        #         # completed.update(new_results)
        #         # print("Collecting results")
        #         for (idx, i), matches in new_results.items():
        #             if (idx, i) not in completed:
        #                 completed.add((idx, i))
        #                 j = len(completed)
        #                 matched[idx].extend(array.array(code, matches))
        #         j = len(completed)
        #         print(f"Chunk: {j/n*100:5.2f}%  {j:8d}/{n}", end="\n")
        # j = len(completed)
        # print(
        #     f"Finished: {j/n*100:5.2f}%  {j:8d}/{n}"
        # )
        # if j < n:
        #     breakpoint()
        #     print("something is wrong")

        for i, row in enumerate(matched):
            matched[i] = array.array(code, sorted(set(row)))
        compute.workqueue_remove_workspace(wq, ws)

        print("Closing workspace")
        ws.close()
        ws = None

        cnd_lines = {}

        keep = []
        all_S = S
        all_shards = shards
        all_matched = matched
        kept = {}

        if splitter.unique:
            indices = set(range(len(selections)))
            keep = {}
            # reverse because usually simplist patterns are first and
            # we want to prioritize those for the same match
            M = 0
            N = 0
            for i, (Sj, bj, match) in enumerate(reversed(list(zip(S, shards, matched))), 1):

                unmatch = sorted(indices.difference(match))
                N = len(unmatch)
                M = len(match)
                key = tuple(match)
                if key not in kept:
                    kept[key] = []
                if splitter.unique_complements:
                    if N < M:
                        key = unmatch
                key = tuple(key)

                if key not in kept:
                    kept[key] = []

                if key in keep:
                    # prefer the one that matched the least
                    if splitter.unique_complements_prefer_min:
                        if len(keep[key][2]) > M:
                            keep[key] = Sj, bj, match
                            kept[key].insert(0, i)
                            # print("insert", i, key)
                        else:
                            kept[key].append(i)
                            # print("append", i, key)
                    elif len(keep[key][2]) < M:
                            keep[key] = Sj, bj, match
                            kept[key].insert(0, i)
                            # print("insert", i, key)
                    else:
                        kept[key].append(i)
                        # print("append", i, key)
                else:
                    keep[key] = Sj, bj, match
                    kept[key].insert(0, i)
                    # print("insert", i, key)

            S, shards, matched = [], [], []

            for s, b, m in keep.values():
                S.append(s)
                shards.append(b)
                matched.append(match)
            unique_hits = len(keep)
            print(
                datetime.datetime.now(),
                f"Unique hits {unique_hits}/{hits}",
                end="\n",
            )
            hits = unique_hits
        if kept is not None:
            sj = [*reversed(all_S)]
            for key, ilst in kept.items():

                if not ilst:
                    continue

                i = ilst[0]
                match = [*reversed(all_matched)][i-1]
                unmatch = sorted(indices.difference(match))
                N = len(unmatch)
                M = len(match)
                sma = gcd.smarts_encode(sj[i-1].H)
                line = f"{i:5d} HIT S0= {N:<5d} -> Sj= {M:<5d} {sma}"
                print(line)
                for i in ilst[1:]:
                    sma = gcd.smarts_encode(sj[i-1].H)
                    line = f"{i:5d}     DUP {sma}"
                    print(line)
        else:
            for i, (Sj, bj, match) in enumerate(
                reversed(list(zip(all_S, all_shards, all_matched))),
                1
            ):
                key = match
                unmatch = sorted(indices.difference(match))
                N = len(unmatch)
                M = len(match)
                is_uniq = "X"
                sma = gcd.smarts_encode(Sj.H)
                line = f"{i:5d} HIT S0= {N:<5d} -> Sj= {M:<5d} UNIQUE={is_uniq} {sma}"
                print(line)

    print(
        datetime.datetime.now(),
        f"Searching atoms done; data={len(selections)} hits={hits}",
        end="\n",
    )

    return S, shards, matched


def split_single_bits(
    topology: structure_topology,
    splitter: smarts_splitter_config,
    S0: subgraph,
    G: Sequence[graphs.graph],
    selections,
    icd,
    Q=None,
):
    # A = [structure(ai.nodes, ai.edges, ai.select, topology) for ai in A]
    S0 = structure(S0.nodes, S0.edges, S0.select, topology)
    # Q = graphs.structure_copy(A[0])

    if Q is None:
        Q = mapper.union_list_parallel(
            G, selections, topology, reference=S0, max_depth=graphs.structure_max_depth(S0), icd=icd
        )

    relabel = {}
    for i, j in enumerate(
        [x for x in Q.select] + [x for x in Q.nodes if x not in Q.select], 1
    ):
        relabel[j] = i

    bes = graphs.structure_relabel_nodes(Q, relabel)

    T = mapper.map_to(bes, S0, add_nodes=2, fill=0)
    bes = T.G
    M = T.map
    bes.select = tuple((relabel[Q.select[i]] for i in topology.primary))
    extend = set(M).difference(bes.select)
    bes.select = tuple(
        [relabel[Q.select[i]] for i in topology.primary] + list(extend)
    )
    print(datetime.datetime.now(), "Generating single splits")

    # if both are set then we are looking for combos of specific and 
    # general, e.g. [#6H2!X3] so we need both bits here since we
    # combine them together before taking a specific/general split
    iter_inverse = splitter.split_general and splitter.split_specific

    primitives = None
    if hasattr(splitter, "primitives"):
        primitives = splitter.primitives

    single_bits = [
        (bit, M)
        for bit in graph_visitors.structure_iter_bits(
            bes, iter_inverse=iter_inverse, skip_ones=True, primitives=primitives
        )
    ]

    N = len(single_bits)

    if splitter.branch_depth_limit > 0 and splitter.branch_limit > 0:
        print(datetime.datetime.now(), "Generating branched splits")
        branch_bits = split_branches(
            S0, G, selections, splitter.branch_depth_limit, splitter.branch_limit, primitives=splitter.primitives, icd=icd
        )
        for i, b in enumerate(branch_bits, 1):
            if (splitter.branch_limit is not None) and len(b[0].select) - len(
                S0.select
            ) > splitter.branch_limit:
                continue
            single_bits.append(b)
            # b = (graphs.structure_invert(b[0]), M)
            # single_bits.append(b)

    return single_bits


def split_structures_distributed(
    splitter: smarts_splitter_config,
    S0: structure,
    G: Sequence[graphs.graph],
    selections,
    wq: compute.workqueue_local,
    icd,
    Q=None,
) -> split_return_type:
    topology = S0.topology
    # sg_list = tuple([subgraph(ai.nodes, ai.edges, ai.select) for ai in A])

    Tsj, shards, matched = split_subgraphs_distributed(
        topology,
        splitter,
        subgraph(S0.nodes, S0.edges, S0.select),
        G,
        selections,
        wq,
        icd,
        Q=Q,
    )

    result = split_return_type(
        tuple((subgraph(T.H.nodes, T.H.edges, T.H.select) for T in Tsj)),
        tuple((subgraph(ai.nodes, ai.edges, ai.select) for ai in shards)),
        tuple(matched),
        tuple((i for i in range(len(selections)) if i not in matched)),
        selections,
        topology,
    )

    return result


def split_structures(
    splitter: smarts_splitter_config,
    S0: structure,
    G: Sequence[structure],
    selections,
    Q=None,
) -> split_return_type:
    topology = S0.topology
    sg_list = tuple([subgraph(ai.nodes, ai.edges, ai.select) for ai in A])

    Tsj, shards, matched = split_subgraphs(
        topology,
        splitter,
        subgraph(S0.nodes, S0.edges, S0.select),
        G,
        selections,
        Q=Q,
    )

    result = split_return_type(
        tuple((subgraph(T.H.nodes, T.H.edges, T.H.select) for T in Tsj)),
        tuple((subgraph(ai.nodes, ai.edges, ai.select) for ai in shards)),
        tuple(matched),
        tuple((i for i in range(len(A)) if i not in matched)),
        selections,
        topology,
    )

    return result
