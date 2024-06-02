import io
from besmarts.codecs import codec_rdkit, codec_native
from besmarts.core import graphs
from besmarts.core import configs
from besmarts.core import mapper
from besmarts.core import topology
from besmarts.core import splits
from besmarts.core import codecs
from besmarts.core import graph_visitors

import unittest

class test_branch_split(unittest.TestCase):

    def setUp(self):
        pass

    def test_atom_branch_one(self):
        smiles = {}
        s = "[#6:1]([#1:2])([#2:3])([#3:4])[#4:5]"
        smiles[s] = [(1,)]
        branches = 1
        depth = 1
        runner(self, smiles, branches, depth)

    def test_atom_branch_two(self):
        smiles = {}
        s = "[#6:1]([#1:2])([#2:3])([#3:4])[#4:5]"
        smiles[s] = [(1,)]
        branches = 2
        depth = 1
        runner(self, smiles, branches, depth)

def runner(self, smiles, branches, depth):

    #### Next
    # print("Next")

    torsions = []
    # 4 and 7
    gcd = codec_rdkit.graph_codec_rdkit()
    icd = codecs.intvec_codec(
        gcd.primitive_codecs,
        gcd.atom_primitives,
        gcd.bond_primitives
    )
    S0 = gcd.smarts_decode("[*:1]")
    topo = topology.atom

    G = {}
    selections = []

    for i, (smi, keep) in enumerate(smiles.items()):
        G[i] = icd.graph_encode(gcd.smiles_decode(smi))
        for k in keep:
            selections.append((i, k))

    cfg = configs.smarts_extender_config(depth, depth, True)
    mapper.mapper_smarts_extend(cfg, torsions)

    config = configs.mapper_config(1, False, "high")
    t_union = mapper.union_list_parallel(G, selections, topo, config, icd=icd)

    t_union.select = tuple(
        list(t_union.select)
        + [x for x in t_union.nodes if x not in t_union.select]
    )

    # print("Union:")
    # print(t_union.select)
    # print(gcd.smarts_encode(t_union))


    B = graphs.structure_branch(
        t_union, {k: k for k in t_union.select[:1]}, branches, depth, verbose=False
    )
    B = list(B)

    # print(f"There are {len(B)} branches:")

    # for b in B:
    #     print(
    #         # len(b.nodes),
    #         # graphs.structure_max_depth(b),
    #         # sorted(b.nodes.keys()),
    #         gcd.smarts_encode(b),
    #     )

    Q = mapper.union_list(
        [graphs.subgraph_to_structure(x, topo) for x in B], reference=graphs.subgraph_to_structure(S0, topo), max_depth=depth
    )

    # print("Union of branches is", gcd.smarts_encode(Q))

    scfg = configs.smarts_splitter_config(1, 1, 0, branches, depth, depth, False, True, 0, split_general=False)

    bits = list(graph_visitors.structure_iter_bits(
        Q, iter_inverse=False, skip_ones=True
    ))
    self.assertTrue(len(set(bits)) == 4)
    bits = list(graph_visitors.structure_iter_bits(
        Q, iter_inverse=True, skip_ones=True
    ))
    self.assertTrue(len(set(bits)) == 8)

    bits = list(graph_visitors.structure_iter_bits(
        Q, iter_inverse=False, skip_ones=False
    ))
    self.assertTrue(len(set(bits)) == 19)
    bits = list(graph_visitors.structure_iter_bits(
        Q, iter_inverse=True, skip_ones=False
    ))
    self.assertTrue(len(set(bits)) == 38)

    if 0:
        uniq = set()
        for s in bits:
            if s in uniq:
                print("X", graphs.structure_max_depth(s), len(s.nodes), f"{hash(s):+0d}", gcd.smarts_encode(s))
                continue
            else:
                uniq.add(s)
                print("O", graphs.structure_max_depth(s), len(s.nodes), f"{hash(s):+0d}", gcd.smarts_encode(s))
        print(len(bits), len(uniq))



