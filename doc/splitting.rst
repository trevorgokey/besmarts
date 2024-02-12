
Splitting
=========

The goal of SMARTS clustering is to find a SMARTS pattern that discriminates 
a group of SMARTS from another group. As such, we are interested in finding
a SMARTS pattern which partitions a set of SMARTS patterns. Since partitioning
is an expensive problem all by itself. We offer two ways to look for SMARTS
clusters:

    1. Given a set of SMARTS, generate the partitions which minimize the number of SMARTS edits
    2. Given an explicit partition, find satisfactory SMARTS pattern that minimizes the number of partition edits

In the first case, we are interested in all of the "adjacent" SMARTS patterns 
and the partitions that they induce. In the second case, we are seeking a
specific partitioning and want to know a SMARTS pattern that will descriminate 
between the two. Note that in both cases we allow a fuzzy result by trying to 
minimize the number of edits; the SMARTS edits in the first case and the 
partition edits in the second case.

Here is how to do the first case

.. code-block:: python

    from besmarts.core.topology import bond_topology
    from besmarts.core.mapper import mapper_smarts_extend
    from besmarts.core.primitives import primitive_key
    from besmarts.core import graphs
    
    from besmarts.codecs.codec_rdkit import graph_codec_rdkit
    
    from besmarts.splitter.perceive import split_return_type, perceive_split_structure
    from besmarts.splitter.configs import smarts_splitter_config, smarts_extender_config
    
    
    gcd = graph_codec_rdkit(*prims)
    gcd.smiles_config.strip_hydrogen = False
    
    smi = "CC(C)=O"
    mol = pcp.smiles_decode(smi)

    ic_list = [s for s in graphs.graph_to_structure_bonds(mol)]
    
    frags = list(ic_list)
    
    branch_limit = 999
    branch_depth = 0
    bit_depth_min = 1
    bit_depth_max = 1
    
    splitter = smarts_splitter_config(bit_depth_min, bit_depth_max, branch_limit, branch_depth, True, True)
    
    # for this to work, we need to extend our graphs to at least the depth of S0
    extender = smarts_extender_config(branch_depth, branch_depth, True)
    mapper_smarts_extend(extender, frags)
    
    S0 = gcd.smarts_decode("[*:1]~[*:2]")
    S0 = graphs.structure(S0.nodes, S0.edges, (1, 2), bond_topology())
    
    for i, f in enumerate(frags):
        print(i, pcp.smarts_encode(f))
    
    print("matched:", len(frags))
    
    frags = [s for s in graphs.graph_to_structure_bonds(mol)]
    mapper_smarts_extend(extender, frags)
    splits: split_return_type = perceive_split_structure(splitter, S0, list(frags))

And here is how to do the second case:

.. code-block:: python

   
    from besmarts.core.graphs import graph_to_structure_bonds
    from besmarts.core.configs import smarts_extender_config
    from besmarts.core.primitives import primitive_key
    from besmarts.core.topology import bond_topology
    from besmarts.core import mapper
    
    from besmarts.codecs.codec_rdkit import graph_codec_rdkit
    
    from besmarts.splitter.configs import smarts_splitter_config, smarts_perception_config
    from besmarts.splitter import perceive
    
    
    
    # want this to be a little nicer:
    
    prims = (primitive_key.ELEMENT, primitive_key.HYDROGEN), (primitive_key.BOND_ORDER,)
    prims = (None, None)
    pcp = graph_codec_rdkit(*prims)
    pcp.smiles_config.strip_hydrogen = False
    
    ###
    
    spec = smarts_perception_config(
        smarts_splitter_config(1, 1, 0, 99, False, True), smarts_extender_config(0, 1, False)
    )
    spec.extender.depth_max = 0
    spec.extender.depth_min = 0
    
    
    ###
    
    
    ###
    smi = "CC(C)(O1)C[C@@H](O)[C@@]1(O2)[C@@H](C)[C@@H]3CC=C4[C@]3(C2)C(=O)C[C@H]5[C@H]4CC[C@@H](C6)[C@]5(C)Cc(n7)c6nc(C[C@@]89(C))c7C[C@@H]8CC[C@@H]%10[C@@H]9C[C@@H](O)[C@@]%11(C)C%10=C[C@H](O%12)[C@]%11(O)[C@H](C)[C@]%12(O%13)[C@H](O)C[C@@]%13(C)CO"
    # smi = "CCCCC"
    beg = pcp.smiles_decode(smi)
    smi = "CCC(C)=O"
    beg2 = pcp.smiles_decode(smi)
    
    ###
    
    sg_list = [s for mol in [beg, beg2] for s in graph_to_structure_bonds(mol)]
    # mapper.mapper_smarts_extend(spec.extender, sg_list)
    # sg_list = [graphs.structure_remove_unselected(s) for s in graph_to_structure_bonds(beg)]
    
    spec.extender.depth_max = 3
    spec.extender.depth_min = 0
    
    ###
    topo = bond_topology()
    
    
    ###
    # matches = (13, 17, 29, 30, 31, 32, 36, 49, 71)
    matches = (0, 1, 2)
    matches = set(
        (
            0,
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21,
            22,
            23,
            24,
            25,
            26,
            27,
            28,
            29,
            30,
            31,
            32,
            33,
            34,
            35,
            36,
            37,
            38,
            39,
            40,
            41,
            42,
            43,
            44,
            45,
            46,
            47,
            48,
            49,
            50,
            51,
            52,
            53,
            54,
            55,
            56,
            57,
            58,
            59,
            60,
            61,
            62,
            63,
            64,
            65,
            66,
            67,
            68,
            69,
            70,
            71,
            72,
            73,
            74,
            75,
            76,
            77,
            78,
            79,
            80,
            81,
            82,
            83,
            84,
            85,
            86,
            87,
            88,
            89,
            90,
            91,
            92,
            93,
            94,
            95,
            96,
            97,
            98,
            99,
            100,
            101,
            102,
            103,
            104,
            105,
            106,
            107,
            108,
            109,
            110,
            111,
            112,
            113,
            114,
            115,
            116,
            117,
            118,
            119,
            120,
            121,
            122,
            123,
            124,
            125,
            126,
            127,
            128,
            129,
            130,
            131,
            132,
            133,
            134,
            135,
            136,
            137,
            138,
            139,
            140,
            141,
            142,
            143,
            144,
            145,
            146,
            147,
            148,
            149,
            150,
        )
    )
    
    matches = [
        2,
        5,
        7,
        17,
        45,
        51,
        53,
        57,
        59,
        64,
        65,
        66,
        76,
        77,
        87,
        131,
        137,
        143,
        151,
    ]
    matches = [2,161]
    matches = [80]
    
    for i in range(len(sg_list)):
        if i not in matches:
            print(i, pcp.smarts_encode(sg_list[i]))
    for i in matches:
        print(i, "->", pcp.smarts_encode(sg_list[i]))
    # matches = tuple(range(10))
    # matches = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 25, 26, 27, 29, 30, 31, 32, 34, 35, 36, 38, 39, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 65, 67, 68, 69, 71, 72, 73, 74, 75, 76, 77, 84, 85, 86, 88, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 112, 113, 114, 115, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 135, 136, 138, 142, 144, 145, 149, 150)
    ret = perceive.perceive_shard_of_partition(topo, spec, sg_list, matches, pcp, 9999)
    shards = ret.value
    
    ###
    removeA = shards[2]
    addA = shards[3]
    nummoves = len(removeA) + len(addA)
    verbose = True
    shard = shards[0]
    matches = [x for x in range(len(sg_list)) if x not in removeA and (x in matches or x in addA)]
    if shard is not None:
        print(f"Matches only the input with {nummoves} swaps:", pcp.smarts_encode(shard))
        if verbose and (removeA or addA):
            print("RemoveA", removeA)
            print("AddA", addA)
            for i in range(len(sg_list)):
                if i not in matches:
                    print(i, pcp.smarts_encode(sg_list[i]))
            for i in range(len(sg_list)):
                if i in matches:
                    print(i, "->", pcp.smarts_encode(sg_list[i]))
    
    shard = shards[1]
    if shard is not None:
        print(f"Matches the input complement with {nummoves} swaps:", pcp.smarts_encode(shard))
        if verbose and (removeA or addA):
            print("RemoveA", removeA)
            print("AddA", addA)
            for i in range(len(sg_list)):
                if i in matches:
                    print(i, pcp.smarts_encode(sg_list[i]))
            for i in range(len(sg_list)):
                if i not in matches:
                    print(i, "->", pcp.smarts_encode(sg_list[i]))

