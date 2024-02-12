import io
from besmarts.codecs import codec_rdkit, codec_native
from besmarts.core import graphs
from besmarts.core import configs
from besmarts.core import mapper
from besmarts.core import topology
from besmarts.core import splits


def build_input():
    inp = """#GRAPH
#ATOM element hydrogen connectivity_total connectivity_ring ring_smallest aromatic formal_charge
#BOND bond_ring bond_order
-21 -21   2   1   2   1   1   1   1
-10 -10  64   2   8   4  16   2   1
 -9  -9  64   2   8   4  16   2   1
 -8  -8  64   2   8   4  16   2   1
 -7  -7  64   1   8   4  16   2   1
-12 -12  64   2   8   4  16   2   1
-11 -11  64   2   8   4  16   2   1
-22 -22   2   1   2   1   1   1   1
-23 -23   2   1   2   1   1   1   1
 -6  -6 128   2   8   1   1   1   1
-18 -18   2   1   2   1   1   1   1
 -4  -4 65536   1   8   1   1   1   1
 -5  -5 256   1   2   1   1   1   1
 -3  -3 256   1   4   1   1   1   1
 -2  -2  64   4  16   1   1   1   1
-16 -16   2   1   2   1   1   1   1
-17 -17   2   1   2   1   1   1   1
 -1  -1  64   8  16   1   1   1   1
-13 -13   2   1   2   1   1   1   1
-14 -14   2   1   2   1   1   1   1
-15 -15   2   1   2   1   1   1   1
-19 -19   2   1   2   1   1   1   1
-20 -20   2   1   2   1   1   1   1
 10  21   1   2
  9  10   2  32
 10  11   2  32
  8   9   2  32
  9  20   1   2
  7   8   2  32
  8  19   1   2
  7  12   2  32
  6   7   1   2
 11  12   2  32
 12  23   1   2
 11  22   1   2
  6  18   1   2
  4   6   1   2
  4   5   1   4
  3   4   1   2
  2   3   1   2
  2  16   1   2
  2  17   1   2
  1   2   1   2
  1  13   1   2
  1  14   1   2
  1  15   1   2
"""
    return io.StringIO(inp)


# check depth of 3 on 21-10-11-22 should get me 6. Interesting test because
# we need to ensure that it doesn't go around the ring the long way

# buf_in = build_input()
# g = codec_native.graph_codec_native_read(buf_in)[0]

if 0:
    smiles = "[H:26][c:6]1[c:7]([c:8]2[c:9]([c:10]([c:11]([c:12]([c:13]2[S:14](=[O:15])(=[O:16])[N:17]([H:31])[N:18]([H:32])[H:33])[H:30])[H:29])[H:28])[c:4]([c:5]1[H:25])[N:2]([C:1]([H:19])([H:20])[H:21])[C:3]([H:22])([H:23])[H:24])[H:27]"
    gcd = codec_rdkit.graph_codec_rdkit()
    g = gcd.smiles_decode(smiles)
    # # codec_native.graph_codec_native_write("g.bg", [g] )
    S0 = gcd.smarts_decode("[#1:1]~[#6H1:2](~[H1])~[H1:3](~[!X4:4])~[x2H0]")
    torsions = graphs.graph_to_structure_dihedrals(g)

    torsions = [
        x for x in torsions if x.select in [(25, 5, 6, 26), (29, 11, 12, 30)]
    ]
    assert len(torsions) == 2

    cfg = configs.smarts_extender_config(2, 2, True)
    mapper.mapper_smarts_extend(cfg, torsions)

    assert 2 in torsions[0].select
    assert 14 in torsions[1].select

    config = configs.mapper_config(1, False, "high")
    t_union = mapper.union_list_parallel(torsions, config)
    print("Union:")
    t_union.select = tuple(
        list(t_union.select)
        + [x for x in t_union.nodes if x not in t_union.select]
    )
    print(t_union.select)
    print(gcd.smarts_encode(t_union))

    branches = 2
    depth = 2

    B = graphs.structure_branch(
        t_union, {k: k for k in t_union.select[:4]}, branches, depth
    )
    B = list(B)
    print(f"There are {len(B)} branches")

    for b in B:
        print(
            len(b.nodes),
            graphs.structure_max_depth(b),
            sorted(b.nodes.keys()),
            gcd.smarts_encode(b),
        )

    if 1:
        scfg = configs.smarts_splitter_config(1, 1, 0, branches, depth, depth)
        spl = splits.split_single_bits(topology.dihedral, scfg, S0, torsions)
        for s, _ in spl:
            print(gcd.smarts_encode(s))

def test(smiles, branches, depth):

    if 1:
        #### Next
        print("Next")

        torsions = []
        # 4 and 7
        gcd = codec_rdkit.graph_codec_rdkit()
        S0 = gcd.smarts_decode("[*:1]~[*:2]~[*:3]~[*:4]")

        for smi, keep in smiles.items():
            a = gcd.smiles_decode(smi)
        # # codec_native.graph_codec_native_write("g.bg", [g] )

            a_torsions = graphs.graph_to_structure_dihedrals(a)
            a_torsions = [
                x
                for x in a_torsions
                if x.select in keep
            ]

            torsions.extend(a_torsions)
        assert len(torsions) == sum((len(x) for x in smiles.values()))

        cfg = configs.smarts_extender_config(depth, depth, True)
        mapper.mapper_smarts_extend(cfg, torsions)

        config = configs.mapper_config(1, False, "high")
        t_union = mapper.union_list_parallel(torsions, config)
        print("Union:")
        t_union.select = tuple(
            list(t_union.select)
            + [x for x in t_union.nodes if x not in t_union.select]
        )
        print(t_union.select)
        print(gcd.smarts_encode(t_union))


        B = graphs.structure_branch(
            t_union, {k: k for k in t_union.select[:4]}, branches, depth
        )
        B = list(B)
        print(f"There are {len(B)} branches")

        for b in B:
            print(
                len(b.nodes),
                graphs.structure_max_depth(b),
                sorted(b.nodes.keys()),
                gcd.smarts_encode(b),
            )

        if 1:
            scfg = configs.smarts_splitter_config(1, 1, 0, branches, depth, depth)
            spl = splits.split_single_bits(topology.dihedral, scfg, S0, torsions)
            uniq = set()
            for s, _ in spl:
                if s in uniq:
                    print("X", graphs.structure_max_depth(s), len(s.nodes), f"{hash(s):+d}", gcd.smarts_encode(s))
                    continue
                else:
                    uniq.add(s)
                    print("O", graphs.structure_max_depth(s), len(s.nodes), f"{hash(s):+d}", gcd.smarts_encode(s))

smiles = {}
s = "[H:25][c:8]1[c:9]([c:10]([c:11]([c:12]([c:7]1[C:5](=[O:6])[N:4]([H:24])[C:2]([H:20])([C:1]([H:17])([H:18])[H:19])[C:3]([H:21])([H:22])[H:23])[H:28])[H:27])[C:13]([H:29])([H:30])[N:14]([H:31])[N:15]([H:32])[C:16]([H:33])([H:34])[H:35])[H:26]"
smiles[s] = [(11, 10, 9, 26), (9, 10, 11, 27)]
s = "[H:11][c:1]1[c:2]([c:3]([c:4]([c:5]([c:6]1[H:15])[H:14])[C:7]([H:16])([H:17])[C:8]([H:18])([H:19])[N:9]([H:20])[N:10]([H:21])[H:22])[H:13])[H:12]"
smiles[s] = [(5, 4, 3, 13), (3, 4, 5, 14)]

branches = 3
depth = 2
#test(smiles, branches, depth)

smiles = {}
s = '[H:26][c:6]1[c:7]([c:8]2[c:9]([c:10]([c:11]([c:12]([c:13]2[S:14](=[O:15])(=[O:16])[N:17]([H:31])[N:18]([H:32])[H:33])[H:30])[H:29])[H:28])[c:4]([c:5]1[H:25])[N:2]([C:1]([H:19])([H:20])[H:21])[C:3]([H:22])([H:23])[H:24])[H:27]'
smiles[s] = [(6, 5, 4, 9), (8, 13, 12, 11)]
branches = 3
depth = 1
# test(smiles, branches, depth)

smiles = {}
s = '[H:25][c:8]1[c:9]([c:10]([c:11]([c:12]([c:7]1[C:5](=[O:6])[N:4]([H:24])[C:2]([H:20])([C:1]([H:17])([H:18])[H:19])[C:3]([H:21])([H:22])[H:23])[H:28])[H:27])[C:13]([H:29])([H:30])[N:14]([H:31])[N:15]([H:32])[C:16]([H:33])([H:34])[H:35])[H:26]' 
smiles[s] = [(2, 4, 5, 6)]
s = '[H:12][C:1]([H:13])([H:14])[N:2]([C:3]([H:15])([H:16])[C:4]([H:17])([H:18])[C:5](=[O:6])[O-:7])[N+:8]([C:9]([H:19])([H:20])[H:21])([C:10]([H:22])([H:23])[H:24])[C:11]([H:25])([H:26])[H:27]'
smiles[s] = [(3, 4, 5, 6)]
branches = 0
depth = 0
test(smiles, branches, depth)

