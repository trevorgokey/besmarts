"""
examples/split_smarts.py
"""

from besmarts.core import graphs
from besmarts.core import configs
from besmarts.core import topology
from besmarts.core import splits
from besmarts.core import codecs
from besmarts.core import configs
from besmarts.core import compute

from besmarts.core.primitives import primitive_key

# use the RDKit plugin
from besmarts.codecs import codec_rdkit

configs.processors = 16

prims = (primitive_key.ELEMENT, primitive_key.HYDROGEN), (
    primitive_key.BOND_ORDER,
)

# use all available primitives
# prims = (None, None)

pcp = codec_rdkit.graph_codec_rdkit(*prims)
pcp.smiles_config.strip_hydrogen = False
icd = codecs.intvec_codec(
    pcp.primitive_codecs,
    pcp.atom_primitives,
    pcp.bond_primitives
)

# smi = "CC(C)(O1)C[C@@H](O)[C@@]1(O2)[C@@H](C)[C@@H]3CC=C4[C@]3(C2)C(=O)C[C@H]5[C@H]4CC[C@@H](C6)[C@]5(C)Cc(n7)c6nc(C[C@@]89(C))c7C[C@@H]8CC[C@@H]%10[C@@H]9C[C@@H](O)[C@@]%11(C)C%10=C[C@H](O%12)[C@]%11(O)[C@H](C)[C@]%12(O%13)[C@H](O)C[C@@]%13(C)CO"
# smi = "CCCO"
# smi = "CCCC"
# trialanine
smi = "C[C@@H](C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O)N"
N = 1000

# G = pcp.smiles_decode(smi)
G = {0: pcp.smiles_decode(smi)}
ic_list = [s for s in graphs.graph_to_structure_torsions(G[0])]
selections = [(i, x) for i in G for x in graphs.graph_torsions(G[i])]
G[0] = icd.graph_encode(G[0])

# set these all to 1 to split on neighbors too
branch_min = 0
branch_limit = 0
branch_depth_min = 0
branch_depth = 0

bit_depth_min = 1
bit_depth_max = 1

splitter = configs.smarts_splitter_config(
    bit_depth_min,
    bit_depth_max,
    branch_min,
    branch_limit,
    branch_depth_min,
    branch_depth,
    unique=False,
    return_matches=True,
    max_splits=0,
    split_general=True,
    split_specific=True,
    unique_compliments=False,
    unique_compliments_prefer_min=True,
)

# for this to work, we need to extend our graphs to at least the depth of S0
extender = configs.smarts_extender_config(branch_depth, branch_depth, True)
graphs.structure_extend(extender, ic_list)

S0 = pcp.smarts_decode("[*:1]~[*:2]~[*:3]~[*:4]")
S0 = graphs.structure(S0.nodes, S0.edges, (1, 2, 3, 4,), topology.torsion)

for i, f in enumerate(ic_list):
    print(i, pcp.smarts_encode(f))

configs.remote_compute_enable = False

wq = compute.workqueue_local("127.0.0.1", 63210)
results: splits.split_return_type = splits.split_structures_distributed(splitter, S0, G, selections, wq, icd)


seen = {}
keep = {}

print("Results:", len(results.splits))
for j, (Sj, matches, bj) in enumerate(
    zip(results.splits, results.matched_idx, results.shards), 1
):
    Sj = graphs.structure(Sj.nodes, Sj.edges, Sj.select, results.topology)
    atoms, bits = len(Sj.select), graphs.graph_bits(Sj, maxbits=True)
    matches = tuple(sorted(matches))
    unmatches = tuple(
        sorted([i for i in range(len(ic_list)) if i not in matches])
    )
    entry = tuple(sorted([matches, unmatches]))
    if len(matches) > 0 and len(ic_list) != len(matches):
        if entry in seen:
            if (atoms, bits) < seen[entry]:
                seen[entry] = (atoms, bits)
                keep[entry] = j
        else:
            seen[entry] = (atoms, bits)
            keep[entry] = j

unique = {}
found = 0
for j, (Sj, matches, bj) in enumerate(
    zip(results.splits, results.matched_idx, results.shards), 1
):

    matches = tuple(matches)
    l = unique.get(matches, list())
    l.append((Sj, bj))
    unique[matches] = l

for j, (matches, params) in enumerate(unique.items(), 1):
    matches = tuple(matches)
    found += 1
    if splitter.return_matches:
        print(
            f"{found:4d}",
            f"{j:4d}",
            "match:",
            f"{len(matches):4d}",
            "unmatched:",
            f"{len(ic_list) - len(matches):4d}",
        )
    else:
        print(
            f"{found:4d}",
            f"{j:4d}",
        )
    for k, (Sj, bj) in enumerate(params, 1):

        Sj = graphs.structure(Sj.nodes, Sj.edges, Sj.select, results.topology)
        print(f"   {k:2d} Sj:", pcp.smarts_encode(Sj), hash(Sj))
        # print(f"   {k:2d} Sj:    ", Sj.nodes)
    if splitter.return_matches:
        print("      ", matches)
        for i, f in enumerate(ic_list):
            if i in matches:
                print(f"{i:4d}", " -> ", pcp.smarts_encode(f), f.select, hash(f))
            else:
                print(f"{i:4d}", pcp.smarts_encode(f), f.select, hash(f))
    print("####################################")
