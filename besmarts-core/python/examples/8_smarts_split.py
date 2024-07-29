"""
examples/split_smarts.py

Enumerate all SMARTS patterns that can split (partition) the bonds of CCO into
two groups.
"""

from besmarts.core import graphs
from besmarts.core import configs
from besmarts.core import splits
from besmarts.core import codecs
from besmarts.core import compute
from besmarts.core.primitives import primitive_key

# Use the RDKit plugin for SMILES perception and substructure search
from besmarts.codecs import codec_rdkit

# Set the compute environment. Use 4 local cores
configs.processors = 4
configs.remote_compute_enable = False
wq = compute.workqueue_local("127.0.0.1", 63210)

prims = (
    (primitive_key.ELEMENT, primitive_key.HYDROGEN),
    (primitive_key.BOND_ORDER,)
)

# Use all available primitives
# prims = (None, None)

gcd = codec_rdkit.graph_codec_rdkit(*prims)

# Packs graphs into a space efficient serialization for sending over networks
# for distributed computing
icd = codecs.intvec_codec(
    gcd.primitive_codecs,
    gcd.atom_primitives,
    gcd.bond_primitives
)

# Currently the splitting functions are designed to handle intvec encoded
# graphs since we assume jobs will be large and distributed
smi = "CCO"
G = {0: gcd.smiles_decode(smi)}
ic_list = [s for s in graphs.graph_to_structure_bonds(G[0])]
selections = [(i, x) for i in G for x in graphs.graph_bonds(G[i])]
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
    unique_complements=False,
    unique_complements_prefer_min=True,
)

# for this to work, we need to extend our graphs to at least the depth of S0
extender = configs.smarts_extender_config(branch_depth, branch_depth, True)
graphs.structure_extend(extender, ic_list)

S0 = gcd.smarts_decode("[*:1]~[*:2]")
S0 = graphs.subgraph_to_structure_bond(S0)

for i, f in enumerate(ic_list):
    print(i, gcd.smarts_encode(f))


# The function that does the actual work
results: splits.split_return_type = splits.split_structures_distributed(
    splitter, S0, G, selections, wq, icd)

# Custom processing of results as an example to use the results

seen = {}
keep = {}

print("Results:", len(results.splits))
for j, (Sj, matches, bj) in enumerate(
    zip(results.splits, results.matched_idx, results.shards), 1
):
    Sj = graphs.subgraph_as_structure(Sj, results.topology)
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

        Sj = graphs.subgraph_as_structure(Sj, results.topology)
        print(f"   {k:2d} Sj:", gcd.smarts_encode(Sj))
        # print(f"   {k:2d} Sj:    ", Sj.nodes)
    if splitter.return_matches:
        print("      ", matches)
        for i, f in enumerate(ic_list):
            if i in matches:
                print(f"{i:4d}", " -> ", f.select, gcd.smarts_encode(f))
            else:
                print(f"{i:4d}", f.select, gcd.smarts_encode(f))
    print("####################################")
