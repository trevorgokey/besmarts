
from besmarts.codecs.codec_rdkit import graph_codec_rdkit

from besmarts.core import codecs, graphs, topology, mapper, arrays

besmarts_codecs = codecs.primitive_codecs_get()
supersmarts_atom = ["unit", "unit_index"]
supersmarts_bond = ["link_src", "link_dst"]

# create the supersmartser

gcd = graph_codec_rdkit()
gcd.smiles_config.strip_hydrogen = True

bb_sma = "[#7:1]-[#6X4:2]-[#6X3:3](=[#8X1])-[#8:4]"
bb_g = gcd.smarts_decode(bb_sma)
bb_g = graphs.graph_invert_all(bb_g)
bb_sg = graphs.graph_to_subgraph(bb_g, (1,2,3,4,5))
# print("BB UNIT SMARTS :\n" + gcd.smarts_encode(bb_sg))

bb_unit = graphs.subgraph_to_structure(bb_sg, topology.dihedral)

print("BB UNIT SMARTS :\n" + gcd.smarts_encode(bb_unit))

smi = "NCC(=O)O" # glycine
g = gcd.smiles_decode(smi)

print("Molecule SMARTS\n", gcd.smarts_encode(graphs.graph_to_subgraph(g, g.nodes)))

for n in g.nodes:
    g.nodes[n].primitives["unit"] = arrays.bitvec()

    g.nodes[n].primitives["unit_index"] = arrays.bitvec()
    g.nodes[n].select = tuple(["unit", "unit_index"] + list(g.nodes[n].select))

for n in g.edges:
    g.edges[n].primitives["link_src"] = arrays.bitvec()
    g.edges[n].primitives["link_dst"] = arrays.bitvec()
    g.edges[n].select = tuple(["link_src", "link_dst"] + list(g.edges[n].select))

gcd.primitive_codecs["unit"] = codecs.primitive_codec_unit()
gcd.primitive_codecs["unit_index"] = codecs.primitive_codec_unit_index()
gcd.primitive_codecs["link_src"] = codecs.primitive_codec_link_src()
gcd.primitive_codecs["link_dst"] = codecs.primitive_codec_link_dst()

gcd.atom_primitives = tuple(["unit", "unit_index"] + list(gcd.atom_primitives))
gcd.bond_primitives = tuple(["link_src", "link_dst"] + list(gcd.bond_primitives))

dihedrals = graphs.graph_to_structure_dihedrals(g)

for d in dihedrals:
    d.select = tuple(list(d.select) + [x for x in sorted(d.nodes) if x not in d.select])
    print("Query", gcd.smarts_encode(bb_unit))
    print("Target", gcd.smarts_encode(d))
    print("Target select", d.select)

    T = mapper.mapper_force_subset(bb_unit, d)
    if T is not None:
        print(T.map)
        for k,v in T.map.items():
            g.nodes[v].primitives["unit"][1] = True
            g.nodes[v].primitives["unit_index"][k] = True

            for v2 in graphs.graph_connection(g, v):
                e = graphs.edge((v,v2))
                g.edges[e].primitives["link_src"][v2] = True
                g.edges[e].primitives["link_dst"][v2] = True
                g.edges[e].primitives["link_src"][v] = True
                g.edges[e].primitives["link_dst"][v] = True

    else:
        print(T)
    
print("SuperSMARTS of molecule")
graphs.graph_set_primitives_atom(g, ["unit", "unit_index"])
graphs.graph_disable_primitives_bond(g, ["link_src", "link_dst"])

print(gcd.smarts_encode(g))

# Now the graph has unit information. Now we might want to compact it into a
# single node. We then set all frontier edges to have the correct link_src set,
# and then I guess set all link_dst?

# I can then encode this as a SMARTS which will be a very compact SMARTS, i.e.
# just [$1:1]. Ideally this can then be used as a coarse representation of a large SMARTS






exit()
# for key in supersmarts_atom:
#     pcp.primitive_codecs[key] = besmarts_codecs[key]

# pcp.atom_primitives = tuple(supersmarts_atom + list(pcp.atom_primitives))
# pcp.bond_primitives = tuple(supersmarts_bond + list(pcp.bond_primitives))


smi = "CC(C)(O1)C[C@@H](O)[C@@]1(O2)[C@@H](C)[C@@H]3CC=C4[C@]3(C2)C(=O)C[C@H]5[C@H]4CC[C@@H](C6)[C@]5(C)Cc(n7)c6nc(C[C@@]89(C))c7C[C@@H]8CC[C@@H]%10[C@@H]9C[C@@H](O)[C@@]%11(C)C%10=C[C@H](O%12)[C@]%11(O)[C@H](C)[C@]%12(O%13)[C@H](O)C[C@@]%13(C)CO"
# from besmarts! It gets rid of chirality
smi = "CC(C)(CC1[OH])OC1(OC2)C(C)C(C32C(=O)CC4C5(C)Cc6nc7C8)CC=C3C4CCC5Cc6nc7CC(C)(C9CC([OH])C%10(C)C%11([OH])C%12C)C8CCC9C%10=CC%11OC%12(C(C%13)[OH])OC%13(C)C[OH]"

print("Input SMILES:\n" + smi)
beg = gcd.smiles_decode(smi)
print("Encoded SMILES:\n" + gcd.smiles_encode(beg))
