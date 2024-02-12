
from besmarts.codecs.codec_rdkit import graph_codec_rdkit
from besmarts.core import graphs

gcd = graph_codec_rdkit()
# g = gcd.smiles_decode("[C:1][C:2][C:3][C:4][C:5][C:6]")
g = gcd.smiles_decode("CCCCCC")
print(gcd.smiles_encode(g))
g = graphs.graph(g.nodes, g.edges)
print(gcd.smiles_encode(g))
structures = set(graphs.graph_to_structure_dihedrals(g))
# structures are 1-based indexing
dihedrals = [tuple((i - 1 for i in x.select)) for x in structures]
unique = set()
for s in structures:
    inner = tuple(sorted(s.select[1:3]))
    if inner not in unique:
        print(s.select, gcd.smarts_encode(s))
        unique.add(inner)
