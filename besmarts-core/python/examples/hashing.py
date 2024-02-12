
from besmarts.core import topology
from besmarts.core import graphs
from besmarts.core import graph_visitors
from besmarts.core import splits
from besmarts.core import configs
from besmarts.codecs import codec_native, codec_rdkit


gcd = codec_rdkit.graph_codec_rdkit()
gcd.atom_primitives = gcd.list_implemented_atom_primitives()

g_lst = codec_native.graph_codec_native_load("./out.bes")
g = g_lst[0]

bonds = graphs.graph_to_structure_bonds(g)

# for bnd in bonds:
#     for bit in graph_visitors.structure_iter_bits(bnd):
#         print(gcd.smarts_encode(bit), hash(bit))

S0 = graphs.structure_copy(bonds[0])
graphs.graph_fill(S0)

splitter = configs.smarts_splitter_config(0,1,0,0,0,0,0,1)

for bit, m in splits.split_single_bits(topology.bond, splitter, S0, bonds):
    print(gcd.smarts_encode(bit), hash(bit))
