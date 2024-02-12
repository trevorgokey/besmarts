
from besmarts.codecs.codec_rdkit import graph_codec_rdkit
from besmarts.core.graphs import graph_to_structure_bonds, graph_to_structure_angles, graph_to_structure_dihedrals, graph_to_structure_impropers

from besmarts.hierarchy.smirnoff import hindex_smirnoff_load_xml
from besmarts.hierarchy.smirnoff_assign import assign
from besmarts.codecs.assign import assign_rdkit
import pprint

gcd = graph_codec_rdkit()

ff = hindex_smirnoff_load_xml("in.offxml", gcd)


smi = "CC(C)(O1)C[C@@H](O)[C@@]1(O2)[C@@H](C)[C@@H]3CC=C4[C@]3(C2)C(=O)C[C@H]5[C@H]4CC[C@@H](C6)[C@]5(C)Cc(n7)c6nc(C[C@@]89(C))c7C[C@@H]8CC[C@@H]%10[C@@H]9C[C@@H](O)[C@@]%11(C)C%10=C[C@H](O%12)[C@]%11(O)[C@H](C)[C@]%12(O%13)[C@H](O)C[C@@]%13(C)CO"

g = gcd.smiles_decode(smi)

print(gcd.smiles_encode(g))

# bonds = next((x for x in ff.entries.values() if x.key == 'Bonds'))
# x = assign_bonds(ff, bonds, make_rdmol(gcd.smiles_config, smi))
x = assign_rdkit(gcd, ff, smi)
pprint.pprint(x)
exit()
