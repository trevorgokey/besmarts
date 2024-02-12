
from besmarts.codecs.codec_rdkit import graph_codec_rdkit

pcp = graph_codec_rdkit()
pcp.smiles_config.strip_hydrogen = True

smi = "CC(C)(O1)C[C@@H](O)[C@@]1(O2)[C@@H](C)[C@@H]3CC=C4[C@]3(C2)C(=O)C[C@H]5[C@H]4CC[C@@H](C6)[C@]5(C)Cc(n7)c6nc(C[C@@]89(C))c7C[C@@H]8CC[C@@H]%10[C@@H]9C[C@@H](O)[C@@]%11(C)C%10=C[C@H](O%12)[C@]%11(O)[C@H](C)[C@]%12(O%13)[C@H](O)C[C@@]%13(C)CO"

print("Input SMILES:\n" + smi)
beg = pcp.smiles_decode(smi)
print("Encoded SMILES:\n" + pcp.smiles_encode(beg))


