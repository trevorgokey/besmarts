

# from besmarts.codecs.codec_rdkit import graph_codec_rdkit
from besmarts.core.codecs import intvec_codec
from besmarts.codecs import codec_native
from timeit import timeit

n = 1 
# gcd = graph_codec_rdkit()

g = codec_native.graph_codec_native_load("./out.bes")[0]
t = timeit('codec_native.graph_codec_native_load("./out.bes")[0]', globals=globals(), number=n)/n

atom_prims = []
for k, v in g.nodes.items():
    atom_prims.extend(list(v.primitives))
    break

bond_prims = []
for k, v in g.edges.items():
    bond_prims.extend(list(v.primitives))
    break

codecs = codec_native.primitive_codecs_get()
gcd = codec_native.graph_codec_native(codecs, atom_prims, bond_prims)

# gcd.atom_primitives = atom_prims
# gcd.bond_primitives = bond_prims


giv = intvec_codec(gcd.primitive_codecs, atom_prims, bond_prims)
print(f"SMILES loaded, time={t} sec/call")


smi_a = gcd.smiles_encode(g)
t = timeit('gcd.smiles_encode(g)', globals=globals(), number=n)/n
print(smi_a)
print(f"SMILES encoded, time={t} sec/call")

# t = timeit('gcd.smiles_decode(smi_a)', globals=globals(), number=n)/n
# print(f"SMILES decoded, time={t} sec/call")

intvec = giv.graph_encode(g)
t = timeit('giv.graph_encode(g)', globals=globals(), number=n)/n
# print(intvec.v)
print(f"intvec Encoded, time={t} sec/call")

g = giv.graph_decode(intvec)
t = timeit('giv.graph_decode(intvec)', globals=globals(), number=n)/n
print(g)
print(f"intvec Decoded, time={t} sec/call")

smi_b = gcd.smiles_encode(g)
t = timeit('gcd.smiles_encode(g)', globals=globals(), number=n)/n
print(smi_b)
print(f"SMILES encoded, time={t} sec/call")

assert smi_a == smi_b
print("SMILES are the same; test passed!")
