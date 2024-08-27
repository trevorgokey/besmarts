"""
examples/01_smiles_io.py

Example using RDKit to decode a SMILES string and save to a native BESMARTS
format.
"""

from besmarts.codecs import codec_rdkit
from besmarts.codecs import codec_native
from besmarts.core import graphs

gcd = codec_rdkit.graph_codec_rdkit()
smi = "CCC"

G = gcd.smiles_decode(smi)

# unmapped
print(gcd.smiles_encode(G))

print(f"Input SMILES:\n{smi}")
# mapped
print("Output SMILES:")
print(gcd.smiles_encode(graphs.graph_as_subgraph(G, tuple(G.nodes))))

# Save to disk instead
# codec_native.graph_codec_native_save("propane.bg", [G])
out = codec_native.graph_save(G)

print("Serialized SMILES:")
for line in out:
    print(line)

bonds = graphs.graph_to_structure_bonds(G)
# codec_native.graph_codec_native_save("propane_bonds.bg", bonds)

# Load from disk
# G = codec_native.graph_codec_native_load("propane.bg")[0]
G = codec_native.graph_load(out)
