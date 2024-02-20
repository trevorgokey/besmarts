"""
Example using RDKit to decode a SMILES string and save to a native BESMARTS format.
"""

from besmarts.codecs import codec_rdkit 
from besmarts.codecs import codec_native
from besmarts.core import graphs

gcd = codec_rdkit.graph_codec_rdkit()
smi = "CCC"

G = gcd.smiles_decode(smi)
codec_native.graph_codec_native_save("propane.bg", [G])


bonds = graphs.graph_to_structure_bonds(G)
codec_native.graph_codec_native_save("propane_bonds.bg", bonds)


G = codec_native.graph_codec_native_load("propane.bg")[0]

# unmapped
print(gcd.smiles_encode(G))

# mapped
print(gcd.smiles_encode(graphs.graph_to_subgraph(G, tuple(G.nodes))))
