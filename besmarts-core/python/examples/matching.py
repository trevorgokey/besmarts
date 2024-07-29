
"""
examples/matching.py
"""

from besmarts.core import primitives
from besmarts.core import graphs
from besmarts.core import mapper
from besmarts.codecs.codec_rdkit import graph_codec_rdkit

gcd = graph_codec_rdkit(
    (primitives.primitive_key.ELEMENT, primitives.primitive_key.HYDROGEN),
    (primitives.primitive_key.BOND_ORDER,)
)
gcd.smiles_config.strip_hydrogen = False

target = "[#6H3:1](-[#1])(-[#1])(-[#1])-[#6H0:2](-[#6H2])(-[#8H0])-[#6H3]"
queries = [
    "[#6H2:1](~[!H3])~[#6:2]",
    "[#6:1](~[!H3])~[#6:2]",
    "[#6:1]~[#6:2]~[#7]"
]

t = graphs.subgraph_to_structure_bond(gcd.smarts_decode(target))

for query in queries:
    S0 = graphs.subgraph_to_structure_bond(gcd.smarts_decode(query))

    print("testing if", gcd.smarts_encode(t))
    print("is a subset (match) to parameter", gcd.smarts_encode(S0))
    print(mapper.mapper_match(t, S0))
    print("\n")
