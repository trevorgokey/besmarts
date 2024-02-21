
from besmarts.core import topology 
from besmarts.core import primitives
from besmarts.core import graphs 
from besmarts.core import mapper
from besmarts.codecs.codec_rdkit import graph_codec_rdkit

gcd = graph_codec_rdkit(
    (primitives.primitive_key.ELEMENT, primitives.primitive_key.HYDROGEN), (primitives.primitive_key.BOND_ORDER,)
)
gcd.smiles_config.strip_hydrogen = False

target = "[#6;H3:1](-[#1;H0])(-[#1;H0])(-[#1;H0])-[#6;H0:2](-[#6;H2])(-[#8;H0])-[#6;H3]"
queries = ["[#6H2:1](~[!H3])~[#6:2]", "[#6:1](~[!H3])~[#6:2]", "[#6:1]~[#6:2]~[#7]"]

t = gcd.smarts_decode(target)
t = graphs.structure(t.nodes, t.edges, tuple(sorted(t.nodes)), topology.bond_topology())

for query in queries:

    S0 = gcd.smarts_decode(query)

    S0 = graphs.structure(
        S0.nodes,
        S0.edges,
        tuple(sorted(S0.nodes)),
        topology.bond_topology()
    )

    print("testing if", gcd.smarts_encode(t))
    print("is a subset (match) to parameter", gcd.smarts_encode(S0))
    print(mapper.mapper_match(t, S0))
    print("\n")


