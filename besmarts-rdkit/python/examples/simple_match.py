
from besmarts.core import topology 
from besmarts.core import primitives
from besmarts.core import graphs 
from besmarts.core import mapper
from besmarts.core import configs
from besmarts.codecs.codec_rdkit import graph_codec_rdkit
from besmarts.assign.hierarchy_assign_rdkit import smarts_hierarchy_assignment_rdkit
import pickle
import pprint

configs.processors = 1
gcd = graph_codec_rdkit(
    # (primitives.primitive_key.ELEMENT, primitives.primitive_key.HYDROGEN), (primitives.primitive_key.BOND_ORDER,)
)
gcd.smiles_config.strip_hydrogen = False
labeler = smarts_hierarchy_assignment_rdkit()

target = r"[c:1]1([H:12])[c:2]([H:13])[c:3]([H:14])[c:4]([S:7](=[O:8])(=[O:9])[N:10]([N:11]([H:18])[H:19])[H:17])[c:5]([H:15])[c:6]1[H:16]"


# target = "[#6;H3:1](-[#1;H0])(-[#1;H0])(-[#1;H0])-[#6;H0:2](-[#6;H2])(-[#8;H0])-[#6;H3]"
#queries = [r'[#1:1]~[#6H1:2](~[H0])~[H1:3](~[*])~[!#1!#16H1:4](~[*])~[*]']
_, cst, _ = pickle.load(open('chk.cst.p', 'rb'))

cst.hierarchy.smarts[30] = '[#1:1]~[#6H1:2](~[H0])~[H1:3](~[*])~[!#1!#16H1:4](~[*])~[*]'
breakpoint()

pprint.pprint(cst.hierarchy.smarts)
labels = labeler.assign(cst.hierarchy, gcd, [target], topology.dihedral)
print(labels.assignments[0].selections)

# t = gcd.smiles_decode(target)
# # t = graphs.structure(t.nodes, t.edges, tuple(sorted(t.nodes)), topology.bond_topology())

# for query in queries:
#     S0 = gcd.smarts_decode(query)
#     S0 = graphs.structure(S0.nodes, S0.edges, tuple(sorted(S0.nodes)), topology.bond_topology())


#     print("testing if", gcd.smarts_encode(t))
#     print("is a subset (match) to parameter", gcd.smarts_encode(S0))
#     print(mapper.mapper_match(t, S0))
#     print("\n")


