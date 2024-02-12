

from besmarts.hierarchy import hierarchy_index
from besmarts.hierarchy import hindex_iters
from besmarts.hierarchy import split
from besmarts.core import mapper


from besmarts.codecs.codec_rdkit import graph_codec_rdkit
from besmarts.codecs import native
from besmarts.core import mapper
from besmarts.core import graphs
from besmarts.core import configs
from pprint import pprint

gcd = graph_codec_rdkit()

g = native.graph_codec_native_load("./big.bes")[0]

atoms = graphs.graph_to_structure_atoms(g)
N = len(atoms)
for i in range(1,2):
    atoms = graphs.graph_to_structure_atoms(g)
    mapper.mapper_smarts_extend(configs.smarts_extender_config(i,i, True), atoms)
# bonds = [graphs.structure_remove_unselected(b) for b in bonds]
    print(f"There are {len(set(atoms))}/{N} unique structures at depth", i)

print("T")

for depth in range(1,3):
    atoms = graphs.graph_to_structure_atoms(g)
    mapper.mapper_smarts_extend(configs.smarts_extender_config(depth, depth, True), atoms)
    T = mapper.structure_mapper(atoms[0], add_nodes=True, fill_new_nodes=False)
    for i, b in enumerate(atoms[1:],1):
        T.add(b)
        print("Added", i, len(atoms)-1, "domain:", len(T.domain.nodes), [len(v.nodes) for v in T.codomain])
    print(f"There are {len(set(atoms))}/{N} unique structures")

# breakpoint()
# for i, b in enumerate(atoms,1):
#     T.add(b)
#     print("Added", i, len(atoms)-1, "domain:", len(T.domain.nodes), [len(v.nodes) for v in T.codomain])
# mapper.mapper_smarts_extend(configs.smarts_extender_config(3, 3, True), atoms)
# for i, b in enumerate(atoms,1):
#     T.add(b)
#     print("Added", i, len(atoms)-1, "domain:", len(T.domain.nodes), [len(v.nodes) for v in T.codomain])
# print("Domain:")
# pprint(T.domain.nodes)
# pprint(T.domain.edges)
assignments = {i: str(i)  for i, _ in enumerate(T.codomain_map)}
for i, a in enumerate(T.codomain_map):
    print(i, gcd.smarts_encode(a))

hidx, match, groups = split.hierarchy_split(T, assignments)
print("match")
print(match)
print("groups")
print(groups)

for e in hindex_iters.hindex_iter_dive(hidx, hidx.root):
    s = " " * hierarchy_index.hindex_node_depth(hidx, e)
    # print(s, e.index, e.key, gcd.smarts_encode(hidx.db[e.key]['graph']))
    print(s, e.index, e.key)
