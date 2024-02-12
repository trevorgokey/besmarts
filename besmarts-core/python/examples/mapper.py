
from besmarts.codecs import codec_native as native
from besmarts.core import mapper
from pprint import pprint

g, h = native.graph_codec_native_load("./mapper.bes")

print("g")
pprint(g.nodes)
pprint(g.edges)

print("h")
pprint(h.nodes)
pprint(h.edges)

m = {1:1, 2:2}
G, H, M = mapper.mapper_compose_graphs(g, h, m, add_nodes=True, fill_new_nodes=True)


print("G")
pprint(G.nodes)
pprint(G.edges)

print("H")
pprint(H.nodes)
pprint(H.edges)

print("Map")
pprint(M)

J = mapper.union(G, H, map=M)
print("UNION")
pprint(J.nodes)
pprint(J.edges)


print("T")
T = mapper.structure_mapper(g, add_nodes=True, fill_new_nodes=True)
T.add(g)
T.add(h)
print("Domain:")
pprint(T.domain.nodes)
pprint(T.domain.edges)


print("Codomains")
for i, (g, cdo) in enumerate(T.codomain_map.items(), 1):
    print(i)
    pprint(cdo.nodes)
    pprint(cdo.edges)
    pprint(T.codomain[cdo])
    print("was built from original:")
    pprint(g.nodes)
    pprint(g.edges)


mapper.mapper_save(T, "mapped.bem")

