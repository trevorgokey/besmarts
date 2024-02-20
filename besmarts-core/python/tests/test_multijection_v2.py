
from besmarts.codecs import native
from besmarts.core import mapper
from besmarts.core import graphs
from besmarts.core import configs
from pprint import pprint

g = native.graph_codec_native_load("./big.bes")[0]

bonds = graphs.graph_to_structure_bonds(g)
mapper.mapper_smarts_extend(configs.smarts_extender_config(1, 1, True), bonds)
# bonds = [graphs.structure_remove_unselected(b) for b in bonds]

T = mapper.structure_mapper(bonds[0], add_nodes=3, fill_new_nodes=0)
print("T nodes", len(T.domain.select))
for i, b in enumerate(bonds,1):
    T.add(b)
    print("Added", i, len(bonds), "domain nodes:", T.domain.select, [len(v.select) for v in T.codomain])
# mapper.mapper_smarts_extend(configs.smarts_extender_config(1, 1, True), bonds)
# for i, b in enumerate(bonds,1):
#     T.add(b)
#     print("Added", i, len(bonds), "domain nodes:", T.domain.select, [len(v.select) for v in T.codomain])
# mapper.mapper_smarts_extend(configs.smarts_extender_config(2, 2, True), bonds)
# for i, b in enumerate(bonds,1):
#     T.add(b)
#     print("Added", i, len(bonds), "domain nodes:", len(T.domain.select), [len(v.select) for v in T.codomain])
print("Domain:")
pprint(T.domain.nodes)
pprint(T.domain.edges)
pprint(T.domain.select)


print("Codomains")
for i, (g, cdo) in enumerate(T.codomain_map.items(), 1):

    print(i, "################################")
    # pprint(cdo.nodes)
    # pprint(cdo.edges)
    pprint((T.domain.select, g.select, T.codomain.get(cdo)))
    # print("was built from original:")
    # pprint(g.nodes)
    # pprint(g.edges)

mapper.mapper_save(T, "mapped_lv2.bem")
T2 = mapper.mapper_load("mapped_lv2.bem")
for i, (g, cdo) in enumerate(T2.codomain_map.items(), 1):

    print(i, "################################")
    # pprint(cdo.nodes)
    # pprint(cdo.edges)
    pprint(T.codomain.get(cdo))
mapper.mapper_save(T2, "mapped_lv2.2.bem")
