
from besmarts.codecs import native
from besmarts.core import mapper
from besmarts.core import graphs
from besmarts.core import configs
from pprint import pprint
from besmarts.core.primitives import primitive_key

from besmarts.codecs.codec_rdkit import graph_codec_rdkit

from besmarts.resolve.enumeration import enumerate_smiles
from besmarts.resolve.rulesets import visitor_ruleset

gcd = graph_codec_rdkit()
gcd.atom_primitives = tuple(list(gcd.atom_primitives) + [primitive_key.VALENCE])
i = 0
j = 20
# g = native.graph_codec_native_load("./big.bes")[0]

smi = "CC(C)(O1)C[C@@H](O)[C@@]1(O2)[C@@H](C)[C@@H]3CC=C4[C@]3(C2)C(=O)C[C@H]5[C@H]4CC[C@@H](C6)[C@]5(C)Cc(n7)c6nc(C[C@@]89(C))c7C[C@@H]8CC[C@@H]%10[C@@H]9C[C@@H](O)[C@@]%11(C)C%10=C[C@H](O%12)[C@]%11(O)[C@H](C)[C@]%12(O%13)[C@H](O)C[C@@]%13(C)CO"
gcd.smiles_config.strip_hydrogen = True

g = gcd.smiles_decode(smi)

bonds = graphs.graph_to_structure_bonds(g)
mapper.mapper_smarts_extend(configs.smarts_extender_config(1, 1, True), bonds)
# bonds = [graphs.structure_remove_unselected(b) for b in bonds]
print(bonds[i].select, bonds[j].select)

def align_list(structs):
    ref = graphs.structure_copy(structs[0])

    for i, s in enumerate(structs[1:]):
        print(i, gcd.smarts_encode(s))

        T = mapper.map_to(ref, s, add_nodes=3, fill=0)
        ref = mapper.union(T.G, T.H, map=T.map)
        print(gcd.smarts_encode(ref))
    for d in range(graphs.structure_max_depth(ref)+1):
        for n in graphs.structure_vertices_at_depth(ref, d):
            print(f"depth={d:d} node={n:d} n_adj={graphs.subgraph_connection(ref, n)}")
    for k in ref.nodes:
        for p, b in ref.nodes[k].primitives.items():
            print(f"{k:3d} {p.value:30s} {str(b):>20}")
    for k, v in ref.nodes.items():
        v.disable(primitive_key.CONNECTIVITY_TOTAL)
        v.disable(primitive_key.CONNECTIVITY_RING)
        v.disable(primitive_key.RING_SMALLEST)
        v.disable(primitive_key.AROMATIC)
        v.primitives[primitive_key.HYDROGEN][:] = True
    for k, v in ref.edges.items():
        v.disable(primitive_key.BOND_RING)
        if v.primitives[primitive_key.BOND_ORDER][5]:
            v.primitives[primitive_key.BOND_ORDER][1] = True
            v.primitives[primitive_key.BOND_ORDER][2] = True
        v.primitives[primitive_key.BOND_ORDER][3:] = False

    ref = graphs.structure_remove_unselected(ref)
    print(gcd.smarts_encode(ref))
    # expand(ref)

def expand(struct):

    mols = []
    smis = set()
    rulesets = [visitor_ruleset(gcd.primitive_codecs)]
    select = [struct.select[i] for i in struct.topology.primary]
    
    for i, P in enumerate(enumerate_smiles(struct, rulesets)):
        if P not in mols:
            mols.append(P)
            smi = gcd.smiles_encode(graphs.structure(P.nodes, P.edges, select, struct.topology))
            if smi not in smis:
                smis.add(smi)
                print(f"{len(smis):4d} {smi}")

def run(A, B, add_nodes=False):
    T = mapper.map_to(bonds[i], bonds[j], add_nodes=add_nodes)
    print(T.map)
    print(bonds[i].select, bonds[j].select)
    print(T.G.select, T.H.select)
    print(gcd.smarts_encode(bonds[i]))
    print(gcd.smarts_encode(bonds[j]))
    for k,v in T.map.items():
        for p, b in T.G.nodes[k].primitives.items():
            print(f"{k:3d} {v:3d} {p.value:30s} {str(b):>20}")
            print(f"{'':38s} {str(T.H.nodes[v].primitives[p]):>20}")
        print("-------------------------------------------------------------------")
    print("G max depth", graphs.structure_max_depth(T.G))
    print("G depths", graphs.structure_node_depths(T.G))
    print("H max depth", graphs.structure_max_depth(T.H))
    print("H depths", graphs.structure_node_depths(T.H))
    print(T.map)
    for d in range(graphs.structure_max_depth(T.G)+1):
        for n in graphs.structure_vertices_at_depth(T.G, d):
            m = T.map.get(n)
            m_adj = None
            if m:
                m_adj = graphs.subgraph_connection(T.H, m)
            print(f"depth={d:d} node={n:d} maps_to={m} n_adj={graphs.subgraph_connection(T.G, n)}, m_adj={m_adj}")
    print("REVERSE")
    M = {v:k for k,v in T.map.items()}
    for d in range(graphs.structure_max_depth(T.H)+1):
        for n in graphs.structure_vertices_at_depth(T.H, d):
            m = M.get(n)
            m_adj = None
            if m:
                m_adj = graphs.subgraph_connection(T.G, m)
            print(f"depth={d:d} node={n:d} maps_to={m} n_adj={graphs.subgraph_connection(T.H, n)}, m_adj={m_adj}")
    print(T.G.select)
    print(T.H.select)

# run(bonds[i], bonds[j], True)
# run(bonds[i], bonds[j], False)
align_list(bonds[:])
