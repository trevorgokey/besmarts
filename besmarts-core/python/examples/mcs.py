from besmarts.core import graphs
from besmarts.core import mapper
from besmarts.core import configs
from besmarts.codecs.codec_rdkit import graph_codec_rdkit
from besmarts.core import primitives

# from besmarts.resolve import enumeration
# from besmarts.resolve import visitor_rulesets


def MCS(g, h):

    g_atoms = graphs.graph_to_structure_atoms(g)
    h_atoms = graphs.graph_to_structure_atoms(h)

    intersect = configs.mapper_config(False, False, "high")

    mcs = {}
    pruned = {gi: [] for gi, _ in enumerate(g_atoms)}

    added = True
    i = -1
    while added:
        added = False
        i += 1
        config = configs.smarts_extender_config(i, i, False)
        modified = mapper.mapper_smarts_extend(config, g_atoms)
        if not modified and i > 0:
            break
        modified = mapper.mapper_smarts_extend(config, h_atoms)
        if not modified and i > 0:
            break
        for gi, ga in enumerate(g_atoms):
            for hi, ha in enumerate(h_atoms):
                if hi in pruned[gi]:
                    continue

                T = mapper.map_to(ga, ha, intersect)
                if not T.map:
                    continue
                r = mapper.intersection(T.G, T.H, map=T.map)
                r = graphs.structure_prune_null_nodes(r)

                if r is None:
                    pruned[gi].append(hi)
                    continue

                n = len(r.nodes)
                if n not in mcs:
                    mcs[n] = []
                mcs[n].append((T.G, T.H, r, T.map))
                added = True

    return mcs


gcd = graph_codec_rdkit()
# gcd.atom_primitives = tuple(x for x in gcd.atom_primitives if x != primitives.primitive_key.HYDROGEN)
gcd.smiles_config.strip_hydrogen = True

g = gcd.smiles_decode(
    "C[C@@H](C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O)N"
)  # trialanine
# h = gcd.smiles_decode("C(C(=O)NCC(=O)NCC(=O)O)N") #triglycine
# h = gcd.smiles_decode("CC=O") #triglycine
h = gcd.smiles_decode("NCC=O")  # triglycine

graphs.graph_disable_primitives_atom(g, [primitives.primitive_key.HYDROGEN])
graphs.graph_disable_primitives_atom(h, [primitives.primitive_key.HYDROGEN])

mcs = MCS(g, h)


if mcs:
    for (gi, hi, r, M) in mcs[max(mcs)]:
        graphs.graph_set_primitives_atom(r, gcd.atom_primitives)
        print(gcd.smarts_encode(r))
        sg = graphs.subgraph(r.nodes, r.edges, tuple(M))
        graphs.graph_set_primitives_atom(sg, gcd.atom_primitives)
        print(gcd.smiles_encode(sg))
        break
