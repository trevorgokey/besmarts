
"""
besmarts.examples.test_transcode
"""

from besmarts.core import transcoders
from besmarts.core import graphs
from besmarts.core import topology
from besmarts.core import primitives
from besmarts.core import arrays
from besmarts.core import codecs
from besmarts.codecs import codec_rdkit

def water_to_bead():
    gcd = codec_rdkit.graph_codec_rdkit()

    gcd.primitive_codecs["unit"] = codecs.primitive_codec_unit()
    gcd.primitive_codecs["unit_index"] = codecs.primitive_codec_unit_index()
    gcd.primitive_codecs["link_src"] = codecs.primitive_codec_link_src()
    gcd.primitive_codecs["link_dst"] = codecs.primitive_codec_link_dst()

    gcd.atom_primitives = tuple(["unit", "unit_index"] + list(gcd.atom_primitives))
    gcd.bond_primitives = tuple(["link_src", "link_dst"] + list(gcd.bond_primitives))

    g = gcd.smiles_decode("[OH2:2]")
    g = graphs.graph_to_structure_angles(g)[0]

    # for n in g.nodes:
    #     g.nodes[n].primitives["unit"] = arrays.bitvec()

    #     g.nodes[n].primitives["unit_index"] = arrays.bitvec()
    #     g.nodes[n].select = tuple(["unit", "unit_index"] + list(g.nodes[n].select))

    # for n in g.edges:
    #     g.edges[n].primitives["link_src"] = arrays.bitvec()
    #     g.edges[n].primitives["link_dst"] = arrays.bitvec()
    #     g.edges[n].select = tuple(["link_src", "link_dst"] + list(g.edges[n].select))

    s = graphs.structure({1: g.nodes[1].copy()}, {}, (1,), topology.atom)
    s.nodes[1].clear()

    # for n in s.nodes:
    #     s.nodes[n].primitives["unit"] = arrays.bitvec()

    #     s.nodes[n].primitives["unit_index"] = arrays.bitvec()
    #     s.nodes[n].select = tuple(["unit", "unit_index"] + list(s.nodes[n].select))

    # for n in s.edges:
    #     s.edges[n].primitives["link_src"] = arrays.bitvec()
    #     s.edges[n].primitives["link_dst"] = arrays.bitvec()
    #     s.edges[n].select = tuple(["link_src", "link_dst"] + list(s.edges[n].select))

    s.nodes[1].primitives["unit"][1] = True
    ret = transcoders.transcode(topology.angle_to_atom, g, (1,2,3), s)

    print(gcd.smarts_encode(ret.H))

def cc_to_bead():
    gcd = codec_rdkit.graph_codec_rdkit()

    gcd.primitive_codecs["unit"] = codecs.primitive_codec_unit()
    gcd.primitive_codecs["unit_index"] = codecs.primitive_codec_unit_index()
    gcd.primitive_codecs["link_src"] = codecs.primitive_codec_link_src()
    gcd.primitive_codecs["link_dst"] = codecs.primitive_codec_link_dst()

    gcd.atom_primitives = tuple(["unit", "unit_index"] + list(gcd.atom_primitives))
    gcd.bond_primitives = tuple(["link_src", "link_dst"] + list(gcd.bond_primitives))

    g = gcd.smiles_decode("[H][C:2]([H])([H])[C:3]([H])([H])[H]")
    # g = graphs.graph_to_structure_angles(g)[0]

    # for n in g.nodes:
    #     g.nodes[n].primitives["unit"] = arrays.bitvec()

    #     g.nodes[n].primitives["unit_index"] = arrays.bitvec()
    #     g.nodes[n].select = tuple(["unit", "unit_index"] + list(g.nodes[n].select))

    # for n in g.edges:
    #     g.edges[n].primitives["link_src"] = arrays.bitvec()
    #     g.edges[n].primitives["link_dst"] = arrays.bitvec()
    #     g.edges[n].select = tuple(["link_src", "link_dst"] + list(g.edges[n].select))

    s = graphs.structure({1: g.nodes[1].copy()}, {}, (1,), topology.atom)
    s.nodes[1].clear()

    # for n in s.nodes:
    #     s.nodes[n].primitives["unit"] = arrays.bitvec()

    #     s.nodes[n].primitives["unit_index"] = arrays.bitvec()
    #     s.nodes[n].select = tuple(["unit", "unit_index"] + list(s.nodes[n].select))

    # for n in s.edges:
    #     s.edges[n].primitives["link_src"] = arrays.bitvec()
    #     s.edges[n].primitives["link_dst"] = arrays.bitvec()
    #     s.edges[n].select = tuple(["link_src", "link_dst"] + list(s.edges[n].select))

    s.nodes[1].primitives["unit"][2] = True
    ret = transcoders.transcode(topology.bond_to_atom, g, (2,3), s)

    assns = transcoders.transcode_assignment_positions(
        topology.bond_to_atom, assn, (2,3)
    )

    print(gcd.smarts_encode(graphs.graph_to_subgraph(ret.H, ret.H.nodes)))

if __name__ == "__main__":

    water_to_bead()
